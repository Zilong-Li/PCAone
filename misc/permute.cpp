#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <vector>
#include <iterator>

using namespace std;

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;

struct Line
{
    std::string data;
    operator std::string const&() const {return data;}
    friend std::istream& operator>>(std::istream& is, Line& line) {
        return std::getline(is, line.data);
    }
};

size_t count_lines(const string& fpath)
{
    std::ifstream in(fpath);
    size_t count = 0;
    string line;
    while (getline(in, line)) {
        count++;
    }
    return count;
}

void permute_plink2(const string& fin, const string& fout,uint blocksize = 100)
{
    uint64 nsnps = count_lines(fin + ".bim");
    uint64 nsamples = count_lines(fin + ".fam");
    uint64 bed_bytes_per_snp = (nsamples+3)>>2;
    uint64 bed_bytes_per_block = bed_bytes_per_snp * blocksize;

    ios_base::sync_with_stdio(false);
    std::ifstream in(fin + ".bed", std::ios::binary);
    std::ofstream out(fout + ".bed", std::ios::binary);
    if (!in.is_open()) {
        throw std::invalid_argument("ERROR: Cannot open bed file.\n");
    }
    uchar header[3];
    in.read(reinterpret_cast<char *> (&header[0]), 3);
    if ( (header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01) ) {
        throw std::invalid_argument("ERROR: Incorrect magic number in bed file.\n");
    }
    out.write(reinterpret_cast<char *> (&header[0]), 3);
    std::ifstream in_bim(fin + ".bim", std::ios::in);
    std::ofstream out_bim(fout + ".bim", std::ios::out);
    vector<string> bims(std::istream_iterator<Line>{in_bim},
                        std::istream_iterator<Line>{});
    vector<uint64> seqs;
    seqs.reserve(nsnps);
    for(uint64 i = 0; i < nsnps; i++) { seqs.push_back(i); }
    auto rng = std::default_random_engine {};
    std::shuffle(std::begin(seqs), std::end(seqs), rng);

    vector<uchar> inbed;
    inbed.resize(bed_bytes_per_block);
    uint64 nblocks = (unsigned int)ceil((double)nsnps / blocksize);
    uint64 c, i, j, idx;
    for(i = 0; i < nblocks; i++)
    {
        in.read(reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_block);
        for(j = 0; j < blocksize && c < nsnps; j++)
        {
            c = i * blocksize + j;
            idx = 3 + seqs[c] * bed_bytes_per_snp;
            out.seekp(idx, std::ios_base::beg);
            out.write(reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_block);
            out_bim << bims[seqs[c]] + "\n";
        }
    }
    in.close(); out.close(); out_bim.close();
    std::ifstream in_fam(fin + ".fam");
    std::ofstream out_fam(fout + ".fam");
    out_fam << in_fam.rdbuf();
}

void permute_plink(const string& fin, const string& fout,uint blocksize = 1)
{
    uint64 nsnps = count_lines(fin + ".bim");
    uint64 nsamples = count_lines(fin + ".fam");
    uint64 bed_bytes_per_snp = (nsamples+3)>>2;
    uint64 bed_bytes_per_block = bed_bytes_per_snp * blocksize;
    uint64 nblocks = (unsigned int)ceil((double)nsnps / blocksize);

    ios_base::sync_with_stdio(false);
    std::ifstream in(fin + ".bed", std::ios::binary);
    std::ofstream out(fout + ".bed", std::ios::binary);
    if (!in.is_open()) {
        throw std::invalid_argument("ERROR: Cannot open bed file.\n");
    }
    uchar header[3];
    in.read(reinterpret_cast<char *> (&header[0]), 3);
    if ( (header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01) ) {
        throw std::invalid_argument("ERROR: Incorrect magic number in bed file.\n");
    }
    out.write(reinterpret_cast<char *> (&header[0]), 3);
    vector<uchar> inbed;
    inbed.resize(bed_bytes_per_block);
    std::ifstream in_bim(fin + ".bim", std::ios::in);
    std::ofstream out_bim(fout + ".bim", std::ios::out);
    vector<string> bims;
    bims.resize(nblocks);
    vector<uint64> seqs;
    seqs.reserve(nblocks);
    string line;
    uint64 i, j;
    for(i = 0; i < nblocks; i++) {
        seqs.push_back(i);
        j = 0;
        while(getline(in_bim, line)) {
            if (j < blocksize) {
                bims[i] += line + "\n";
                j++;
                if (j == blocksize) break;
            }
        }
    }
    auto rng = std::default_random_engine {};
    std::shuffle(std::begin(seqs), std::end(seqs), rng);
    uint64 idx;
    uint64 len = 3 + bed_bytes_per_snp * nsnps;
    string sites = "";
    for(i = 0; i < nblocks; i++)
    {
        idx = 3 + (seqs[i] + 1) * bed_bytes_per_snp * blocksize;
        if (idx >= len -1) {
            bed_bytes_per_block = len - 3 - (nblocks - 1) * bed_bytes_per_snp * blocksize;
        } else {
            bed_bytes_per_block = blocksize * bed_bytes_per_snp;
        }
        idx = 3 + seqs[i] * bed_bytes_per_snp * blocksize;
        in.read(reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_block);
        out.seekp(idx, std::ios_base::beg);
        out.write(reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_block);
        sites += bims[seqs[i]];
    }
    out_bim << sites;
    in.close(); out.close(); out_bim.close();
    std::ifstream in_fam(fin + ".fam");
    std::ofstream out_fam(fout + ".fam");
    out_fam << in_fam.rdbuf();
}

int main(int argc, char *argv[])
{
    if (argc < 3) {
        cerr << "usage: permute prefix_to_bfile output_prefix [blocksize=1]" << endl;
        exit(EXIT_FAILURE);
    }
    string filein(argv[1]);
    string fileout(argv[2]);
    uint blocksize = 1;
    if (argc == 4) blocksize = atoi(argv[3]);
    permute_plink2(filein, fileout, blocksize);
}
