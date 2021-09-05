#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <string>
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

void permute_plink(const string& fin, const string& fout)
{
    uint64 nsnps = count_lines(fin + ".bim");
    uint64 nsamples = count_lines(fin + ".fam");
    uint64 bed_bytes_per_snp = (nsamples+3)>>2;
    vector<uint64> seqs;
    seqs.reserve(nsnps);
    for(uint64 i = 0; i < nsnps; i++) { seqs.push_back(i); }
    // std::random_device r;
    // auto rng = std::default_random_engine {r()};
    auto rng = std::default_random_engine {};
    std::shuffle(std::begin(seqs), std::end(seqs), rng);

    std::ifstream in(fin + ".bed", std::ios::in | std::ios::binary);
    std::ofstream out(fout + ".bed", std::ios::out | std::ios::binary);
    if (!in.is_open())
    {
        cerr << "ERROR: Cannot open bed file.\n";
        exit(EXIT_FAILURE);
    }
    uchar header[3];
    in.read(reinterpret_cast<char *> (&header[0]), 3);
    if ( (header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01) ) {
        cerr << "ERROR: Incorrect magic number in bed file.\n";
        exit(EXIT_FAILURE);
    }
    out.write(reinterpret_cast<char *> (&header[0]), 3);
    vector<uchar> inbed;
    inbed.resize(bed_bytes_per_snp);
    std::ifstream in_bim(fin + ".bim", std::ios::in);
    std::ofstream out_bim(fout + ".bim", std::ios::out);
    vector<string> bims(std::istream_iterator<Line>{in_bim},
                        std::istream_iterator<Line>{});
    uint64 idx;
    for(uint64 i = 0; i < nsnps; i++)
    {
        idx = 3 + seqs[i] * bed_bytes_per_snp;
        in.seekg(idx, std::ios_base::beg);
        in.read(reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_snp);
        out.write(reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_snp);
        out_bim << bims[seqs[i]] + "\n";
    }
    in.close(); out.close(); out_bim.close();
    std::ifstream in_fam(fin + ".fam");
    std::ofstream out_fam(fout + ".fam");
    out_fam << in_fam.rdbuf();
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        cerr << "usage: permute prefix_to_bfile output_prefix" << endl;
        exit(EXIT_FAILURE);
    }
    string filein(argv[1]);
    string fileout(argv[2]);
    permute_plink(filein, fileout);
}