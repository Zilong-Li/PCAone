#include "FileBinary.hpp"

using namespace std;

bool isZstdCompressed(const char * filename)
{
    FILE * file = fopen(filename, "rb");
    if(!file) return false;

    char magicNumber[4];
    if(fread(magicNumber, 1, 4, file) != 4)
    {
        fclose(file);
        return false;
    }

    bool isCompressed = (ZSTD_isFrame(magicNumber, 4) != 0);

    fclose(file);
    return isCompressed;
}

void FileBin::check_file_offset_first_var()
{
    setlocale(LC_ALL, "C");
    ios_base::sync_with_stdio(false);
    // magic += missing_points.size() * sizeof(uint64);
    long long offset = magic + nsnps * bytes_per_snp;
    if(ifs_bin.tellg() == offset)
    {
        // reach the end of bed, reset the position to the first variant;
        ifs_bin.seekg(magic, std::ios_base::beg);
    }
    else if(ifs_bin.tellg() == magic)
    {
        ;
    }
    else
    {
        ifs_bin.seekg(magic, std::ios_base::beg);
        if(params.verbose) cao.warning("make sure you are runing PCAone algorithm2");
    }
}

void FileBin::read_all()
{
    check_file_offset_first_var();
    G = MyMatrix::Zero(nsamples, nsnps);
    Eigen::VectorXf fg(nsamples);
    for(size_t i = 0; i < G.cols(); i++)
    {
        ifs_bin.read((char *)fg.data(), bytes_per_snp);
        G.col(i) = fg.cast<double>();
        G.col(i).array() -= G.col(i).mean();
    }
}

void FileBin::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize)
{
    // magic += missing_points.size() * sizeof(uint64);
    // check where we are
    long long offset = magic + start_idx * bytes_per_snp;
    if(ifs_bin.tellg() != offset) throw std::runtime_error("Error: something wrong with read_snp_block!\n");
    uint actual_block_size = stop_idx - start_idx + 1;
    G = MyMatrix(nsamples, actual_block_size);
    Eigen::VectorXf fg(nsamples);
    for(size_t i = 0; i < G.cols(); i++)
    {
        ifs_bin.read((char *)fg.data(), bytes_per_snp);
        G.col(i) = fg.cast<double>();
        G.col(i).array() -= G.col(i).mean();
    }
}
