#include "FileBinary.hpp"

using namespace std;

void FileBin::check_file_offset_first_var()
{
    setlocale(LC_ALL, "C");
    ios_base::sync_with_stdio(false);
    // magic += missing_points.size() * sizeof(uint64);
    long long offset = magic + nsnps * bytes_per_snp;
    if (ifs_bin.tellg() == offset)
    {
        // reach the end of bed, reset the position to the first variant;
        ifs_bin.seekg(magic, std::ios_base::beg);
    }
    else if (ifs_bin.tellg() == magic)
    {
        ;
    }
    else
    {
        ifs_bin.seekg(magic, std::ios_base::beg);
        if (params.verbose)
            std::cout << colwarn + "make sure you are runing PCAone algorithm2" + colend << std::endl;
    }
}

void FileBin::read_all()
{
    uint64 i, j;
    check_file_offset_first_var();
    G = MyMatrix::Zero(nsamples, nsnps);
    for (i = 0; i < nsnps; ++i)
    {
        for (j = 0; j < nsamples; j++)
        {
            ifs_bin.read(reinterpret_cast<char*>(&G(j, i)), sizeof(double));
        }
    }
}

void FileBin::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize)
{
    // check where we are
    if (params.verbose)
    {
        // magic += missing_points.size() * sizeof(uint64);
        long long offset = magic + start_idx * bytes_per_snp;
        if (ifs_bin.tellg() != offset)
            throw std::runtime_error("Error: something wrong with read_snp_block!\n");
    }
    uint actual_block_size = stop_idx - start_idx + 1;
    G = MyMatrix(actual_block_size, nsamples);
    ifs_bin.read(reinterpret_cast<char*>(&G[0]), nsamples * actual_block_size);
    G.transposeInPlace(); // we want column-major
}

void FileBin::read_block_update(uint64 start_idx, uint64 stop_idx, const MyMatrix& U, const MyVector& svals, const MyMatrix& VT, bool standardize)
{
}
