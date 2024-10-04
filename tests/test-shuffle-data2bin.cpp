#define _DECLARE_TOOLBOX_HERE
#include "../src/Utils.hpp"
#include "../src/FileBgen.hpp"
#include "../src/FileCsv.hpp"
#include "zstd.h"
#include "catch.hh"

using namespace std;

TEST_CASE("convert CSV ZSTD to binary", "[test-shuffle-csvzstd]")
{
    string fin = "../example/test.csv.zst";
    string fout = "test.bin";
    auto perm = shuffle_csvzstd_to_bin(fin, fout, 2, true);
    std::cout << perm.indices().transpose();
}

TEST_CASE("convert BGEN to binary", "[test-shuffle-bgen]")
{
    string fin = "../example/test.bgen";
    string fout = "test.bin";
    int ec = shuffle_bgen_to_bin(fin, fout, 2, true);
    REQUIRE(ec == 1);
}

TEST_CASE("permute BGEN to BGEN", "[test-permute-bgen]")
{
    string fin = "../example/test.bgen";
    string fout = "test.bgen";
    auto perm = permute_bgen(fin, fout, 2);
    std::cout << perm.indices().transpose();
}
