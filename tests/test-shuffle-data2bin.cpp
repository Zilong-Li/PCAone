#include "../src/FileBgen.hpp"
#include "../src/FileCsv.hpp"
#include "catch.hh"

using namespace std;

TEST_CASE("convert CSV ZSTD to binary", "[test-shuffle-csvzstd]")
{
    int ec = shuffle_csvzstd_to_bin("../examples/test.csv.zst", "test.bin", 2, true, true);
    REQUIRE(ec == 1);
    ec = shuffle_csvzstd_to_bin("../examples/test.csv.zst", "test.bin", 2, false, true);
    REQUIRE(ec == 1);
}

TEST_CASE("convert BGEN to binary", "[test-shuffle-bgen]")
{
    int ec = shuffle_bgen_to_bin("../examples/test.bgen", "test.bin", 2, true);
    REQUIRE(ec == 1);
}

TEST_CASE("permute BGEN to BGEN", "[test-permute-bgen]")
{
    string fin = "../examples/test.bgen";
    string fout = "test.bgen";
    permute_bgen(fin, fout);
}
