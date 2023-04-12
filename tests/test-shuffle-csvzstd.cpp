#include "../src/FileCsv.hpp"
#include "catch.hh"

using namespace std;

TEST_CASE("convert csv zstd to binary", "[test-shuffle-csvzstd]")
{
    int ec = shuffle_csvzstd_to_bin("../examples/test.csv.zst", "test.bin", 2, true);
    REQUIRE(ec == 1);
    ec = shuffle_csvzstd_to_bin("../examples/test.csv.zst", "test.bin", 2, false);
    REQUIRE(ec == 1);
}
