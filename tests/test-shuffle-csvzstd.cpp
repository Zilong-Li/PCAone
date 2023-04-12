#include "../src/FileCsv.hpp"
#include "catch.hh"


TEST_CASE("convert csv zstd to binary", "[test-shuffle-csvzstd]")
{
    shuffle_csvzstd_to_bin("../examples/test.csv.zst", "test.bin", 100, true, true);
}
