#!/usr/bin/env python3

import argparse
import random
from bgen import BgenReader
from bgen import BgenWriter


def main():
    parser = argparse.ArgumentParser(
        description="shuffle the order of variants in BGEN file"
    )
    parser.add_argument(
        "bgen", metavar="BGEN_IN", help="input BGEN file to be shuffled"
    )
    parser.add_argument("out", metavar="BGEN_OUT", help="output BGEN file")
    args = parser.parse_args()

    bgenr = BgenReader(args.bgen)
    print(bgenr.header.nvariants)
    print(bgenr.header.nsamples)
    bgenw = BgenWriter(args.out, n_samples=len(bgenr.samples))
    idx = list(range(0, bgenr.header.nvariants))
    random.shuffle(idx)
    for ri in idx:
        var = bgenr[ri]
        bgenw.add_variant(
            varid=var.varid,
            rsid=var.rsid,
            chrom=var.chrom,
            pos=var.pos,
            alleles=var.alleles,
            genotypes=var.probabilities.astype("float64"),
        )
    bgenw.close()


if __name__ == "__main__":
    main()
