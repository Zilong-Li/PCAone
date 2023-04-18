#!/usr/bin/env python3

import argparse
from bgen import BgenReader
from bgen import BgenWriter


def main():
    parser = argparse.ArgumentParser(description="concat bgen files into a big one.")
    parser.add_argument("--bgen", nargs="+", help="List of BGEN files")
    parser.add_argument("--out", help="output BGEN file")
    args = parser.parse_args()

    bgenr_fp = [BgenReader(fn) for fn in args.bgen]
    bgenw = BgenWriter(args.out, n_samples=len(bgenr_fp[0].samples))

    for bgenr in bgenr_fp:
        for var in bgenr:
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
