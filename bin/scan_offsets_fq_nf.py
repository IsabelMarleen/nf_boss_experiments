#!/usr/bin/env python

import sys
sys.path.insert(0,"/")
sys.path.insert(0,"../")
sys.path.insert(0,".")
from boss.sampler import FastqStream_mmap


def scan_fq_offsets(fq_file):
    FastqStream_mmap(source=str(fq_file), shuffle=True)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python scan_offsets_fq_nf.py <fq_file>")
        sys.exit(1)

    scan_fq_offsets(fq_file=sys.argv[1])