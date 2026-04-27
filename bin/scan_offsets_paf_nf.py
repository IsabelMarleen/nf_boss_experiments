#!/usr/bin/env python

import sys
sys.path.insert(0,"/")
sys.path.insert(0,"../")
from boss.sampler import PafStream

def scan_paf_offsets(paf_full, paf_trunc):
    PafStream(paf_full=paf_full, paf_trunc=paf_trunc)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python scan_offsets_paf_nf.py <paf_file_full> <paf_file_trunc>")
        sys.exit(1)

    scan_paf_offsets(paf_full=sys.argv[1], paf_trunc=sys.argv[2])