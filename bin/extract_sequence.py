#!/usr/bin/env python

import sys
from os import path
PARENT_DIR = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(PARENT_DIR)
sys.path.append('.')
from otu_info import otus_no_plasmids


def extract_sequence(fasta_file: str, header: str) -> str:
    '''
    Extracts the sequence of a specific header in a FASTA file.
    :param fasta_file: path to fasta file
    :param header: header of the sequence to extract
    :return: the requested sequence
    '''
    with open(fasta_file, 'r') as file:
        sequence = ""
        record = False
        for line in file:
            if line.startswith(">"):
                if record:
                    break
                # use only the first bit of the header, i.e. whatever is before the first space
                if line[1:].split(' ')[0].strip() == header:
                    record = True
            elif record:
                sequence += line.strip()
    return sequence



if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_sequence.py <fasta_file> <header>")
        sys.exit(1)
    
    otu = sys.argv[2]
    otus_no_plasmids_rev = {v: k for k, v in otus_no_plasmids.items()}
    otu_raw=(lambda w: otus_no_plasmids_rev[w])(otu)

    header = otu_raw
    sequence = extract_sequence(fasta_file=sys.argv[1], header=header)
    if sequence:
        print(f'>{header}\n{sequence}')
    else:
        print('no sequence found')
        sys.exit(1)



