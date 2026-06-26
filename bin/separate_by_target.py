#!/usr/bin/env python

import sys
from os import path
import re
PARENT_DIR = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(PARENT_DIR)
sys.path.append('.')

import numpy as np

from otu_info import otus_clean_names   # noqa
from readfq import readfq


"""
this script is used to separate a fq file into multiple fqs given the mappings in a paf or sam file 
- input: fq & paf/sam
- output: multiple fqs
"""


def parse_aln(aln_path: str, target_dict: dict) -> dict:
    '''
    parse the paf or sam file and return a dictionary with read ids and their target organism
    :param aln_path: path to alignment file
    :param target_dict: dictionary with target reference sequences
    :return: dict of read ids and their target sequences
    '''
    # open the paf and record the source for each mapped read
    read_targets = {}
    target_barcodes = {}
    read_barcodes = {}
    rid_aln_scores = {}

    ext = aln_path.split('.')[-1]
    if ext == "paf":
        tpos = 5
    elif ext == "sam":
        tpos = 2
    else:
        print("unknown extension, exit")
        sys.exit()

    with open(aln_path, 'r') as aln:
        for line in aln:
            if line.startswith('@'):
                continue
            ll = line.split('\t')

            target = ll[tpos]
            barcode_search = re.search(r"(unclassified|barcode([0-9]+))", ll[0])
            if barcode_search is not None:
                rid = re.sub(r".(unclassified|barcode([0-9]+))", "", ll[0])
                read_barcodes[rid] = int(barcode_search[0].split('barcode')[1])
                try:
                    target_dict[target] # Check whether the target is one we are interested in
                except KeyError:  
                    continue
                try:
                    target_barcodes[target_dict[target]].append(int(barcode_search[0].split('barcode')[1]))
                except KeyError:
                    target_barcodes[target_dict[target]] = []
                    target_barcodes[target_dict[target]].append(int(barcode_search[0].split('barcode')[1]))
                # if the pattern is not in the header, skip the read
            else:
                rid = ll[0]
            
            aln_score = int(ll[14].split(':')[-1])
            # if the rid is not recorded yet, just add it
            if rid not in read_targets.keys():
                try:
                    read_targets[rid] = target_dict[target]
                except KeyError:
                    continue
                rid_aln_scores[rid] = aln_score
            else:
                # if the rid is already recorded, check whether this one has higher aln score
                # and only overwrite if its higher
                if aln_score > rid_aln_scores[rid]:
                    try:
                        read_targets[rid] = target_dict[target]
                    except KeyError:
                        continue
    return read_targets, read_barcodes, target_barcodes


def write_target_reads(read_targets: dict, fq: str, target_dict: dict, read_barcodes: dict, target_barcodes:dict) -> None:
    '''
    open a fq file and write reads to individual files based on the target
    :param read_targets: dict of targets for each read
    :param fq: file with the read data
    :param target_dict: dict of names for output files
    :return: None
    '''
    # open a fq file for each target & barcode
    targets = np.unique(list(target_dict.values()))

    fq_base = fq.split('/')[-1].split('.')[0]

    # open a file for each target
    target_files = dict()
    for t in targets:
        if target_barcodes == {}:
            target_files[t] = open(f'{fq_base}_{t}.fq', 'w')
        else:
            for b in np.unique(target_barcodes[t]):
                target_files[t+str(b)] = open(f'{fq_base}_bc{b}_{t}.fq', 'w')

    # iterate over file and write the fq into the correct file
    with open(fq, 'r') as fastq:
        for desc, name, seq, qual in readfq(fastq):
            # check which target the read comes from
            if target_barcodes == {}:
                try:
                    t = read_targets[name]
                except KeyError:
                    continue
            else:
                try:
                    rid = re.sub(r".(unclassified|barcode([0-9]+))", "", name)
                    t = read_targets[rid]
                    b = read_barcodes[rid]
                except KeyError:
                    continue
            # write to file
            if not qual:
                qual = 'z' * len(seq)
            
            if target_barcodes == {}:
                target_files[t].write(f'@{desc}\n{seq}\n+\n{qual}\n')
            else:
                target_files[t+str(b)].write(f'@{desc}\n{seq}\n+\n{qual}\n')

    # close all files
    for tf in target_files.values():
        tf.close()



if __name__ == "__main__":
    # load a dictionary with targets and read ids
    read_targets, read_barcodes, target_barcodes = parse_aln(aln_path=sys.argv[1], target_dict=otus_clean_names)
    # open the fastq, iterate through it and write reads to individual files
    write_target_reads(read_targets=read_targets, fq=sys.argv[2], target_dict=otus_clean_names, read_barcodes=read_barcodes, target_barcodes=target_barcodes)


