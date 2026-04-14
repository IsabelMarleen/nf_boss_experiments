#!/usr/bin/env python

import numpy as np
import scipy.integrate as integrate
import subprocess
import sys

from os import path
PARENT_DIR = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(PARENT_DIR)
sys.path.append('.')
sys.path.append('/hps/nobackup/goldman/ipoetzsch/boss_experiments/work/67/4ce889ccbb44dc05e8e07ef9630530')
sys.path.append('/hps/nobackup/goldman/ipoetzsch/boss_experiments/work/b4/1901098fd43d5a735be8d048edd350')

from otu_info import genome_sizes   # noqa


def calc_lowcov(cov: np.array, otu_size: int, lowcov_threshold: int = 5) -> float:
    '''
    Calculate the proportion of low coverage sites in the OTU.
    Attention the coverage array does not contain sites without any coverage

    :param cov: array of coverage counts
    :param otu_size: genome size
    :return: proportion of sites with low coverage
    '''
    lowcov_n = np.where(cov < lowcov_threshold)[0].shape[0]
    # to account for sites without coverage
    # calc difference in genome length and cov array
    nocov_n = otu_size - cov.shape[0]
    if nocov_n > 0:
        lowcov_n += nocov_n
    lowcov_p = lowcov_n / otu_size
    return lowcov_p

def calc_evenness(cov: np.array, otu_size: int) -> float:
    '''
    Calculate the evenness of coverage sites in the OTU
    as described in Mokry et al, 2010, Nucleic Acids Research.

    :param cov: array of coverage counts
    :param otu_size: genome size
    :return: evenness of coverage
    '''
    def _calc_evn(i, cov, otu_size) -> float:
        avg_cov = np.average(cov)
        cov_norm = cov/avg_cov
        p_i = np.count_nonzero(cov_norm >= (i/avg_cov))

        evn = p_i/otu_size
        return evn

    evenness = integrate.quad(_calc_evn, 0, 1, args=(cov, otu_size))
    
    return evenness[0]


        # mean_cov = np.sum(cov[1]) / otu_size
    # # calculate the proportion of low coverage sites
    # lowcov_p = calc_lowcov(cov=cov[1], otu_size=otu_size)
    # # calculate the evenness of coverage of the OTU
    # evenness = calc_evenness(cov=cov[1], otu_size=otu_size)
    # if len(meta) > 3:
    #     output = f'{cond},{time},{otu},{bc},{mean_cov},{lowcov_p},{evenness}'
    # else:
    #     output = f'{cond},{time},{otu},{mean_cov},{lowcov_p},{evenness}'

def complete_pup(pup, otu_size):
    p = np.stack([np.loadtxt(p) for p in pup])

    p = np.cumsum(p)

    mean_cov = np.sum(p, axis=(1,2))/otu_size
    print(mean_cov)
        # mean_cov = np.sum(cov[1]) / otu_size
    # # calculate the proportion of low coverage sites
    # lowcov_p = calc_lowcov(cov=cov[1], otu_size=otu_size)
    # # calculate the evenness of coverage of the OTU
    # evenness = calc_evenness(cov=cov[1], otu_size=otu_size)
    # if len(meta) > 3:
    #     output = f'{cond},{time},{otu},{bc},{mean_cov},{lowcov_p},{evenness}'
    # else:
    #     output = f'{cond},{time},{otu},{mean_cov},{lowcov_p},{evenness}'


if __name__ == "__main__":
    # input to this script is the output of samtools pileup
    # also pass in the genome size of the OTU to calculate the mean coverage
    otu = sys.argv[1]
    otu_size = (lambda w: genome_sizes[w])(otu)

    # data = 
    complete_pup(pup = sys.argv[2:], otu_size=int(otu_size))
    # print(data)