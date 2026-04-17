#!/usr/bin/env python

import numpy as np
import scipy.integrate as integrate
import subprocess
import sys

from os import path
PARENT_DIR = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(PARENT_DIR)
sys.path.append('.')

from otu_info import genome_sizes   # noqa



def parse_pup(pup_path: str, pup_full:str, count: int) -> np.array:
    '''
    Parse the pileup file to extract the coverage values as an array
    :param pup_path: path to the pileup file
    :return: array of coverage
    '''
    # determine which columns need to be extracted
    # 2 is the position, then the 4th column is the first coverage and then there are two columns in-between all the remaining coverages
    col_string = "2,"+','.join(map(str, range(4, 4+(3*count), 3)))

    # extract all relevant columns and output with comma separation
    running2 = subprocess.Popen(
        args=f"cut -f{col_string} --output-delimiter=',' {pup_path} > tmp.csv",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='utf-8',
        shell=True
    )
    stdout2, stderr2 = running2.communicate()

    cov = np.loadtxt("tmp.csv", delimiter=",", dtype="int")
    # stdout4_interim = [subprocess.Popen(
    #     args=f"cut -f{i} {pup_path}",
    #     stdout=subprocess.PIPE,
    #     stderr=subprocess.PIPE,
    #     encoding='utf-8',
    #     shell=True
    # ) for i in range(4, 4+(3*count), 3)]
    # stdout4 = [[j for j in j.communicate()[0].split("\n") if j] for j in stdout4_interim ]

    # [[subprocess.Popen(
    #     args=f"cut -f{i} {pup_path}",
    #     stdout=subprocess.PIPE,
    #     stderr=subprocess.PIPE,
    #     encoding='utf-8',
    #     shell=True
    # ) for i in range(4, 4+(3*count), 3)] for j in j.communicate()[0].split("\n") if j]
    

    # running4 = subprocess.Popen(
    #     args=f"cut -f4 {pup_path}",
    #     stdout=subprocess.PIPE,
    #     stderr=subprocess.PIPE,
    #     encoding='utf-8',
    #     shell=True
    # )

    # stdout4, stderr4 = running4.communicate()
    # parse to array of coverage
    # pos = np.array([i for i in stdout2.split("\n") if i], dtype="int")
    # coverage = np.array(stdout4, dtype="int")

    # cov = np.array([pos, coverage], dtype="int")
    
    return cov #.transpose()


def calc_lowcov(cov: np.array, otu_size: int, lowcov_threshold: int = 5) -> float:
    '''
    Calculate the proportion of low coverage sites in the OTU.
    Attention the coverage array does not contain sites without any coverage

    :param cov: array of coverage counts
    :param otu_size: genome size
    :return: proportion of sites with low coverage
    '''
    # lowcov_n = np.where(cov < lowcov_threshold)[0].shape[0]
    lowcov_n = np.sum(cov < lowcov_threshold, axis=0)
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
        avg_cov = np.average(cov, axis=0)
        cov_norm = cov/avg_cov
        p_i = np.count_nonzero(cov_norm >= (i/avg_cov))

        evn = p_i/otu_size
        return evn

    evenness = [integrate.quad(_calc_evn, 0, 1, args=(cov[:,i], otu_size))[0] for i in range(0, cov.shape[1])]
    
    return evenness



def process_pup(pup: str, otu_size: int, pup_full: str, count: int):
    '''
    Process the pileup file to extract the mean coverage and proportion of low coverage sites
    :param pup: path of pileup file
    :param otu_size: genome size
    :return: csv data
    '''
    # get metadata
    meta = pup.split("/")[-1].split(".")[0].split('_')
    if len(meta) > 2:
        cond, bc, otu = meta[0], meta[1], meta[2]
    else:
        cond, otu = meta[0], meta[1]
    # get an array of coverage from the pileup file
    cov = parse_pup(pup_path=pup, pup_full=pup_full, count=count)
    cov_cumsum = np.empty_like(cov[:,1:], dtype="int")
    cov[:, 1:].cumsum(axis=1, out=cov_cumsum)
    print("cumsum complete")
    # np.savetxt(pup_full, cov_cumsum)
    # pup_full_2 = "/hps/nobackup/goldman/ipoetzsch/boss_experiments/work/de/87802f7c6910861995d08dc607c94f/hschr21"
    # cov = np.loadtxt(pup_full)
    # calculate the mean coverage of the OTU
    # mean_cov = np.sum(cov[:, 1:], axis=0) / otu_size
    mean_cov = np.sum(cov_cumsum, axis=0) / otu_size
    # calculate the proportion of low coverage sites
    # lowcov_p = calc_lowcov(cov=cov[1], otu_size=otu_size)
    lowcov_p = calc_lowcov(cov=cov_cumsum, otu_size=otu_size)

    # calculate the evenness of coverage of the OTU
    evenness = calc_evenness(cov=cov_cumsum, otu_size=otu_size)
    with open(pup_full, "w") as output:
        if len(meta) > 2:
            output.write(f'{cond},0,{otu},{bc},0,1,0\n')
            for i in range(1, count):
                output.write(f'{cond},{i},{otu},{bc},{mean_cov[i]},{lowcov_p[i]},{evenness[i]}\n')
        else:
            output.write(f'{cond},0,{otu},0,1,0\n')
            for i in range(1, count):
                output.write(f'{cond},{i},{otu},{mean_cov[i]},{lowcov_p[i]},{evenness[i]}\n')
    return




if __name__ == "__main__":
    # input to this script is the output of samtools pileup
    # also pass in the genome size of the OTU to calculate the mean coverage
    otu = sys.argv[2]
    otu_size = (lambda w: genome_sizes[w])(otu)

    process_pup(pup=sys.argv[1], otu_size=otu_size, pup_full=sys.argv[3], count=int(sys.argv[4]))







