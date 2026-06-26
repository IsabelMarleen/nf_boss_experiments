#!/usr/bin/env python
import sys
import random
from contextlib import ExitStack


def combine_barcoded(output_path, num_reads, input_paths):
    with ExitStack() as stack:
        files = [
        stack.enter_context(open(filename, 'r'))
        for filename in input_paths
    ]
        with open(output_path, 'w') as out_file:
            r=0
            file_indices = list(range(0, len(input_paths)))
            while r < num_reads:
                f = random.sample(file_indices, 1)[0]
                try:
                    out_file.writelines(next(files[f]))
                    out_file.writelines(next(files[f])) 
                    out_file.writelines(next(files[f])) 
                    out_file.writelines(next(files[f]).strip()+'\n')
                except StopIteration:
                    del file_indices[f]
                    continue
                else:
                    r = r+1        

            


if __name__ == "__main__":
    random.seed(sys.argv[-1])
    combine_barcoded(output_path=sys.argv[1], num_reads=int(float(sys.argv[2])), input_paths=sys.argv[3:-1])