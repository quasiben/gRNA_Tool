import os
import re
import sys
import csv
import argparse
from glob import glob
from collections import defaultdict
from itertools import izip

import fr_lib
import wtlib_lib

descr = "CLI tool for gRNA analysis"


def parse_args():
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("wtlib", type=str, help="wtlib file")
    parser.add_argument("--fwd", type=str, help="forward read file")
    parser.add_argument("--rev", type=str, help="reverse read file")
    parser.add_argument("-d", "--dir",
                        type=str,
                        help="read directory of files.  "
                             "Files must be of match: *R1*.gz, *R2*.gz",
                        )
    parser.add_argument("--parallel", action="store_true", 
		        help="multicore support with Dask")
    # parser.add_argument("-o", "--output",
    #                     default="output.csv",
    #                     type=str,
    #                     help="output csv file name")

    return parser.parse_args()


def process_gRNA(forward_path, reverse_path):
    """process gRNA

    Args:
      forward_path (str): forward file path
      reverse_path (str): reverse file path

    Output:
      None: write to disk with fwd_name.csv

    """
    print "Processing Files: ", forward_path, reverse_path

    wtlib_dict = wtlib_lib.get_wtlib_index()
    counter_dict = defaultdict(int)

    MIs = ["ATCACG", "CGATGT", "TTAGGC", "TGACCA"]

    for f_seq, r_seq in fr_lib.streaming_open(forward=forward_path,
                                              reverse=reverse_path):

        spacer = fr_lib.get_spacer(f_seq)
        barcode, mi = fr_lib.get_barcode_mi(r_seq)

        if wtlib_dict.get(spacer):

            # skip if mi is invalid
            if mi not in MIs:
                continue

            hash_name = '$'.join([spacer, barcode, mi])
            counter_dict[hash_name] += 1

    # output file clean up :)
    output_file = os.path.basename(forward_path)  # get file name only
    output_file = output_file.split('.')[0]  # remove extensions
    #output_file = re.sub("_R1.*", "", output_file)  # section out _R1*
    output_file = output_file.replace("_R1", "")  # section out _R1*
    output_file = output_file+".csv"

    with open(output_file, 'w') as csvfile:
        print "\tWriting Output: ", output_file

        header = ["wtLibname", "spacer", "MI", "barcode", "count"]
        frwriter = csv.writer(csvfile, delimiter=',')
        frwriter.writerow(header)

        for k, count in counter_dict.iteritems():
            spacer, barcode, mi = k.split("$")

            wtLibname, full_seq = wtlib_dict[spacer]
            # print([wtLibname, spacer, mi, barcode, count])
            frwriter.writerow([wtLibname, spacer, mi, barcode, count])


def main():
    args = parse_args()

    if args.dir and (args.fwd or args.rev):
        sys.exit("Can define both dir and file options.  Please choose one")

    elif args.dir:
        files = glob(args.dir + "/*R1*.gz")
        files.sort()

    else:
        fwd_file = args.fwd
        rev_file = args.rev
        files = [fwd_file, rev_file]
    
    files = [(f, f.replace("_R1", "_R2")) for f in files] 

    if args.parallel:
        from dask.bag import Bag
        b = Bag.from_sequence(files)
        list(b.map(lambda f, r: process_gRNA(f,r)))
    else:
        for f,r  in files:
            process_gRNA(f, r)

if __name__ == '__main__':
    main()
