import os
import csv
import argparse
from collections import defaultdict

import fr_lib
import wtlib_lib

descr = "CLI tool for gRNA analysis"


def parse_args():
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("wtlib", type=str, help="wtlib file")
    parser.add_argument("--fwd", type=str, help="forward read file")
    parser.add_argument("--rev", type=str, help="reverse read file")
    parser.add_argument("-o", "--output",
                        default="output.csv",
                        type=str,
                        help="output csv file name")

    return parser.parse_args()


def main():
    args = parse_args()

    forward_path = args.fwd
    reverse_path = args.rev

    wtlib_dict = wtlib_lib.get_wtlib_index()
    counter_dict = defaultdict(int)

    MIs = ["ATCACG", "CGATGT", "TTAGGC", "TGACCA"]

    for f_seq, r_seq in fr_lib.streaming_open(forward=forward_path, reverse=reverse_path):

        spacer = fr_lib.get_spacer(f_seq)
        barcode, mi = fr_lib.get_barcode_mi(r_seq)
        spacer_lookup = spacer[1:] # offset because N precedes many spacers

        if wtlib_dict.get(spacer_lookup):

            # skip if mi is invalid
            if mi not in MIs:
                continue

            hash_name = '$'.join([spacer, barcode, mi])
            counter_dict[hash_name] += 1

    output_file = args.output
    with open(output_file, 'w') as csvfile:

        header = ["wtLibname", "spacer", "MI", "barcode", "count"]
        frwriter = csv.writer(csvfile, delimiter=',')
        frwriter.writerow(header)

        for k, count in counter_dict.iteritems():
            spacer, barcode, mi = k.split("$")
            spacer_lookup = spacer[1:] # offset because N precedes many spacers

            wtLibname, full_seq = wtlib_dict[spacer_lookup]
            # print([wtLibname, spacer, mi, barcode, count])
            frwriter.writerow([wtLibname, spacer, mi, barcode, count])


if __name__ == '__main__':
    main()
