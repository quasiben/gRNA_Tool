"""
small library for handling forward.fastq and reverse.fastq files
and directories of files
"""

import os
import gzip


def streaming_open(forward=None, reverse=None):
    """stream forward file in chunks of 4 lines
    stream reverse file in chunks of 3 lines

    Args:
      forward (str): foward file gzipped
      reverse (str): reverse file gzipped
    """

    with gzip.open(forward, 'rb') as f:
        with gzip.open(reverse, 'rb') as r:
            while True:
                name, f_seq, _, comment = [f.readline().strip() for i in range(4)]
                name, r_seq, _, something = [r.readline().strip() for i in range(4)]
                if name == "":
                    break
                yield (f_seq, r_seq)


def reverse_complement(seq):
    """reverse complement
    Args:
      seq (str): sequence of nucleotides

    Returns:
      str: reversed sequence of nucleotides
    """

    seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return "".join([seq_dict[base] for base in reversed(seq)])


def get_spacer(seq):
    """Function which returns 20 nt, know

    Args:
      seq (str): sequence of nucleotides

    Returns:
      str: 20 nt -- spacer

    Examples:
      >>> get_spacer("NGGGGCCCTCAGATCTCTCGTGTTTAAGAGCTATGCTGGAAACAGCATAGC")
      'NGGGGCCCTCAGATCTCTCG'

    """

    return seq[:20]


def get_barcode_mi(seq):
    """Function which returns barcode 20 nt, unknow and
    MI (know)

    Args:
      seq (str): sequence of nucleotides

    Returns:
      tuple: (20 nt, 6n) -- (barcode, mi)

    Examples:
      >>> get_barcode_mi("GCTATCCAAGTGCCTACCAATTAGGCTGCATGGCGGTAATACGGTTATCCA")
      'GGGGCCCTCAGATCTCTCGT'

    """

    barcode = seq[:20]
    mi = seq[20:26]
    return (barcode, mi)

if __name__ == '__main__':
    ROOT = "./"
    forward_path = os.path.join(ROOT, "forward.fastq.gz")
    reverse_path = os.path.join(ROOT, "reverse.fastq.gz")

    print("One Chunk")
    f_seq, r_seq = streaming_open(forward=forward_path,
                                  reverse=reverse_path).next()

    print(f_seq, r_seq)
    assert get_barcode_mi(r_seq) == ('GCTATCCAAGTGCCTACCAA', 'TTAGGC')
    assert get_spacer(f_seq) == 'NGGGGCCCTCAGATCTCTCG'
