import os
from tables import *


def create_wtLibHDF5():
    class wtLib(IsDescription):
        full_seq  = StringCol(35)      # full sequence 35 chars
        spacer    = StringCol(20)      # space 20 chars
        name      = StringCol(40)      # file name 33 chars

    h5file = open_file("wtLib.h5", mode="w", title="Indexed wtLib File")
    group = h5file.create_group("/", 'data', 'all data')
    table = h5file.create_table(group, 'seq', wtLib, "wtLib Look Up")
    wtlib = table.row

    with open("wtLib.fa") as f:
        print("Storing wtLib in HDF5")
        while True:
            name, full_seq = [f.readline().strip() for i in range(2)]
            if name == "":
                break
            spacer = full_seq[1:21]

            wtlib['name'] = name
            wtlib['full_seq'] = full_seq
            wtlib['spacer'] = spacer
            wtlib.append()

    table.flush()
    table.cols.spacer.createCSIndex()
    condition = '(spacer == "GAGGTGTATTGAACACGTGG")'
    row = [(x['name'], x['full_seq']) for x in table.where(condition)]
    print(row[0])
    h5file.close()


def get_wtlib_index():
    """Convenience function to return in memory
    dictionary of wtlib data. index is the spacer

    Args:
      None

    Returns:
      dict:
    """

    d = {}
    with open("wtLib.fa") as f:
        print("Storing wtLib in dict")
        while True:
            name, full_seq = [f.readline().strip() for i in range(2)]
            if name == "":
                break

            # offset because space in forward file
            # is often preceded with N
            spacer = full_seq[1:20]

            d[spacer] = (name, full_seq)

    return d


def spacer_lookup(spacer, table=None):
    """Lookup spacer in HDF5 file and return row

    Args:
      spacer (str): 20 nt string
      table (PyTable Object): pytables table of indexed wtLib.fa file

    Returns:
      tuple: name, full_seq
    """

    condition = '(spacer == "{}")'.format(spacer)

    # row should only have one element
    row = [(x['name'], x['full_seq']) for x in table.where(condition)]

    return row[0]


if __name__ == '__main__':

    if not os.path.exists("wtLib.h5"):
        create_wtLibHDF5()

    h5file = open_file("wtLib.h5", mode="r")
    table = h5file.root.data.seq
    data = spacer_lookup("CCGCCCCCAAACCTCGAGCG", table)
    print(data)

    test_name = ">hsa-mir-219a-1_hsa_EDIT_wgcod_1__015026"
    name, full_seq = data
    assert name == test_name

    name, full_seq = get_wtlib_index()['CCGCCCCCAAACCTCGAGCG']
    print(name)
    assert name == test_name
