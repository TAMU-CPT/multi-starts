#!/usr/bin/env python

# import os
import argparse
from Bio import SeqIO
import glob

def sort_and_join(path):
    good_starts = []
    bad_starts = []
    # for fn in os.listdir(path):
    for fn in glob.glob(path+'*'):
        if not fn.endswith('fasta'):
            continue

        for record in SeqIO.parse(fn, 'fasta'):
            if record.seq.startswith('ATG') or record.seq.startswith('TTG') or record.seq.startswith('GTG'):
                good_starts.append(record)
            else:
                bad_starts.append(record)

        with open('good_starts.fa', 'w') as good_starts_file:
            SeqIO.write(good_starts, good_starts_file, "fasta")
        with open('bad_starts.fa', 'w') as bad_starts_file:
            SeqIO.write(bad_starts, bad_starts_file, "fasta")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='sort fasta files based on good start codon. combine them.')
    parser.add_argument('path', type=str, help='dir with fastas')
    args = parser.parse_args()

    sort_and_join(args.path)
