#!/usr/bin/env python

import argparse
from Bio import SeqIO


def repair(fasta):
    for record in SeqIO.parse(fasta, 'fasta'):
        print record.id


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='change bases o get rid of internal starts/SDs')
    parser.add_argument('fasta', type=argparse.FileType("r"), help='fasta of gene seqs')
    args = parser.parse_args()

    repair(args.fasta)
