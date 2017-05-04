#!/usr/bin/env python

import argparse
from Bio import SeqIO
from BCBio import GFF
import sys

def tableify(gff3, fasta):
    names = {}
    for fasta_rec in SeqIO.parse(fasta, 'fasta'):
        names[fasta_rec.id] = []

    for gff_rec in GFF.parse(gff3):
        names[gff_rec.id].append(str(len(gff_rec.features)))
        starts = []
        for feat in gff_rec.features:
            starts.append(str(feat.location.start + 1))
        names[gff_rec.id].append(starts)

    for n in sorted(names):
        if len(names[n]):
            print '\t'.join([n, names[n][0], ', '.join(names[n][1])])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='outputs table of genomes and start data')
    parser.add_argument('gff3', type=argparse.FileType("r"), help='gff3 file')
    parser.add_argument('fasta', type=argparse.FileType("r"), help='fasta file')
    args = parser.parse_args()

    tableify(args.gff3, args.fasta)
