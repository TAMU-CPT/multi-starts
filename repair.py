#!/usr/bin/env python

import argparse
from Bio import SeqIO
from BCBio import GFF
import sys

def get_CDS_and_SD(feat):
    for sf in feat.sub_features:
        if sf.type == 'Shine_Dalgarno_sequence':
            sd = sf
        elif sf.type == 'mRNA':
            for sfmrna in sf.sub_features:
                if sfmrna.type == 'CDS':
                    cds = sfmrna

    return cds, sd

def repair(fasta, gff3):
    recs = {}
    for record in GFF.parse(gff3):
        recs[record.id] = record
    for seq in SeqIO.parse(fasta, 'fasta'):
        if seq.id not in recs:
            continue

        if seq.id == 'GFP':
            current = recs[seq.id]
            for feat in current.features:
                cds, sd = get_CDS_and_SD(feat)
                print cds.location.start, sd.location.start


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='change bases o get rid of internal starts/SDs')
    parser.add_argument('fasta', type=argparse.FileType("r"), help='fasta of gene seqs')
    parser.add_argument('gff3', type=argparse.FileType("r"), help='gff3 of starts in each seq')
    args = parser.parse_args()

    repair(args.fasta, args.gff3)
