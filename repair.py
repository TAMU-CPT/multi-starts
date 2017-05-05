#!/usr/bin/env python

import argparse
from Bio import SeqIO
from BCBio import GFF
from Bio.Data import CodonTable
import sys

def get_CDS_and_SD(feat):
    """
        return the CDS and SD feature in gene
    """
    for sf in feat.sub_features:
        if sf.type == 'Shine_Dalgarno_sequence':
            sd = sf
        elif sf.type == 'mRNA':
            for sfmrna in sf.sub_features:
                if sfmrna.type == 'CDS':
                    cds = sfmrna

    return cds, sd

def break_start(seq):
    """
        change the start sequence while keeping the same
        amino acid translation, if possible
    """
    if seq == 'ATG':    # if methionine, no other option but ATG
        return 'ATG'
    elif seq == 'GTG':  # if valine, return GTA (still valine)
        return 'GTA'
    elif seq == 'TTG':  # if leucine, return TTA (still leucine)
        return 'TTA'

def break_sd(sd):
    # table = CodonTable.unambiguous_dna_by_id[11]
    pass

def repair(fasta, gff3):
    recs = {}
    for record in GFF.parse(gff3):
        recs[record.id] = record

    seqs = []
    for seq in SeqIO.parse(fasta, 'fasta'):
        if seq.id not in recs:
            continue

        # if seq.id == 'GFP':

            # with open('test1.fa', 'w') as seqfile1:
                # SeqIO.write(seq, seqfile1, 'fasta')

        current = recs[seq.id]
        for feat in current.features:
            cds, sd = get_CDS_and_SD(feat)
            cds_start = seq.seq[cds.location.start:cds.location.start+3]
            broken_start = break_start(cds_start)
            if cds_start != broken_start:
                # try to break start sequence while keeping amino acid the same
                seq.seq = seq.seq[0:cds.location.start] + broken_start + seq.seq[cds.location.start+3:]
            else:
                # if couldn't change start, must break SD
                break_sd(sd)
            seqs.append(seq)

    with open('out.fa', 'w') as seqfile:
        SeqIO.write(seqs, seqfile, 'fasta')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='change bases o get rid of internal starts/SDs')
    parser.add_argument('fasta', type=argparse.FileType("r"), help='fasta of gene seqs')
    parser.add_argument('gff3', type=argparse.FileType("r"), help='gff3 of starts in each seq')
    args = parser.parse_args()

    repair(args.fasta, args.gff3)
