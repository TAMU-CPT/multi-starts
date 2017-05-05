#!/usr/bin/env python

import sys
import argparse
import svgwrite
from BCBio import GFF


def draw(recs):
    for i, rec in enumerate(list(recs)):
        height = 30 + 15*(len(list(rec.features)))
        dwg = svgwrite.Drawing(filename='svgs/'+ rec.id + '.svg', size=("1500px", "%spx" % height), debug=True)
        genes = dwg.add(dwg.g(id='genes', fill='white'))

        y = 10

        # draw full genome first in black
        full = 0
        for r in rec.annotations['sequence-region']:
            if rec.id in r:
                full = r[2]

        genes.add(dwg.rect(insert=(10, y), size=(full, 10), fill='black'))
        y += 15

        for j, gene in enumerate(rec.features):
            # draw each gene in blue
            x = 10 + gene.location.start
            length = gene.location.end - gene.location.start
            genes.add(dwg.rect(insert=(x, y), size=(length, 10), fill='blue'))

            for sf in gene.sub_features:
                if sf.type == 'Shine_Dalgarno_sequence':
                    x = 10 + sf.location.start
                    length = sf.location.end - sf.location.start
                    genes.add(dwg.rect(insert=(x, y), size=(length, 10), fill='blue'))

            y += 15

            dwg.save()
        y += 20


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Intron detection')
    parser.add_argument('gff3', type=argparse.FileType("r"), help='gff3 file')
    args = parser.parse_args()

    records = GFF.parse(args.gff3)
    draw(records)
