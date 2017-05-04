#!/usr/bin/env python
# Copyright (c) 2015-6, Eric Rasche <esr@tamu.edu>
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   * Neither the name of the <organization> nor the
#     names of its contributors may be used to endorse or promote products
#     derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


import StringIO
import binascii
import argparse
import sys
import os

format_dictionary = {
    0: 'DNA',
    8: 'Additional Sequence Properties',  # XML
    10: 'Features',
    5: 'Primers',
    6: 'Notes',  # XML
    17: 'Alignable Sequences',
    16: 'Aignable Sequence Segment',
    18: 'Sequence Trace',
}


def b2i(bits):
    return int(binascii.hexlify(bits), 16)


def parse_segments(handle):
    data = {'meta': {}, 'aln_seq': [], 'seq_trace': []}
    data['meta']['filetype'] = 'unknown' if binascii.hexlify(handle.read(2)) == '0000' else 'DNA'
    data['meta']['export_version_number'] = b2i(handle.read(2))
    data['meta']['import_version_number'] = b2i(handle.read(2))

    while True:
        segment_type_peek = handle.read(1)
        if segment_type_peek == '':
            break

        try:
            (key, keydata) = parse_segment(handle, segment_type_peek)
            if key not in ('aln_seq, seq_trace'):
                data[key] = keydata
            else:
                if key == 'aln_seq':
                    (seq_trace, ztr_data) = parse_segment(
                        StringIO.StringIO(keydata['data'][1:]),
                        keydata['data'][0]
                    )
                    if seq_trace != 'seq_trace':
                        raise Warning("Unexpected data...")
                    data['aln_seq'].append(ztr_data)
        except:
            pass
    return data


def parse_segment(handle, segment_type_peek):
    segment_type = b2i(segment_type_peek)
    segment_size = b2i(handle.read(4))
    segment_data = handle.read(segment_size)

    if segment_type == 0:
        dna_sequence = segment_data[1:]
        dna_meta = b2i(segment_data[0])
        return 'dna', {
            'sequence': dna_sequence,
            'topology': 'circular' if dna_meta & 0b1 else 'linear',
            'strandedness': 'double' if (dna_meta & 0b10 == 0b10) else 'single',
            'Dam methylated': dna_meta & 0b100 == 0b100,
            'Dcm methylated': dna_meta & 0b1000 == 0b1000,
            'EcoK1 methylated': dna_meta & 0b10000 == 0b10000,
        }
    elif segment_type == 8:
        pass
        # data['sequence_properties'] = segment_data
    elif segment_type == 10:
        # data['features'] = segment_data
        pass
    elif segment_type == 5:
        return 'primers', segment_data
    elif segment_type == 6:
        return 'notes', segment_data
    elif segment_type == 17:
        return 'aln_seqs', segment_data
    elif segment_type == 16:
        segment_seq_id = segment_data[0:4]
        segment_seq_data = segment_data[4:]
        return 'aln_seq', {
            'id': segment_seq_id,
            'data': segment_seq_data
        }
    elif segment_type == 18:
        return 'seq_trace', {
            'data': segment_data,
        }


def parse_file(path, outdir):
    handle = open(path, 'rb')

    ch = handle.read(1)
    assert ch == '\x09'

    ch = binascii.hexlify(handle.read(4))
    assert ch == '0000000e'

    ch = handle.read(8)
    assert ch == 'SnapGene'

    data = parse_segments(handle)
    # import pprint; pprint.pprint(data)
    gene_name = os.path.basename(path)[:-4]
    with open(os.path.join(outdir, gene_name + '.fasta'), 'w') as outfile:
        outfile.write('>' + gene_name.replace(' ', '_') + '\n')
        outfile.write(data['dna']['sequence'])
    handle.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='convert .dna to .fasta')
    parser.add_argument('dna', type=str, help='dna file')
    parser.add_argument('outdir', type=str, help='output directory')
    args = parser.parse_args()

    parse_file(args.dna, args.outdir)
