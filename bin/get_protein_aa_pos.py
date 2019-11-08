#! /usr/bin/env python3
import sys
import re
import argparse
from Bio import SeqIO

####################

parser = argparse.ArgumentParser(prog='get_protein_aa_pos')

parser.add_argument(
    '-f', '--fasta', required=True,
    help='input fasta file for protein AA sequences.'
)

parser.add_argument(
    '-a', required=True, help='output file for protein AA sequences.'
)

parser.add_argument(
    '-o', '--output',
    help='output file for protein AA sequences, [Gene Name, Amino Acid (one character), Protein Position].'
)

args = vars(parser.parse_args())

####################
__console__ = sys.stdout
if 'output' in args.keys():
    sys.stdout = open(args['output'], 'w')

aa = args['a'].upper()
aap = re.compile(aa)

proteins = dict()

for record in SeqIO.parse(args["fasta"], "fasta"):
    name = [x for x in record.description.split(' ') if x.find('GN') >= 0]
    record.name = name[0].replace('GN=','') if len(name) == 1 else ''
    splitid = record.id.split('|')
    proteins[record.name] = record
    [
        sys.stdout.write(
            '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(
                splitid[1], splitid[2], record.name, aa, pos
            )
        )
        for pos in [x.span()[1] for x in aap.finditer(str(record.seq))]
    ]

sys.stdout = __console__
################################################################################
