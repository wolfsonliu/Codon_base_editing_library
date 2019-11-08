#! /usr/bin/env python3
####################
# First A site  |  sgRNA length
#       19      |      21
#       18      |      20
#       17      |      19
#       16      |      19
#       15      |      19
#       14      |      20
#       13      |      20
####################

import sys
import re
import numpy as np
import pandas as pd
from Bio import SeqIO
import argparse

####################

parser = argparse.ArgumentParser(prog='get_sgrna_pos')

parser.add_argument(
    '-r', '--reference', required=True,
    help='reference genome fasta file.'
)

parser.add_argument(
    '-i', '--input', required=True,
    help="""input tab separated data file with columns: [
    seqname, strand, aa_genome_pos_start, aa_genome_pos_end, aa_genome_seq].
    with header"""
)

parser.add_argument(
    '-a', '--aa-sgrna', required=True, nargs='*', type=int,
    help="""First A site on the sgRNA, counting from PAM.
    For instance: -a 14 15 16"""
)

parser.add_argument(
    '-l', '--spacer-length', required=True, type=int,
    help="""sgRNA spacer length"""
)

parser.add_argument(
    '-o', '--output', help='output file for sgrna position, add columns  [Gene Name, Amino Acid (one character), Protein Position].'
)

args = vars(parser.parse_args())

####################

def complement(seq):
    # calculate the complement sequence string
    # input a string
    # return a string
    iupac = {
        "A": "T", "G": "C", "C": "G", "T": "A", "Y": "R", "R": "Y",
        "W": "W", "S": "S", "K": "M", "M": "K", "D": "H", "V": "B",
        "H": "D", "B": "V", "N": "N", "X": "X", "-": "-",
        "a": "t", "g": "c", "c": "g", "t": "a", "y": "r", "r": "y",
        "w": "w", "s": "s", "k": "m", "m": "k", "d": "h", "v": "b",
        "h": "d", "b": "v", "n": "n", "x": "x", "-": "-",
    }
    return ''.join(iupac[i] for i in seq)


def reverse_complement(seq):
    # calculate the reverse complement sequence string
    # input a string
    # return a string
    return complement(seq)[::-1]

def gc_content(seq):
    # calculate the sequence GC content
    # input a string
    # return a number
    return len([i for i in seq if i in ['C','G','S']])/len(seq)

####################
__console__ = sys.stdout
if 'output' in args.keys():
    sys.stdout = open(args['output'], 'w')

colnames = [
    'seqname', 'strand', 'aa_genome_pos_start', 'aa_genome_pos_end', 'aa_genome_seq',
    'sgrna_genome_start', 'sgrna_genome_end',
    'sgrna_genome_seq', 'genome_pam', 'sgrna_seq', 'pam', 'codon_1stnt_pos', 'gc'
]

data = pd.read_table(args['input'], header=None)

data.columns = ['seqname', 'strand', 'aa_genome_pos_start', 'aa_genome_pos_end', 'aa_genome_seq']

result = list()

for a in args['aa_sgrna']:
    result.append(
        data.loc[data['strand'] == '+'].assign(
            sgrna_genome_start=lambda x: x.aa_genome_pos_start + a - args['spacer_length'],
            sgrna_genome_end=lambda x: x.aa_genome_pos_start + a - 1,
            codon_1stnt_pos=lambda x: a
        )
    )
    result.append(
        data.loc[data['strand'] == '-'].assign(
            sgrna_genome_start=lambda x: x.aa_genome_pos_end - a + 1,
            sgrna_genome_end=lambda x: x.aa_genome_pos_end - a + args['spacer_length'],
            codon_1stnt_pos=lambda x: a
        )
    )

del data
output = pd.concat(result, axis=0).reset_index(drop=True)

vseqname = [''.join(['chr', str(i)]) for i in list(range(1, 23)) + ['X', 'Y']]
fa = dict()
for record in SeqIO.parse(args["reference"], "fasta"):
    if record.id in vseqname:
        fa[record.id] = record

output.loc[:,'sgrna_genome_seq'] = output.apply(
    lambda x: str(
        fa[x['seqname']].seq[x['sgrna_genome_start'] - 1: x['sgrna_genome_end']]
    ).upper(),
    axis=1
)
output.loc[:,'genome_pam'] = output.apply(
    lambda x: str(
        fa[x['seqname']].seq[x['sgrna_genome_end']: x['sgrna_genome_end'] + 3]
    ).upper()
    if x['strand'] == '+'
    else str(
        fa[x['seqname']].seq[x['sgrna_genome_start']-1-3: x['sgrna_genome_start']-1]
    ).upper(),
    axis=1
)
output.loc[:,'sgrna_seq'] = output.apply(
    lambda x: x['sgrna_genome_seq']
    if x['strand'] == '+' else reverse_complement(x['sgrna_genome_seq']) ,
    axis=1
)
output.loc[:,'pam'] = output.apply(
    lambda x: x['genome_pam']
    if x['strand'] == '+' else reverse_complement(x['genome_pam']) ,
    axis=1
)

output.loc[:,'gc'] = output.apply(lambda x: gc_content(x['sgrna_seq']), axis=1)

output.loc[:,'keep'] = (output['pam'].str[1:] == 'GG')

output[colnames].loc[output['keep']].to_csv(
    sys.stdout, index=False, quoting=0, sep='\t'
)
sys.stdout = __console__
################################################################################
