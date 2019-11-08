#! /usr/bin/env python3
import sys
import re
import numpy as np
import pandas as pd
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(prog='get_sgrna_pos')

parser.add_argument(
    '-r', '--reference', required=True,
    help='reference genome fasta file.'
)

parser.add_argument(
    '-i', '--input', required=True,
    help="""input tab separated data file with columns: [
       name, strand, seqname, start, seq, quality, n, mismatch
    ]. (No title) Start should be 0 based."""
)

parser.add_argument(
    '-o', '--output',
    help="""output file with columns: [seqname, start, end, strand, ref, seq].
    (No title) for + strand, ref is equal to seq, for the - strand,
    seq is the reverse complementary of ref"""
)

args = vars(parser.parse_args())


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

####################
__console__ = sys.stdout
if 'output' in args.keys():
    sys.stdout = open(args['output'], 'w')

colnames = [
    'name', 'strand', 'seqname', 'start', 'seq', 'quality', 'n', 'mismatch'
]

data = pd.read_table(args['input'], header=None)
data.columns = colnames

# read genomes
vseqname = [''.join(['chr', str(i)]) for i in list(range(1, 23)) + ['X', 'Y']]
fa = dict()
for record in SeqIO.parse(args["reference"], "fasta"):
    fa[record.id] = record

data.loc[:,'pstart'] = data.apply(
    lambda x: x['start'] + len(x['name']) if x['strand'] == '+' else x['start'] - 3,
    axis=1
)

data.loc[:,'pend'] = data.apply(
    lambda x: x['start'] + len(x['name']) + 3 if x['strand'] == '+' else x['start'],
    axis=1
)

data.loc[:,'pam'] = data.apply(
    lambda x: str(fa[x['seqname']].seq[x['pstart']: x['pend']]).upper()
    if x['strand'] == '+' else reverse_complement(
            str(fa[x['seqname']].seq[x['pstart']: x['pend']]).upper()
    ),
    axis=1
)

data[colnames + ['pstart', 'pend', 'pam']].to_csv(
    sys.stdout, index=False, header=False, sep='\t'
)

sys.stdout = __console__
################################################################################
