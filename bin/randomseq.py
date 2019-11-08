#! /usr/bin/env python3 
import sys
import random as rd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    'length', metavar='N', type=int,
    help='length of the random sequence'
)

parser.add_argument(
    'numbers', metavar='N', type=int,
    help='numbers of sequences'
)

args = vars(parser.parse_args())

rd.seed(rd.randint(1000, 1000000))

for i in range(args['numbers']):
    sys.stdout.write(''.join(rd.choices(['A', 'T', 'G', 'C'], k=args['length'])) + '\n')

################################################################################
