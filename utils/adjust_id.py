#!/usr/bin/env python
"""
Script to update sequence ids in multi-fasta file
"""
import argparse
import sys
from Bio import SeqIO


def main(argv):
    parser = argparse.ArgumentParser(description="Update sequence ids in multi-fasta file")
    parser.add_argument('input_file')
    args = parser.parse_args()
    records = list(SeqIO.parse(args.input_file, 'fasta'))
    for i, rec in enumerate(records):
        rec.id = str(i)
    SeqIO.write(records, args.input_file, 'fasta')


if __name__ == '__main__':
    main(sys.argv[1:])
