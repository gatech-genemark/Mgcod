#!/usr/bin/env python
"""
Script to split sequences of annotated coding sequences by genetic code
"""
from Bio import SeqIO
import argparse
import sys


def read_coding_seq(fasta_file):
    """
    Read sequences of annotated coding sequences
    :param fasta_file: str, path to file with annotated coding sequences
    :return: list, Bio.SeqRecords of annotated coding regions
    """
    records = list(SeqIO.parse(fasta_file, 'fasta'))
    return records


def split_by_genetic_code(seqs):
    """
    Split sequences by genetic code which is inferred from SeqRecords
    :param seqs: list, Bio.SeqRecords of annotated coding regions
    :return: dict, mapping of genetic code to sequences
    """
    genetic_codes = {}
    for seq in seqs:
        code = seq.description.split(',')[1].split('_')
        if len(code) == 2:
            code = code[1]
        elif len(code) == 3:
            code = code[1] + "_" + code[2]
        if not code in genetic_codes:
            genetic_codes[code] = [seq]
        else:
            genetic_codes[code].append(seq)
    return genetic_codes


def write_coding_seq(seqs_by_code, seq_type, output_base):
    """
    Write sequences with different genetic codes to different files
    :param seqs_by_code: dict, mapping of genetic code to sequences
    :param seq_type: str, indicating sequence type either nt for nucleotide or aa for amino acid
    :param output_base: str, directory where to save results
    """
    for code, seqs in seqs_by_code.items():
        if seq_type == 'aa':
            SeqIO.write(seqs, f'{output_base}/proteins_{code}.faa', 'fasta')
        elif seq_type == 'nt':
            SeqIO.write(seqs, f'{output_base}/proteins_{code}.fna', 'fasta')


def main(argv):
    parser = argparse.ArgumentParser(description="Split sequences of annotated coding sequences by genetic code")
    parser.add_argument('--input_fasta', help='Fasta file with proteins of different genetic code outputed by AGCA', required=True)
    parser.add_argument('--seq_type', required=True)
    args = parser.parse_args()
    output_base = '/'.join(args.input_fasta.split('/')[:-1])
    records = read_coding_seq(args.input_fasta)
    seqs_by_code = split_by_genetic_code(records)
    write_coding_seq(seqs_by_code, args.seq_type, output_base)


if __name__ == '__main__':
    main(sys.argv[1:])
