#!/usr/bin/env python
"""
Script to split genome into equal sized contigs
"""
import sys
import argparse
from Bio import SeqIO
import os


def split_multifasta(genome, contig_size, output_dir, suffix):
    """
    Split genome into equal sized contigs
    :param genome: str, path to genome
    :param contig_size: int, contig size
    :param output_dir: str, directory where to write contigs to
    :param suffix: str, suffix of input file (e.g., fna, fasta, fa)
    """
    # load the sequence
    genome_name = genome.split('/')[-1].split(f".{suffix}")[0]
    records = list(SeqIO.parse(genome, 'fasta'))
    header = records[0].description
    if len(records) == 1:
        seqs = records[0].seq
    else:
        # if genome is assembled concatenate all the reads                                                                                                                                              
        seqs = records[0].seq
        for i in range(1, len(records)):
            seqs += records[i].seq
    contigs = [seqs[contig_size * i: contig_size * i + contig_size] for i in range(0, len(seqs) // contig_size)]
    for i, contig in enumerate(contigs):
        with open(f"{output_dir}{genome_name}_{i}.fna", 'w+') as c:
            c.writelines(f">{header}")
            c.writelines('\n')
            c.writelines(contig)
        c.close()


def main(argv):
    parser = argparse.ArgumentParser(description="Split genome into equal sized contigs")
    parser.add_argument('-i', '--input', help='input genomes which is meant to split', required=True, dest='genome')
    parser.add_argument('-c', '--contig-size', help='size of contigs to split genome in', type=int, required=True, dest='contig_size')
    args = parser.parse_args()
    args = vars(args)
    genome = args['genome']
    contig_size = args['contig_size']
    suffix = genome.split('.')[-1]
    output_dir = f"{'/'.join(genome.split('/')[:-1])}/contigs_{contig_size}/"
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir) 
    split_multifasta(genome, contig_size, output_dir, suffix)


if __name__=='__main__':
    main(sys.argv[1:])
