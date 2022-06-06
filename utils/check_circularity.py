#!/usr/bin/env python
"""
Script to check if a contig is circular based on overlap between the ends
"""
from Bio import SeqIO
import sys
import argparse


def check_circularity(queries, output, min_overlap=25, max_allowed_overlap=200):
    """
    Determines circularity of query contigs based on overlap  between the ends
    :param queries: list, contigs to check for circularity
    :param output: str, output file to write results to
    :param min_overlap: int, minimum overlap between ends
    :param max_allowed_overlap: int, maximum overlap between ends
    """
    with open(queries, 'r') as queries:
        with open(output, 'w+') as out:
            for query in queries:
                seq = list(SeqIO.parse(query.strip(), 'fasta'))[0].seq
                max_overlap = 0
                i = 1
                while i <= max_allowed_overlap:
                    if seq[:i] == seq[-i:]:
                        max_overlap = i
                    i += 1
                if max_overlap >= min_overlap:
                    out.write(f'{query.strip()}\tcircular\t{max_overlap}\n')
                else:
                    out.write(f'{query.strip()}\tnot-circular\t{max_overlap}\n')
        out.close()
    queries.close()


def main(argv):
    parser = argparse.ArgumentParser(description="Check if a contig is circular based on overlap between the ends")
    parser.add_argument('--queries', help='File with query genomes, one per line', required=True)
    parser.add_argument('--output', help='Output file', required=True)
    parser.add_argument('--min_overlap', help='Minimal overlap between ends. [25]', type=int, default=25)
    parser.add_argument('--max_overlap', help='Maximal overlap between ends. [200]', type=int, default=200)
    args = parser.parse_args()
    check_circularity(args.queries, args.output, args.min_overlap, args.max_overlap)


if __name__ == '__main__':
    main(sys.argv[1:])
