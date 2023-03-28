#!/usr/bin/env python
import pandas as pd
from Bio import SeqIO
import sys
import argparse


def read_gff(filepath):
    """
    Reads gene annotations in GFF format, only keep those that were annotated as isoforms, i.e., dual-coding
    :param filepath: str, path to file
    :return: pd.DataFrame, gene annotations
    """
    df = pd.read_csv(filepath, sep='\t',
                     header=None, comment='#',
                     names=['sequence', 'source', 'feature', 'start',
                            'end', 'logodd', 'strand', 'phase', 'attributes'])
    df['gene_id'] = [attr.split(';')[0] for attr in df.attributes.values]
    df['isoform'] = [True if '.' in gid else False for gid in df.gene_id.values]
    df = df[df.isoform] 
    df.reset_index(drop=True, inplace=True)
    return df


def read_sequence(filepath):
    """
    Parse genomic sequence stored in Fasta format
    :param filepath: str, path to fasta file
    :return: dict, genomic sequences
    """
    records = SeqIO.to_dict(SeqIO.parse(filepath, 'fasta'))
    return records


def extract_sequences_following_stop(df, sequence):
    """
    Extract 30 bp following a stop codon
    :param df: pd.DataFrame, gene annotations
    :param sequence: dict, corresponding genomic sequences
    :return: list, list, list, sequences following stop codon predicted by both code models, code model 11,
                               alternative code model with stop codon reassignment
    """
    sequences_11 = []
    sequences_alt = []
    sequences_all = []
    for i in range(df.shape[0]):
        start = df.loc[i, 'start']
        end = df.loc[i, 'end']
        strand = df.loc[i, 'strand']
        contig = sequence[df.loc[i, 'sequence']]
        source = df.loc[i, 'source']
        if strand == '+':
            seq = contig[end -3: end+30]
        elif strand == '-':
            seq = contig[start -31: start +2].reverse_complement()
        if len(seq) < 33:
           continue
        sequences_all.append(seq)
        if source.endswith('11'):
            sequences_11.append(seq)
        else:
            sequences_alt.append(seq)
    return sequences_all, sequences_11, sequences_alt


def write_sequences_to_file(sequences, output_file):
    """
    Writes genomic sequences to a txt file
    :param sequences: list, DNA sequences
    :param output_file: str, path to output file
    """
    with open(output_file, 'w') as out:
        for seq in sequences:
            out.write(str(seq.seq) + '\n')
    out.close()    


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--genomes', help='File specifying one genome per line that should be included in analysis')
    parser.add_argument('--mgcod_results', help='Path to directory containing MgCod annotations')
    args = parser.parse_args()
    genomes  = []
    gffs = []
    with open(args.genomes, 'r') as fp:
        for line in fp:
            genomes.append(line.strip())
            gffs.append(args.mgcod_results + line.strip().split('/')[-1].replace('.fna', '.gff')) 
    fp.close()
    sequences_11 = []
    sequences_alt = []
    sequences_all = []

    for gff, genome in zip(gffs, genomes):
        df = read_gff(gff)
        records = read_sequence(genome)
        seq_all, seq_11, seq_alt = extract_sequences_following_stop(df, records)
        sequences_all.extend(seq_all)
        sequences_11.extend(seq_11)
        sequences_alt.extend(seq_alt)
    write_sequences_to_file(sequences_all, 'sequences_following_stop_all.fasta')
    write_sequences_to_file(sequences_11, 'sequences_following_stop_11.fasta')
    write_sequences_to_file(sequences_alt, 'sequences_following_stop_alt.fasta')


if __name__ == '__main__':
    main(sys.argv[1:])
