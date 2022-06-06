#!/usr/bin/env python
"""
Script to simulate genomes with stop codon reassignments
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
import sys
import argparse
import pandas as pd
import glob
import os
import subprocess


def read_data(genome, path_to_predictions, genetic_code):
    """
    Read template genome, and corresponding predictions
    :param genome: str, path to template genome
    :param path_to_predictions: str, path to corresponding mgm predictions
    :param genetic_code: int, genetic code from which to parse mgm predictions
    :return: Bio.Seq, pd.DataFrame, genomic sequences and corresponding gene annotation
    """
    records = list(SeqIO.parse(genome, 'fasta'))
    if not os.path.isfile(path_to_predictions):
        cwd = sys.path[0]
        cwd = '/'.join(cwd.split('/')[:-1])
        mgm = [f'{cwd}/lib/gmhmmp2', '-M', f'{cwd}/lib/mgm_{genetic_code}.mod', '-s', genome, '-f', 'gff', '-o', path_to_predictions]
        call_mgm = subprocess.check_call(mgm)
    coords = pd.read_csv(path_to_predictions, sep='\t', header=None, comment='#',
                         names=['seq_id', 'source', 'feature', 'start',
                         'end', 'logodd', 'strand', 'phase', 'attributes'])
    if len(records) == 1:
        seqs = records[0].seq
    else:
        # if genome is assembled concatenate all the reads
        seqs = records[0].seq
        for i in range(1, len(records)):
            seqs += records[i].seq
        # update the coordinates by adding length of all previous reads to obtain the correct positions in the concatenated sequence
        to_add = 0
        #iterate over all the ids of the reads (sometimes they are short so that there are no genes --> cannot use the dataframe)
        for i, seq_id in enumerate([seq.id for seq in records]):
            if i == 0:
                continue
            else:
                to_add += len(records[i - 1].seq)
                current_inds = coords[coords.seq_id == seq_id].index.values
                coords.loc[current_inds, 'start'] += to_add
                coords.loc[current_inds, 'end'] += to_add
    coords = coords.reset_index()
    return seqs, coords


def infer_reverse_codon(codon):
    """
    Infer codon on complementary strand
    :param codon: str, input codon
    :return: str, complementary codon
    """
    reverse_codon = ''
    for i in range(2, -1, -1):
        if codon[i] == 'A':
            reverse_codon += 'T'
        elif codon[i] =='T':
            reverse_codon += 'A'
        elif codon[i] =='G':
            reverse_codon += 'C'
        elif codon[i] =='C':
            reverse_codon += 'G'
        else:
            reverse_codon += 'N'
    return reverse_codon    


def substitute_stop_codon(genome, predictions, codon_to_substitute, to_substitute_with):
    """
    Substitute stop codon that is simulated to be reassigned with one of the other stop codons
    :param genome: Bio.Seq, template genomic sequence
    :param predictions: pd.DataFrame, gene annotation
    :param codon_to_substitute: str, stop codon to replace
    :param to_substitute_with: list, remaining stop codons
    :return: Bio.Seq, modified template genomic sequences with out codon regions ending stop codon that is
                      simulated to be reassignd
    """
    # convert to mutable object to allow editiing
    genome = genome.tomutable()
    # get starts ends and strands values
    starts, ends, strands = predictions.start.values, predictions.end.values, predictions.strand.values
    # infer reverse codons
    reverse_codon_to_substitute = infer_reverse_codon(codon_to_substitute)
    reverse_to_substitute_with = [infer_reverse_codon(codon) for codon in to_substitute_with]
    for start, end, strand in zip(starts, ends, strands):
        # replace stop codon of all genes ending reassigned codon
        if strand == "+":
            if len(to_substitute_with) == 2:
                substitution = to_substitute_with[np.random.choice([0, 1], p=[0.35, 0.65])]
            else:
                substitution = to_substitute_with[0]
            if genome[end - 3: end] == codon_to_substitute:
                genome[end - 3:end] = substitution
        if strand == "-":
            if len(reverse_to_substitute_with) == 2:
                substitution = reverse_to_substitute_with[np.random.choice([0, 1], p=[0.35, 0.65])]
            else:
                substitution = reverse_to_substitute_with[0]
            if genome[start - 1: start + 2] == reverse_codon_to_substitute:
                genome[start - 1: start + 2] = substitution
    return genome


def substitute_coding_codon(genome, predictions, codons_to_substitute, to_substitute_with, probability):
    """
    Insert stop codon that is simulated to be reassigned as sense codon into annotated codon regions
    :param genome: Bio.Seq, template genomic sequence
    :param predictions: pd.DataFrame, gene annotation
    :param codons_to_substitute: list, codons coding for a particular amino acid for which the reassigned stop codon
                                       is supposed to code for
    :param to_substitute_with: str, stop codon that is simulated to be reassigned
    :param probability: float, probability with that a codon is replace with a former stop codon
    :return: Bio.Seq, modified template genomic sequence
    """
    starts, ends, strands = predictions.start.values, predictions.end.values, predictions.strand.values
    # infer reverse codons for all codons coding for amino acid to which the stop codon is reassigned
    reverse_codons_to_substitute = []
    for codon in codons_to_substitute:
        reverse_codons_to_substitute.append(infer_reverse_codon(codon))
    reverse_to_substitute_with = infer_reverse_codon(to_substitute_with)
    # get variable telling if next appearance should be changed or not --> determined by probability
    #substitute = np.random.choice([0, 1], p=[1 - probability, probability])
    for start, end, strand in zip(starts, ends, strands):
        i = start - 1
        substitute = np.random.choice([0, 1], p=[1 - probability, probability])
        while i + 3 <= end:
            if strand == "+":
                if genome[i: i + 3] in codons_to_substitute:
                    if substitute == 1:
                        genome[i: i + 3] = to_substitute_with
                    # get variable telling if next appearance should be changed or not --> determined by probability                
                    substitute = np.random.choice([0, 1], p=[1 - probability, probability])
            if strand == "-":
                if genome[i: i + 3] in reverse_codons_to_substitute:
                    if substitute == 1:
                        genome[i: i + 3] = reverse_to_substitute_with
                    # get variable telling if next appearance should be changed or not --> determined by probability
                    substitute = np.random.choice([0, 1], p=[1 - probability, probability])

            i += 3
    return genome
        

def main(argv):
    parser = argparse.ArgumentParser(description="Simulate genomes with stop codon reassignments")
    parser.add_argument("-d", '--directory', help="Directory to input genomes", required=True)
    parser.add_argument("-p", '--mgm_predictions', help="Directory to mgm predictions", required=True)
    parser.add_argument('-g', '--genetic_code', help="Original genetic code", required=True)
    parser.add_argument('-s', '--stop_codon_to_replace', help='Stop codon to replace')
    parser.add_argument('-r', '--stop_codons_to_replace_with', nargs='+', help='Stop codons to use')
    parser.add_argument('-c', '--coding_codons_to_replace', nargs='+', help='Coding codons to replace')
    parser.add_argument('--prob', help='Probablility with that to replace coding codons. [0.8]', default=0.8, type=float)
    parser.add_argument('-o', '--output_dir', help='Directory where to save simulated genomes')
    args = parser.parse_args()
    path_to_genomes = args.directory
    if not path_to_genomes.endswith('/'):
        path_to_genomes += '/'
    path_to_predictions = args.mgm_predictions
    if not path_to_predictions.endswith('/'):
        path_to_predictions += '/'
    path_to_outputs = args.output_dir
    if not path_to_outputs.endswith('/'):
        path_to_outputs += '/'
    genomes = glob.glob('{}*.fna'.format(path_to_genomes))
    predictions = ['{}mgm_{}_{}.gff'.format(path_to_predictions, args.genetic_code, genome.split('/')[-1].split('.fna')[0]) for genome in genomes]
    outputs = ['{}{}'.format(path_to_outputs, genome.split('/')[-1]) for genome in genomes]    

    stop_codon_to_substitute = args.stop_codon_to_replace
    to_substitute_with = args.stop_codons_to_replace_with
    coding_codon_to_substitute = args.coding_codons_to_replace
    probability = args.prob

    for genome, prediction, output in zip(genomes, predictions, outputs):
        seq, coords = read_data(genome, prediction, args.genetic_code)
        genome = substitute_stop_codon(seq, coords, stop_codon_to_substitute, to_substitute_with)
        genome = substitute_coding_codon(genome, coords, coding_codon_to_substitute, stop_codon_to_substitute, probability)
        genome = SeqRecord(genome)
        #pdb.set_trace()
        SeqIO.write(genome, output, "fasta")


if __name__=="__main__":
    main(sys.argv[1:])
