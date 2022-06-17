#!/usr/bin/env python
"""
Script to simulate contigs with multiple genetic codes
"""
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import subprocess
import sys
import argparse


def read_genome(genome, coords, genetic_code):
    """
    Read template genome, and corresponding predictions
    :param genome: str, path to template genome
    :param path_to_predictions: str, path to corresponding mgm predictions
    :param genetic_code: int, genetic code from which to parse mgm predictions
    :return: Bio.Seq, pd.DataFrame, genomic sequences and corresponding gene annotation
    """
    records = list(SeqIO.parse(genome, 'fasta'))
    # load the predicted coordinates of corresponding genome
    cwd = sys.path[0]
    cwd = '/'.join(cwd.split('/')[:-1])
    if not os.path.isfile(coords):
            mgm = [f'{cwd}/lib/gmhmmp2', '-M', f'{cwd}/lib/mgm_{genetic_code}.mod', '-s',
                      genome, '-f', 'gff', '-o', coords]
            call_mgm = subprocess.check_call(mgm)
    coords = pd.read_csv(coords, sep='\t', header=None, comment='#', 
                            names=['seq_id', 'source', 'feature', 'start', 'end', 'logodd', 'strand', 'phase', 'attributes'])
    coords['length'] = coords['end'] - coords['start'] + 1
    coords['seq_id'] = coords['seq_id'].astype(str)
    if len(records) == 1:
        seqs = records[0].seq
    else:
        # if genome is assembled concatenate all the reads                                                                                                                                              
        seqs = records[0].seq
        for i in range(1, len(records)):
            seqs += records[i].seq

        # update the coordinates by adding length of all previous reads to obtain the correct positions in the concatenated sequence
        to_add = 0

        # iterate over all the ids of the reads (sometimes they are short so that there are no genes --> cannot use the dataframe)
        for i, seq_id in enumerate([seq.id for seq in records]):
            if i == 0:
                continue
            else:
                to_add += len(records[i - 1].seq)
                current_inds = coords[coords.seq_id == seq_id].index.values
                coords.loc[current_inds, 'start'] += to_add
                coords.loc[current_inds, 'end'] += to_add

    coords = coords.reset_index()
    seqs = seqs
    return seqs, coords


def simulate_switch_point(genome_11, genome_4, mgm_predictions_11, mgm_predictions_4, path_to_genomes,
                          path_to_annotations, output_name):
    """
    Simulate a genome with multiple genetic codes by merging random segments from two genomes with different segments
    :param genome_11: str, genome with standard genetic code
    :param genome_4: str, genome with genetic code 4
    :param mgm_predictions_11: str, mgm code 11 predictions corresponding to genome_11
    :param mgm_predictions_4: str, code 4 predictions corresponding to genome_4
    :param path_to_genomes: str, directory containing genomes
    :param path_to_annotations: str, directory containing mgm predictions
    :param output_name: str, path to save simulated genome to
    """
    # select random switch coordinates --> shorter segment must be at least 25kb
    splitting_point = np.random.randint(25000, 175000)
    # select random genome to start with
    start_with_11 = np.random.randint(2, size=1)[0]
    # read template genomes and predictions
    name_11 = genome_11.split('/')[-1].split('.fna')[0]
    name_4 = genome_4.split('/')[-1].split('.fna')[0]
    path_coords_11 = f"{mgm_predictions_11}/mgm_11_{name_11}.gff"
    path_coords_4 = f"{mgm_predictions_4}/mgm_4_{name_4}.gff"
    seq_11, coords_11 = read_genome(genome_11, path_coords_11, 11)
    seq_4, coords_4 = read_genome(genome_4, path_coords_4, 4)
    if start_with_11 == 1:
        coords1 = coords_11
        coords2 = coords_4
        seq1 = seq_11
        seq2 = seq_4
        label1 = 11
        label2 = 4
        name1 = name_11
        name2 = name_4
    else:
        coords1 = coords_4
        coords2 = coords_11
        seq1 = seq_4
        seq2 = seq_11
        label1 = 4
        label2 = 11
        name1 = name_4
        name2 = name_11
    # select coordinates of random segments to be spliced out from template genomes
    start_1 = np.random.randint(0, len(seq1) - splitting_point - 1)
    start_2 = np.random.randint(0, len(seq2) - (200000 - splitting_point) - 1)
    # ensure sequences are separated by intergenic region
    margin1 = 0
    margin2 = 0
    # ensure slice sides fall into intergenic regions
    if coords1[(coords1['start'] <= start_1) & (coords1['end'] < start_1)].index.values.shape[0] >= 1:
        ind = coords1[(coords1['start'] <= start_1) & (coords1['end'] < start_1)].index.values[-1]
        # move to intergenic region between current gene and next gene
        start_1 = int((coords1.loc[ind, 'end'] + coords1.loc[ind + 1, 'start']) / 2)
     
    if coords1[(coords1['start'] <= start_1 + splitting_point) &
               (coords1['end'] > splitting_point + start_1)].index.values.shape[0] >= 1:
        ind = coords1[(coords1['start'] <= start_1 + splitting_point) &
                      (coords1['end'] > splitting_point + start_1)].index.values[-1]
        # move to intergenic region between current gene or end of sequence
        try:
            margin1 = int((coords1.loc[ind, 'end'] + coords1.loc[ind + 1, 'start']) / 2)  - start_1 - splitting_point - 1
        except KeyError:
            margin1 = len(seq1) - start_1 - splitting_point - 1
    if coords2[(coords2['start'] <= start_2) & (coords2['end'] < start_2)].index.values.shape[0] >= 1:
        ind = coords2[(coords2['start'] <= start_2) & (coords2['end'] < start_2)].index.values[0]
        # move to previous intergenic region or start of sequence
        try:
            start_2 = int((coords2.loc[ind, 'start'] + coords2.loc[ind - 1, 'end']) / 2)
        except KeyError:
            start_2 = 0
    if coords2[(coords2['start'] <= start_2 + (200000 - splitting_point)) &
               (coords2['end'] > start_2 + (200000 - splitting_point))].index.values.shape[0] >= 1:
        ind = coords2[(coords2['start'] <= start_2 + (200000 - splitting_point)) &
                      (coords2['end'] > start_2 + (200000 - splitting_point))].index.values[0]
        # move to previous intergenic region
        margin2 = int((coords2.loc[ind, 'start'] + coords2.loc[ind - 1, 'end']) / 2) - start_2 - (200000 - splitting_point) - 1

    # simulated genome with multiple genetic codes
    mixed_seq = seq1[start_1: start_1 + splitting_point + margin1] + seq2[start_2: start_2 + (200000 - splitting_point) + margin2]
    switch_zone = ((coords1['end'][coords1['end'] < start_1 + margin1 + splitting_point].values[-1] - start_1),
                   (coords2['start'][coords2['start'] >= start_2].values[0] - start_2 + margin1 + splitting_point))
    with open(f"{path_to_genomes}{output_name}.fna", 'w+') as f:
        f.writelines(f">{name1}, {name2}, {label1} to {label2} between {switch_zone}")
        f.writelines('\n')
        f.writelines(mixed_seq)
    f.close()
    with open(f"{path_to_annotations}switch_points.txt", 'a+') as f:
        f.writelines(f"{output_name}\t[{label1}, {label2}]\t{switch_zone}")
        f.writelines('\n')
    f.close()
    coords2['index'] = coords2['index'].values + (coords1.shape[0] + 1)
    coords1 = coords1.append(coords2)
    coords1.reset_index(drop=True, inplace=True)
    coords1.to_csv(f"{path_to_annotations}{output_name}.csv", sep='\t', header=True, index=False)


def main(argv):
    parser = argparse.ArgumentParser(description="Simulate contigs with multiple genetic codes")
    parser.add_argument('-d1', '--directory_genomes_11', help='Directory to genomes with standard genetic code',
                        required=True)
    parser.add_argument('-d2', '--directory_genomes_4', help='Directory to genome with genetic code 4', required=True)
    parser.add_argument('-p1', '--mgm_predictions_11', help='Directory to MGM predictions of genomes with standard genetic code',
                        required=True)
    parser.add_argument('-p2', '--mgm_predictions_4', help='Directory to MGM predictions of genomes with genetic code 4',
                        required=True)
    parser.add_argument('-g', '--genomes', help='Output directory for simulated genomes', required=True)
    parser.add_argument('-a', '--annotations', help='Output directory for annotations', required=True)
    args = parser.parse_args()

    genomes_11 = [os.path.join(args.directory_genomes_11, f) for f in os.listdir(args.directory_genomes_11)
                  if os.path.isfile(os.path.join(args.directory_genomes_11, f))]
    genomes_4 = [os.path.join(args.directory_genomes_4, f) for f in os.listdir(args.directory_genomes_4)
                 if os.path.isfile(os.path.join(args.directory_genomes_4, f))]

    nr_genomes_11 = len(genomes_11)
    nr_genomes_4 = len(genomes_4)

    if not os.path.isdir(args.genomes):
        os.makedirs(args.genomes)
    if not os.path.isdir(args.annotations):
        os.makedirs(args.annotations)
    mgm_predictions_11 = args.mgm_predictions_11
    if not mgm_predictions_11.endswith('/'):
        mgm_predictions_11 += '/'
    
    mgm_predictions_4 = args.mgm_predictions_4
    if not mgm_predictions_4.endswith('/'):
        mgm_predictions_4 += '/'

    for i in range(0, 1000, 1):
        current_genome_11 = genomes_11[np.random.randint(0, nr_genomes_11)]
        current_genome_4 = genomes_4[np.random.randint(0, nr_genomes_4)]
        output_name = f"switch_point_{i}"
        simulate_switch_point(current_genome_11, current_genome_4, mgm_predictions_11, mgm_predictions_4, args.genomes, args.annotations, output_name)


if __name__ == '__main__':
    main(sys.argv[1:])
