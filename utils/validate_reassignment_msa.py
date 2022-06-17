#!/usr/bin/env python
"""
Script to validate predicted stop codon reassignments using MSA
"""
import sys
import argparse
import pandas as pd
import glob
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import logging


def read_queries(input_file):
    """
    Parse paths to input genomes
    :param input_file: str, file containing paths to query genomes
    :return: list, query sequences
    """
    queries = []
    with open(input_file, 'r') as infile:
        for line in infile:
            queries.append(line.strip())
    infile.close()
    return queries


def predict_genes(input_files, output_dir, mgcod_results, mgm_results, mgcod):
    """
    Run Mgcod on target genomes
    :param input_files: list, paths to target genomes
    :param output_dir: str, directory where to save results
    :param mgcod_results: str, directory where to save Mgcod results
    :param mgm_results: str, directory where to save MGM results
    :param mgcod: str, path to Mgcod script
    :return: dictionary, mapping gene ids to gene descriptions
    """
    cwd = output_dir
    written_genes = 0
    gene_id_mapping = {}
    if os.path.isfile(f'{cwd}to_cluster.fna'):
        os.remove(f'{cwd}to_cluster.fna')
    for path_to_file in input_files:
        if not os.path.isfile(f"{mgcod_results}proteins_nt_{path_to_file.split('/')[-1][:-4]}.fasta"):
            subprocess.call(['python', mgcod, '-i', f'{path_to_file}',
                             '-o', mgcod_results, '-p', mgm_results, '-r', '--isoforms', '-NT', '-AA'])
        # prepare for clustering
        # read predictions
        records_nt = list(SeqIO.parse(f"{mgcod_results}proteins_nt_{path_to_file.split('/')[-1][:-4]}.fasta", 'fasta'))

        #records_capitalized = []
        for rec_nt in records_nt:
            rec_nt.seq = rec_nt.seq.upper()
            rec_nt.id = str(int(rec_nt.id) + written_genes)
            rec_nt.description = ' '.join(rec_nt.description.split(' ')[1:])
            gene_id_mapping[rec_nt.id] = rec_nt.description

        # write genes
        with open(f'{cwd}to_cluster.fna', 'a+') as output_handle:
            SeqIO.write(records_nt, output_handle, 'fasta')
        output_handle.close()

        written_genes += len(records_nt)
    return gene_id_mapping


def cluster_genes(gene_id_contig_mapping, output_dir, usearch):
    """
    Cluster nucleotide sequences of predicted coding regions using UCLUST
    :param gene_id_contig_mapping: dict, gene id to contig of origin mapping
    :param output_dir: str, directory where to save clustering results
    :param usearch: str, path to usearch binary
    :return: list, list paths to clusters of nucleotide sequences;
                        pd.DataFrame mapping sequences in cluster to contigs of origin
    """
    cwd = output_dir
    #if os.path.isfile(f'{cwd}clustered_genes.fna'):
    #    os.remove(f'{cwd}clustered_genes.fna')
    if not os.path.isdir(f'{cwd}cluster/'):
        os.makedirs(f'{cwd}cluster/')
    #else:
    previous_cluster = glob.glob(f'{cwd}cluster/*')
    for cluster in previous_cluster:
        os.remove(cluster)

    subprocess.call([usearch, '--cluster_fast', f'{cwd}to_cluster.fna', '-id', '0.8', '-sort', 'size', '-clusters', f'{cwd}cluster/'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) 
    clusters = glob.glob(f'{cwd}cluster/*')
    cluster_contig_mapping = pd.DataFrame(index=set(gene_id_contig_mapping.values()))
    cluster_for_msa = []
    for cluster in clusters:
        cluster_id = cluster.split('/')[-1]
        records = list(SeqIO.parse(cluster, 'fasta'))
        for rec in records:
            cluster_contig_mapping.loc[gene_id_contig_mapping[rec.id], cluster_id] = 1
        if len(records) > 1:
            cluster_for_msa.append(cluster)
    return set(cluster_for_msa), cluster_contig_mapping


def do_msa(cluster_nrs, output_dir, clustalo):
    """
    Perform MSA using ClustalOmega
    :param cluster_nrs: list, cluster ids
    :param output_dir: str, output directory
    :param clustalo: str, path to clustalo binary
    """
    cwd = output_dir
    if not os.path.isdir(f'{cwd}msa_results/'):
        os.makedirs(f'{cwd}msa_results/')
    else:
        previous_msa = glob.glob(f'{cwd}msa_results/*.clu')
        for msa in previous_msa:
            os.remove(msa)
    for cluster in cluster_nrs:
        clustalo_cl = ClustalOmegaCommandline(clustalo,
                                              infile=cluster, outfile=f"{cwd}msa_results/{cluster.split('/')[-1].split('.')[0]}.clu",
                                              outfmt='clustal', outputorder='input-order', force=True, threads=16)
        clustalo_cl()


def get_other_codons_at_reassigned_positions(records, stop_codon):
    """
    Get other codons to which putatively reassigned stop codon aligns in MSA
    :param records: list, paths to MSA results for all clusters
    :param stop_codon: str, target stop codon
    :return: list, codons to which putatively reassigned stop codon aligns in MSA
    """
    aligned_codons = []
    seen_positions = []
    for rec in records:
        length = len(rec.seq)
        i = 0
        while i <= length -4 and rec.seq[i:].find(stop_codon) != -1:
            tag = rec.seq[i:].find(stop_codon)
            if (i + tag) % 3 == 0 and tag not in seen_positions:
                for rec1 in records:
                    aligned_codons.append(str(rec1.seq[i + tag: i + tag + 3]))
                seen_positions.append(tag)
            i += tag + 3
    return aligned_codons


def capitalize_in_frame_stop(seq, stop_codon):
    """
    Formats sequence by capitalizing in frame stop codons
    :param seq: str, nucleotide sequence
    :param stop_codon: str, target stop codon
    :return: str, formatted nucleotide sequence
    """
    length = len(seq)
    i = 0
    while i <= length - 4 and seq[i:].find(stop_codon) != -1:
        tag = seq[i:].find(stop_codon)
        if (i + tag) % 3 == 0:
            seq = seq[: i + tag] + seq[i + tag: i + tag + 3].upper() + seq[i + tag + 3:]
        i += tag + 3
    return seq


def write_to_file(records, filename):
    """
    Format MSA  by capitalizing in frame stop codons
    :param records: list, MSA results
    :param filename: list, path of MSA result
    """
    filename = filename.replace('.clu', '_formatted.clu')
    with open(filename, 'w+') as handle:
        for rec in records:
            handle.write(rec.description + '\n')
        handle.write('\n\n')
        alignment_length = len(records[0].seq)
        i = 0
        while i < alignment_length:
            ref = np.array([n for n in str(records[0].seq[i:i+50])])
            matches = np.array(['*'] * len(ref))
            for rec in records:
                nuc_rec = np.array([n for n in str(rec.seq[i:i+50])])
                if np.all(nuc_rec == ref):
                    pass
                else:
                    matches[np.where(nuc_rec != ref)[0]] = ' '
                handle.write(rec.id +'\t' + str(rec.seq[i:i+50]) + '\n')
            handle.write('\t' + ''.join(matches.tolist()))
            handle.write('\n\n')
            i += 50
    handle.close()
        

def summarize_aligned_codons(aligned_codons, stop_codon):
    """
    Convert counts of number of times a stop codon aligned to another codon to amino acids
    :param aligned_codons: list, of aligned codons
    :param stop_codon: str, target stop codon
    :return: dict, counts of number of times a stop codon aligned a particular amino acid
    """
    table = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
             'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
             'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C', 'TGC': 'C',
             'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L',
             'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
             'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
             'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
             'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
             'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T',
             'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
             'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R',
             'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
             'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A',
             'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E',
             'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
             'GGG': 'G'}
    counts = {}
    for codon in aligned_codons:
        if '-' in codon:
            c_cod = 'gap'
        elif codon == stop_codon:
            c_cod = codon
        elif codon == 'TAG' or codon == 'TGA' or codon == 'TAA':
            c_cod = 'Other stop codon'
        else:
            c_cod = table[codon]
        if c_cod in counts.keys():
            counts[c_cod] += 1
        else:
            counts[c_cod] = 1
    return counts
    

def read_msa_results(output_dir, reassigned_stop_codon):
    """
    Parse all MSA results in specified directory
    :param output_dir: str, path to directory containing MSA results
    :param reassigned_stop_codon: str, reassigned stop codon
    :return: dict, counts of number of times a stop codon aligned a particular amino acid
    """
    cwd = output_dir
    to_cluster = SeqIO.to_dict(SeqIO.parse(f'{cwd}to_cluster.fna', 'fasta'))
    results = glob.glob(f'{cwd}msa_results/*[0-9].clu')
    aligned_codons = []
    for res in results:
        records = list(AlignIO.read(res, 'clustal'))
        aligned_codons.extend(get_other_codons_at_reassigned_positions(records, reassigned_stop_codon))
        for rec in records:
            rec.seq = rec.seq.lower()
            rec.seq = capitalize_in_frame_stop(rec.seq, reassigned_stop_codon)
            rec.description = to_cluster[rec.id].description
        write_to_file(records, res)
    amino_acid_counts = summarize_aligned_codons(aligned_codons, reassigned_stop_codon)
    for aa, count in amino_acid_counts.items():
        print(f"{aa}: {count}")
    return amino_acid_counts


def plot_amino_acid_counts(counts):
    """
    Visualize how frequencies with that the putatively reassigned stop codon aligns to codons coding for a
    particular amino acid.
    :param counts: dict, counts of number of times a stop codon aligned a particular amino acid
    """
    np.random.seed(13)
    fig, ax = plt.subplots()
    possible_amino_acids = ['TGA', 'TAG', 'TAA', "Other stop codon", "A", 'R', 'N', 'D', 'C', 
                            'Q', 'E', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'G', 'Other AA < 2%']
    color_mapping = {}
    color_mapping['TGA'] = 'green'
    color_mapping['TAG'] = 'red'
    color_mapping['TAA'] = 'orange'
    color_mapping['Other stop codon'] = "deepskyblue"
    color_mapping['Q'] = "magenta"
    color_mapping["L"] = "greenyellow"
    color_mapping['W'] = "darkcyan"
    color_mapping["Other AA < 2%"] = 'black'
    for paa in possible_amino_acids:
        if paa not in color_mapping:
            color_mapping[paa] =(np.random.random(), np.random.random(), np.random.random())
    present_aa = []
    handles = []
    for stop_codon, c_counts in counts.items():
        amino_acids = [x for x in c_counts.keys()]
        amino_acid_counts = [x for x in c_counts.values()]
        sorted_index = np.argsort(amino_acid_counts)[::-1]
        amino_acids = np.array(amino_acids)[sorted_index]
        amino_acid_counts = np.sort(amino_acid_counts)[::-1]
        # exclude gaps
        amino_acid_counts = amino_acid_counts[np.where(amino_acids != 'gap')[0]]
        amino_acids = amino_acids[np.where(amino_acids != 'gap')[0]]
        amino_acid_counts = amino_acid_counts / amino_acid_counts.sum() 
        cumsum = 0
        other_aa = 0
        for aa, c in zip(amino_acids, amino_acid_counts):
            if c >= 0.02:
                ax.bar(stop_codon, c, color=color_mapping[aa], bottom=cumsum)
                cumsum += c
                if aa not in present_aa:
                    present_aa.append(aa)
                    handles.append(Line2D([0], [0], c=color_mapping[aa], lw=5))
            else:
                other_aa += c
        if other_aa > 0:
            ax.bar(stop_codon, other_aa, color=color_mapping['Other AA < 2%'], bottom=cumsum)
            if 'Other AA < 2%' not in present_aa:
                present_aa.append('Other AA < 2%')
                handles.append(Line2D([0], [0], c=color_mapping['Other AA < 2%'], lw=5))
    ax.set_xlabel("Reassigned stop codon")
    ax.set_ylabel("Frequency")
    # sort legend items
    indices_present_aa = np.argsort([possible_amino_acids.index(aa) for aa in present_aa])
    sorted_handles = [handles[i] for i in indices_present_aa]
    sorted_aa = [present_aa[i] for i in indices_present_aa]
    ax.legend(sorted_handles, sorted_aa, bbox_to_anchor=(1.025, 0.75), ncol=1)
    fig.savefig('alignment_stop_codons.pdf', bbox_inches='tight', dpi=600)


def main(argv):
    parser = argparse.ArgumentParser(description="Validate predicted stop codon reassignments using MSA")
    parser.add_argument('--input_genomes', nargs="+", help='txt files containing one genome per line', required=True)
    parser.add_argument('--mgcod_results', nargs="+", help='Paths to Mgcod results', required=True)
    parser.add_argument('--mgm_results', nargs="+", help='Paths to MGM results', required=True)
    parser.add_argument('--output_dir', nargs="+", help='Paths to output directory', required=True)
    parser.add_argument('--reassigned_stop_codon', nargs="+", help='Stop codons that are reassigned', required=True)
    parser.add_argument('--clustalo_path', help='Specify clustalo path', default='clustalo')
    parser.add_argument('--mgcod_path', help='Specify Mgcod path', default='bin/mgcod.py')
    parser.add_argument('--usearch_path', help='Specify sumaclust path', default='usearch11.0.667_i86linux32')
    parser.add_argument('--log', default='./confirming_reassignment_msa.log', help='Path to Log file [./confirming_reassignment_msa.log]',
                        required=False)
    args = parser.parse_args()
    amino_acid_counts = {}
    logging.basicConfig(filename=args.log, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M:%S')
    for input_genomes, mgcod_results, mgm_results, output_dir, stop_codon in zip(args.input_genomes, args.mgcod_results,
                                                                               args.mgm_results, args.output_dir, args.reassigned_stop_codon):
        logging.info("Reading queries")
        queries = read_queries(input_genomes)
        logging.info("Predicting genes")
        gene_id_mapping = predict_genes(queries, output_dir, mgcod_results, mgm_results, args.mgcod_path)
        gene_id_contig_mapping = {k: v.split(' ')[0] for (k, v) in gene_id_mapping.items()}
        logging.info("Clustering predicted genes")
        cluster_for_msa, cluster_contig_mapping = cluster_genes(gene_id_contig_mapping, output_dir, args.usearch_path)
        logging.info("Performing MSA")
        do_msa(cluster_for_msa, output_dir, args.clustalo_path)
        logging.info("Processing MSA results")
        aa_counts = read_msa_results(output_dir, stop_codon)
        amino_acid_counts[stop_codon] = aa_counts
    logging.info("Plotting figure")
    plot_amino_acid_counts(amino_acid_counts)


if __name__ == '__main__':
    main(sys.argv[1:])
