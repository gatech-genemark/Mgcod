#!/usr/bin/env python
"""
Script to summarize Mgcod predictions on metagenomic contigs
"""
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
import argparse
import sys
from matplotlib.lines import Line2D
from tqdm import tqdm

def get_size_and_gc(path):
    """
    Get size and GC% of a particular contig
    :param path: str, path to contig
    :return: int, float, contig size and GC%
    """
    records = list(SeqIO.parse(path, 'fasta'))
    seq = records[0].seq
    if len(records) > 1:
        for i in range(1, len(records)):
            seq += records[i].seq
    size = len(seq)
    gc = GC(seq)
    return size, gc


def read_results(result_paths, genome_paths):
    """
    Parse Mgcod results
    :param result_paths: list, paths to Mgcod results (*.tab files)
    :param genome_paths: dict, mapping of contigs to contig paths
    :return: pd.DataFrame, Mgcod results
    """
    results = pd.DataFrame(columns=['genetic_code', 'segment_length', 'genome_size', 'gc', 'genes_11',
                                    'genes_4', 'genes_15', 'genes_101', 'genes_equivalent']) 
    for rp in tqdm(result_paths):
        c_result = pd.read_csv(rp, sep='\t')
        contig = c_result.iloc[0, 0]
        gp = genome_paths[contig]
        genes_11 = np.nan
        genes_15 = np.nan
        genes_4 = np.nan
        genes_101 = np.nan
        genes_equivalent = np.nan

        if not 'Segment 1' in c_result.columns:
            genetic_code = str(c_result.loc[0, 'Sequence'].astype(int))
            segment_length = np.nan
        else:
            genetic_code = c_result.iloc[1, 2:-1].values.astype(int)
            segment_length = np.array([int(x[1:-1].split(', ')[1]) - int(x[1:-1].split(', ')[0]) for x in c_result.iloc[0, 2:-1].values.tolist()])
            if c_result.iloc[2,1].endswith('11'):
                genes_11 = c_result.iloc[2,2:-1].values.astype(int)
            elif c_result.iloc[2,1].endswith('4'):
                genes_4 = c_result.iloc[2,2:-1].values.astype(int)
            elif c_result.iloc[2,1].endswith('15'):
                genes_15 = c_result.iloc[2,2:-1].values.astype(int)
            elif c_result.iloc[2,1].endswith('101'):
                genes_101 = c_result.iloc[2,2:-1].values.astype(int)
            genes_equivalent = (c_result.iloc[3,2:-1].values.astype(int) + c_result.iloc[4,2:-1].values.astype(int))
            if c_result.iloc[6,1].endswith('11'):
                genes_11 = c_result.iloc[6,2:-1].values.astype(int)
            elif c_result.iloc[6,1].endswith('4'):
                genes_4 = c_result.iloc[6,2:-1].values.astype(int)
            elif c_result.iloc[6,1].endswith('15'):
                genes_15 = c_result.iloc[6,2:-1].values.astype(int)
            elif c_result.iloc[6,1].endswith('101'):
                genes_101 = c_result.iloc[6,2:-1].values.astype(int)

        size, gc = get_size_and_gc(gp)
        results.loc[contig] = [genetic_code, segment_length, size, gc, genes_11, genes_4, genes_15, genes_101, genes_equivalent]
    return results


def filter_contigs_based_on_segment_length(results):
    """
    Filter out contigs whose shorter block is shorter than 20% of the contig
    :param results: pd.DataFrame, Mgcod results
    :return: pd.DataFrame, filtered Mgcod results
    """
    to_keep = []
    segment_lengths = results.segment_length.values
    genetic_codes = results.genetic_code.values
    genome_size = results.genome_size.values
    for sl, gc, gs in zip(segment_lengths, genetic_codes, genome_size):
        if np.all(sl / gs >= 0.2):
            to_keep.append(True)
        else:
            if len(gc) == 3 and gc[0] == gc[2] and (sl / gs)[0] + (sl /gs)[2] >= 0.2 and (sl / gs)[1] >= 0.2:
                to_keep.append(True)
            else:
                to_keep.append(False)
    results = results.loc[to_keep]
    return results


def find_contigs_with_dual_coding(results):
    """
    Extract results of metagenomic contigs with dual coded blocks
    :param results: pd.DataFrame, Mgcod results
    :return: pd.DataFrame, dict, extended dataframe indicating dual coding, mapping of genetic code and frequency with
                                 that it is dual coded in contigs with different types of multiple genetic codes
    """
    genetic_codes = results.genetic_code.values
    genes_11 = results.genes_11.values
    genes_4 = results.genes_4.values
    genes_15 = results.genes_15.values
    genes_101 = results.genes_101.values
    genes_equivalent = results.genes_equivalent.values
    dual_coded_contigs = []
    dual_coded_segments = {'11/4': [], '11/15': [], '11/101': [], '4/15': [], '4/101': [], '15/101': []}
    number_of_genes_per_genetic_code = {}
    number_of_genes_per_genetic_code[11] = {11: [], 15: [], 4: [], 101: [], 'equivalent': []}
    number_of_genes_per_genetic_code[4] = {11: [], 15: [], 4: [], 101: [], 'equivalent': []}
    number_of_genes_per_genetic_code[15] = {11: [], 15: [], 4: [], 101: [], 'equivalent': []}
    number_of_genes_per_genetic_code[101] = {11: [], 15: [], 4: [], 101: [], 'equivalent': []}
    # count how often genes are coded by which genetic code
    for contig, gc, g11, g4, g15, g101, geq in zip(results.index.values, genetic_codes, genes_11, genes_4, genes_15, genes_101, genes_equivalent):
        to_stack = []
        dual_coding = False
        if not np.all(np.isnan(g11)):
            number_of_genes_per_genetic_code[11][11].extend(g11[np.where(gc == 11)[0]].tolist())
            number_of_genes_per_genetic_code[4][11].extend(g11[np.where(gc == 4)[0]].tolist())
            number_of_genes_per_genetic_code[15][11].extend(g11[np.where(gc == 15)[0]].tolist())
            number_of_genes_per_genetic_code[101][11].extend(g11[np.where(gc == 101)[0]].tolist())
            to_stack.append(g11)

        if not np.all(np.isnan(g4)):
            number_of_genes_per_genetic_code[11][4].extend(g4[np.where(gc == 11)[0]].tolist())
            number_of_genes_per_genetic_code[4][4].extend(g4[np.where(gc == 4)[0]].tolist())
            number_of_genes_per_genetic_code[15][4].extend(g4[np.where(gc == 15)[0]].tolist())
            number_of_genes_per_genetic_code[101][4].extend(g4[np.where(gc == 101)[0]].tolist())
            to_stack.append(g4)

        if not np.all(np.isnan(g15)):
            number_of_genes_per_genetic_code[11][15].extend(g15[np.where(gc == 11)[0]].tolist())
            number_of_genes_per_genetic_code[4][15].extend(g15[np.where(gc == 4)[0]].tolist())
            number_of_genes_per_genetic_code[15][15].extend(g15[np.where(gc == 15)[0]].tolist())
            number_of_genes_per_genetic_code[101][15].extend(g15[np.where(gc == 101)[0]].tolist())
            to_stack.append(g15)

        if not np.all(np.isnan(g101)):
            number_of_genes_per_genetic_code[11][101].extend(g101[np.where(gc == 11)[0]].tolist())
            number_of_genes_per_genetic_code[4][101].extend(g101[np.where(gc == 4)[0]].tolist())
            number_of_genes_per_genetic_code[15][101].extend(g101[np.where(gc == 15)[0]].tolist())
            number_of_genes_per_genetic_code[101][101].extend(g101[np.where(gc == 101)[0]].tolist())
            to_stack.append(g101)

        if not np.all(np.isnan(geq)):
            number_of_genes_per_genetic_code[11]['equivalent'].extend(geq[np.where(gc == 11)[0]].tolist())
            number_of_genes_per_genetic_code[4]['equivalent'].extend(geq[np.where(gc == 4)[0]].tolist())
            number_of_genes_per_genetic_code[15]['equivalent'].extend(geq[np.where(gc == 15)[0]].tolist())
            number_of_genes_per_genetic_code[101]['equivalent'].extend(geq[np.where(gc == 101)[0]].tolist())
            to_stack.append(geq)
            gene_counts = np.stack(to_stack)
            if np.any(np.argmax(gene_counts, axis=0) == gene_counts.shape[0] - 1):
                dual_coding = True
                if 11 in gc and 4 in gc:
                    dual_coded_segments['11/4'].extend(gc[np.argmax(gene_counts, axis=0) == gene_counts.shape[0] - 1].tolist())
                elif 11 in gc and 15 in gc:
                    dual_coded_segments['11/15'].extend(gc[np.argmax(gene_counts, axis=0) == gene_counts.shape[0] - 1].tolist())
                elif 11 in gc and 101 in gc:
                    dual_coded_segments['11/101'].extend(gc[np.argmax(gene_counts, axis=0) == gene_counts.shape[0] - 1].tolist())
                elif 4 in gc and 15 in gc:
                    dual_coded_segments['4/15'].extend(gc[np.argmax(gene_counts, axis=0) == gene_counts.shape[0] - 1].tolist())
                elif 4 in gc and 101 in gc:
                    dual_coded_segments['4/101'].extend(gc[np.argmax(gene_counts, axis=0) == gene_counts.shape[0] - 1].tolist())
                elif 15 in gc and 101 in gc:
                    dual_coded_segments['15/101'].extend(gc[np.argmax(gene_counts, axis=0) == gene_counts.shape[0] - 1].tolist())
        dual_coded_contigs.append(dual_coding)
    results['dual_coded'] = dual_coded_contigs
    return results, dual_coded_segments


def summarize_results(results, output_path):
    """
    Summarize Mgcod results on metagenomic contigs
    :param results: pd.DataFrame, Mgcod predictions
    :param output_path: str, path to where to save summary file
    :return: pd.DataFrame, dict, extended dataframe indicating dual coding, mapping of genetic code and frequency with
                                 that it is dual coded in contigs with different types of multiple genetic codes
    """
    print(f"Total {results.shape[0]} samples")
    print(f"11\t{results[results.genetic_code.values == '11'].shape[0]}")
    print(f"15\t{results[results.genetic_code == '15'].shape[0]}")
    print(f"101\t{results[results.genetic_code == '101'].shape[0]}")
    print(f"4\t{results[results.genetic_code == '4'].shape[0]}")

    # drop samples without switch point
    inds_to_keep = [ind for ind, genetic_code in zip(results.index, results.genetic_code) if isinstance(genetic_code, np.ndarray)]
    results = results.loc[inds_to_keep]
    print(f"11/4\t{len([x for x in results.genetic_code.values if 11 in x and 4 in x])}")
    print(f"11/15\t{len([x for x in results.genetic_code.values if 11 in x and 15 in x])}")
    print(f"11/101\t{len([x for x in results.genetic_code.values if 11 in x and 101 in x])}")
    print(f"101/15\t{len([x for x in results.genetic_code.values if 101 in x and 15 in x])}")
    print(f"101/4\t{len([x for x in results.genetic_code.values if 101 in x and 4 in x])}")
    print(f"15/4\t{len([x for x in results.genetic_code.values if 15 in x and 4 in x])}")

    # remove predictions with more than 2 switch points
    nr_switches = [x.shape[0] for x in results.genetic_code.values]
    keep = np.where(np.array(nr_switches) <= 3)[0]
    results = results.iloc[keep]
    #results.reset_index(drop=True, inplace=True)
    print("\nRemoved samples with more than 2 code switches:")
    print(f"11/4\t{len([x for x in results.genetic_code.values if 11 in x and 4 in x])}")
    print(f"11/15\t{len([x for x in results.genetic_code.values if 11 in x and 15 in x])}")
    print(f"11/101\t{len([x for x in results.genetic_code.values if 11 in x and 101 in x])}")
    print(f"101/15\t{len([x for x in results.genetic_code.values if 101 in x and 15 in x])}")
    print(f"101/4\t{len([x for x in results.genetic_code.values if 101 in x and 4 in x])}")
    print(f"15/4\t{len([x for x in results.genetic_code.values if 15 in x and 4 in x])}")

    results = filter_contigs_based_on_segment_length(results)
    print("\nFiltered out samples with shortest block being less than 1/5 of genome:")
    print(f"11/4\t{len([x for x in results.genetic_code.values if 11 in x and 4 in x])}")
    print(f"11/15\t{len([x for x in results.genetic_code.values if 11 in x and 15 in x])}")
    print(f"11/101\t{len([x for x in results.genetic_code.values if 11 in x and 101 in x])}")
    print(f"101/15\t{len([x for x in results.genetic_code.values if 101 in x and 15 in x])}")
    print(f"101/4\t{len([x for x in results.genetic_code.values if 101 in x and 4 in x])}")
    print(f"15/4\t{len([x for x in results.genetic_code.values if 15 in x and 4 in x])}")

    results, dual_coded_segments = find_contigs_with_dual_coding(results)
    print("\nSamples with dual coding:")
    print(f"11/4\t{len([x for x, y in zip(results.genetic_code.values, results.dual_coded.values) if 11 in x and 4 in x and y])}")
    print(f"11/15\t{len([x for x, y in zip(results.genetic_code.values, results.dual_coded.values) if 11 in x and 15 in x and y])}")
    print(f"11/101\t{len([x for x, y in zip(results.genetic_code.values, results.dual_coded.values) if 11 in x and 101 in x and y])}")
    print(f"101/15\t{len([x for x, y in zip(results.genetic_code.values, results.dual_coded.values) if 101 in x and 15 in x and y])}")
    print(f"101/4\t{len([x for x, y in zip(results.genetic_code.values, results.dual_coded.values) if 101 in x and 4 in x and y])}")
    print(f"15/4\t{len([x for x, y in zip(results.genetic_code.values, results.dual_coded.values) if 15 in x and 4 in x and y])}")
    results.to_csv(output_path, sep='\t', index=True, header=True)
    return results, dual_coded_segments


def plot_frequency_with_that_a_genetic_code_is_distinctly_encoded_in_one_block(dual_coded_segments, plot_dir):
    """
    Plots the frequency with that a black is distinctly encoded in a particular genetic code in phages with
    different types of multiple genetic codes
    :param dual_coded_segments: dict, mapping of genetic code and frequency with
                                      that it is dual coded in contigs with different types of multiple genetic codes
    :param plot_dir: str, directory where to save figure (frequency_of_blocks_with_distinct_code)
    """
    colors = {11: 'blue', 4: 'green', 15: 'red', 101: 'orange'}
    fig, ax = plt.subplots()
    x_labels = []
    for i, (key, values) in enumerate(dual_coded_segments.items()):
        if len(values) == 0:
            continue
        x_labels.append(key)
        if key == '11/4':
            ax.bar(i - 0.2, 1 - values.count(11) / len(values), width=0.35, color=colors[11])
            ax.bar(i + 0.2, 1 - values.count(4) / len(values), width=0.35, color=colors[4])
        elif key == '11/15':
            ax.bar(i - 0.2, 1 - values.count(11) / len(values), width=0.35, color=colors[11])
            ax.bar(i + 0.2, 1 - values.count(15) / len(values), width=0.35, color=colors[15])
        elif key == '11/101':
            ax.bar(i - 0.2, 1 - values.count(11) / len(values), width=0.35, color=colors[11])
            ax.bar(i + 0.2, 1 - values.count(101) / len(values), width=0.35, color=colors[101])
        elif key == '4/15':
            ax.bar(i - 0.2, 1 - values.count(4) / len(values), width=0.35, color=colors[4])
            ax.bar(i + 0.2, 1 - values.count(15) / len(values), width=0.35, color=colors[15])
        elif key == '4/101':
            ax.bar(i - 0.2, 1 - values.count(4) / len(values), width=0.35, color=colors[4])
            ax.bar(i + 0.2, 1 - values.count(101) / len(values), width=0.35, color=colors[101])
        elif key == '15/101':
            ax.bar(i - 0.2, 1 - values.count(15) / len(values), width=0.35, color=colors[15])
            ax.bar(i + 0.2, 1 - values.count(101) / len(values), width=0.35, color=colors[101])
    ax.set_xticks([i for i in range(len(x_labels))])
    ax.set_xticklabels(x_labels)
    ax.set_ylabel("Frequency of blocks with distinct genetic code")
    ax.set_xlabel('Type of multiple genetic codes')
    legend_elements = []
    genetic_codes = sorted(set([y for x in x_labels for y in x.split('/')]))
    for gc in genetic_codes:
        if gc == '11':
            legend_elements.append(Line2D([0], [0], color='blue', label='Code 11'))
        if gc == '4':
            legend_elements.append(Line2D([0], [0], color='green', label='Code 4'))
        if gc == '15':
            legend_elements.append(Line2D([0], [0], color='red', label='Code 15'))
        if gc == '101':
            legend_elements.append(Line2D([0], [0], color='orange', label='Code 101'))
    ax.legend(legend_elements, ["Code " + gc for gc in genetic_codes], loc='upper center', bbox_to_anchor=(0.5, -.13), ncol=4)
    fig.savefig(f'{plot_dir}frequency_of_blocks_with_distinct_code.pdf', bbox_inches='tight', dpi=600)


def main(argv):
    parser = argparse.ArgumentParser(description="Summarizes Mgcod predictions on metagenomic contigs")
    parser.add_argument("-r", "--results", help='File with paths to result files, one path per line', required=True)
    parser.add_argument("-g", "--genomes", help="File with paths to genome files, one path per line", required=True)
    parser.add_argument("-o", "--output", help="Path to file to write summarized dataframe to", required=True)
    parser.add_argument('-p', '--plot_dir', help='Directory where to save plots', required=True)
    args = parser.parse_args()

    result_paths = []
    with open(args.results, 'r') as results:
        for line in results:
            result_paths.append(line.strip())
    results.close()
    
    genome_paths = {}
    with open(args.genomes, 'r') as genomes:
        for line in genomes:
            contig = line.strip().split('_')[-1].split('.fna')[0]
            genome_paths[contig] = line.strip()
    genomes.close()
    results = read_results(result_paths, genome_paths)   
    results, dual_coded_segments = summarize_results(results, args.output)
    plot_frequency_with_that_a_genetic_code_is_distinctly_encoded_in_one_block(dual_coded_segments, args.plot_dir)


if __name__ == '__main__':
    main(sys.argv[1:])
