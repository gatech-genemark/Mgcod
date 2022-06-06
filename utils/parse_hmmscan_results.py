#!/usr/bin/env python
"""
Script to parse hmmscan output
"""
import pandas as pd
from Bio import SeqIO
import sys
import argparse
import numpy as np


def read_hmmscan_results(path_result, nr_queries):
    """
    Script for parsing hmmscan output files. E-values are bonferroni corrected using number queries
    :param path_result: str, path to hmmscan output file
    :param nr_queries: int, number of queries
    :return: pd.DataFrame, hmmscan results
    """
    columns = ['target', 't_acc', 'query','q_acc', 'evalue_full', 'score_full', 'bias_full',
               'evalue_domain', 'score_domain', 'bias_domain', 'exp', 'reg', 'clu', 'ov',
               'env', 'dom', 'rep', 'inc', 'description']
    results = pd.read_csv(path_result, sep='\s+', names=columns, header=None, comment='#')
    # Bonferroni correction
    results.evalue_full *= nr_queries
    results = results[results.evalue_full < 0.0001]
    results = results[results.evalue_domain < 0.01]
    results.sort_values(['evalue_domain', 'evalue_full'], inplace=True)
    results.drop_duplicates(subset='query', inplace=True)
    results.reset_index(drop=True, inplace=True)
    results.drop(['t_acc', 'q_acc', 'description'], axis=1, inplace=True)
    return results


def read_vog_table(vog_table_path, vog):
    """
    Read annotation table of virus orthologous groups
    :param vog_table_path: str, path to annotation table
    :param vog: str, virus orthologous group
    :return: pd.DataFrame, annotation of target VOG
    """
    columns = ['phylum', 'species', 'genome_accession', 'protein_accession',
               'position', 'nr_aa', 'function']
    vog_table = pd.read_csv(f'{vog_table_path}{vog}.txt', sep='\t',
                            names=columns, header=None, comment='#')
    vog_table['vog'] = vog
    return vog_table


def analyze_functions(df):
    """
    Summarize results of annotated function of hits
    :param df: pd.DataFrame, hmmscan results with db annotations
    :return: pd.DataFrame, datafrane summarizing most commonly annotated function for each query
    """
    queries = df['query'].unique()
    new_df = pd.DataFrame(columns=df.columns, dtype=object)
    print("Query\tfunction\tcount\ttotal_nr_hits")
    df['function'] = df['function'].astype(str)
    for query in queries:
        functions = df.loc[df['query'] == query, 'function'].values
        unique_func, counts = np.unique(functions, return_counts=True)
        sorted_func = unique_func[np.argsort(counts)]
        counts = np.sort(counts)
        if np.all(['hypothetical' in x or 'unnamed' in x or 'uncharacterized' in x or 'predicted protein' in x or 'unknown' in x  or 'DUF' in x or 'Uncharacterised' in x or 'predicted product' in x for x in sorted_func]):
            continue
        else:
            mode_func = sorted_func[-1]
            if not 'hypothetical' in mode_func and not 'unnamed' in mode_func and not 'uncharacterized' in mode_func and not 'predicted protein' in mode_func and not 'unknown' in mode_func and not 'DUF' in mode_func and not 'Uncharacterised' in mode_func and not 'predicted product' in mode_func:
                pass
            else:
                mode_func = sorted_func[[not 'hypothetical' in x and not 'unnamed' in x and not 'uncharacterized' in x and not 'predicted protein' in x  and not 'unknown' in x and not 'DUF' in x and not 'Uncharacterised' in x and not 'predicted product' in x for x in sorted_func]][-1]
                counts = counts[[not 'hypothetical' in x and not 'unnamed' in x for x in sorted_func]]
            print(f'{query}\t{mode_func}\t{counts[-1]}\t{len(functions)}')
            new_df = new_df.append([df[(df.function == mode_func) & (df['query'] == query)].sort_values('evalue_full').iloc[0, :]])
    new_df = new_df.reset_index(drop=False).loc[:, ["query", "function", "protein_accession", "index"]]
    new_df = new_df.rename({"query": "Cluster ID", "function" : "Predicted function", 
                            "protein_accession": "Protein accession", "index": "PVOG ID"}, axis=1)
    return new_df


def main(argv):
    parser = argparse.ArgumentParser(description="Parse hmmscan output")
    parser.add_argument('--hmmscan_result', help='hmmscan tblout file', required=True)
    parser.add_argument('--queries', help='Fasta file with queries', required=True)
    parser.add_argument('--vog_path', help='Path to VOG tables', required=True)
    parser.add_argument('--genetic_code', help='Genetic code of current protein sequences', required=True)
    args = parser.parse_args()
    queries = list(SeqIO.parse(args.queries, 'fasta'))
    results =  read_hmmscan_results(args.hmmscan_result, len(queries))
    vog_tables = []
    for vog in results.target.values:
        vog_tables.append(read_vog_table(args.vog_path, vog))
    vog_tables = pd.concat(vog_tables)
    results = results.set_index('target').join(vog_tables.set_index('vog'))
    new_df = analyze_functions(results)
    new_df.to_excel(args.hmmscan_result.replace('.tab', '_filtered.xlsx'))
    output_base = '/'.join(args.hmmscan_result.split('/')[:-1])
    with open(f"{output_base}/annotated_queries_{args.genetic_code}_pvog.csv", 'w') as annotated:
        for query in new_df['Cluster ID'].values:
            annotated.write(str(query) + ',')
    annotated.close()
    with open(f"{output_base}/annotations_{args.genetic_code}_pvog.csv", 'w') as annotation:
        for func in new_df["Predicted function"].values:
            annotation.write(func + ',')
    annotation.close()


if __name__ == '__main__':
   main(sys.argv[1:])
