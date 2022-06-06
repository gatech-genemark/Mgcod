#!/home/aaron/anaconda3/envs/geneticcode/bin/python3.7
"""
Parse fastANI results and summarize average ANI values
"""
import pandas as pd
import sys
import argparse
import numpy as np 


def parse_fastani_results(path_to_file):
    """
    Parse fastANI result to pd.DataFrame
    :param path_to_file: str, path to fastANI output file
    :return: pd.DataFrame, fastANI results
    """
    result = pd.read_csv(path_to_file, sep='\t',
                         header=None, names=['query', 'reference', 'ani', 'mapped_frags', 'total_frags'])
    return result


def main(argv):
    parser = argparse.ArgumentParser(description="Parse fastANI results and summarize average ANI values")
    parser.add_argument('--fastani_result', help='FastANI result to use', required=True)
    args = parser.parse_args()
    anis = parse_fastani_results(args.fastani_result)
    # Drop queries that were references
    anis = anis[~np.isin(anis['query'].values, anis['reference'].values)] 
    # sort by query and then by ani
    anis.sort_values(["query", "ani"], ascending=[1, 0], inplace=True)
    # keep only best match
    anis.drop_duplicates("query", inplace=True)
    # drop matches with ANI < 80%
    anis = anis[anis.ani >= 80]
    anis.to_csv(args.fastani_result + "_best_match_only", header=False, index=False, sep='\t')
    print("Mean ANI: {:.3f}".format(anis.ani.mean() / 100))
    print("STD ANI: {:.3f}".format(anis.ani.std() / 100))


if __name__ == '__main__':
    main(sys.argv[1:])
