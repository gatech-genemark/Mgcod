#!/usr/bin/env python
# -------------------------------------------------
# Mgcod
# Aaron Pfennig
# Georgia Institute of Technology
# -------------------------------------------------

import multiprocessing as mp
import sys 
import argparse
from itertools import islice, takewhile, repeat
import subprocess
import logging
import os
from multiprocessing_logging import install_mp_handler
import time


def get_command(args):
        cwd = sys.path[0]
        command = ['{}/mgcod.py'.format(cwd), '-o', args.path_to_output, '-p', args.path_to_predictions, '-n',
                   str(args.consecutive_windows), '-g', str(args.consecutive_gene_labels), '-w', str(args.window_size),
                   '-t', str(args.tolerance)]
        if args.path_to_plots:
            command.append('-m')
            command.append(args.path_to_plots)
        if args.circular:
            command.append('-r')
        if args.isoforms:
            command.append('--isoforms')
        if args.delete:
            command.append('-d')
        if args.amino_acids:
            command.append('-AA')
        if args.nucleotides:
            command.append('-NT')
        if args.short_contigs:
            command.append('--short_contigs')
        if args.stride:
            command.append('-st')
            command.append(str(args.stride)) 
        logging.info(f"Base command: {' '.join(command)}")
        return command


def call_mgcod(args):
    genomes, base_command = args
    for g in genomes:
        command = base_command.copy()
        command.append('-i')
        command.append(g)
        output = subprocess.check_call(command, stderr=subprocess.STDOUT)
        if output != 0:
            print('{}: {}'.format(g, output))   


def main(argv):
    cwd = sys.path[0]
    parser = argparse.ArgumentParser(description='Runs Mgcod on many contigs in parallel.')
    optional_args = parser._action_groups.pop()

    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('--genomes', help='path to txt file with genomes', required=True)
    optional_args.add_argument('-p', '--path_to_mgm_predictions', help="Directory where to save MGM predictions so that"
                                                                       " they can be re-used. If path does not exist, "
                                                                       "it will be created. [./mgm_results/]",
                               required=False, dest='path_to_predictions', default=f'{cwd}/mgm_results/')
    optional_args.add_argument('-o', '--path_to_output', help='Directory where to save final gene annotations and '
                                                              'supporting outputs. If path does not exist, it will be '
                                                              'created. If -AA or -NT flag is set, sequences will be '
                                                              'saved here, too. [./results/]',
                               required=False, dest='path_to_output', default=f'{cwd}/results/')
    optional_args.add_argument('-m', '--path_to_plots', help='Directory where to save plots. Plots logodd ratio per '
                                                             'window for different MGM models. Only available with '
                                                             'isoform prediction. If path does not exist, it will be '
                                                             'created',
                               required=False, dest='path_to_plots', default=None)
    optional_args.add_argument('-r', '--circular', help='Set if sequence is circular. Only relevant for isoform '
                                                        'prediction',
                               action='store_true', required=False, default=False, dest='circular')
    optional_args.add_argument('--isoforms', help='Enable prediction of isoforms', required=False, default=False,
                               action='store_true', dest='isoforms')
    optional_args.add_argument('-n', '--consecutive_windows', help='Number of consecutive windows to be required with '
                                                                   'same genetic code to keep. Only relevant for '
                                                                   'isoform prediction. Minimum is 2. [3]',
                               default=3, required=False, dest='consecutive_windows', type=int)
    optional_args.add_argument('-g', '--consecutive_gene_labels', help='Number of consecutive gene labels to be '
                                                                       'required with same genetic code to keep. Only '
                                                                       'relevant for isoform prediction. Minimum is 2. '
                                                                       '[5]',
                               default=5, required=False, dest='consecutive_gene_labels', type=int)
    optional_args.add_argument('-w', '--window_size', help="Window size in bp applied to search for isoform. Only "
                                                           "relevant for isoform prediction. [5000]",
                               default=5000, dest='window_size', required=False, type=int)
    optional_args.add_argument('-st', '--stride', help='Step size in bp, with which window will be moved along sequence'
                                                       '. Only relevant for isoform prediction. [If sequence <= 100 kb 2500 bp else 5000 bp]',
                               default=None, required=False, type=int, dest='stride')
    optional_args.add_argument('-t', '--tolerance', help='The maximally tolerated difference in prediction of gene '
                                                         'start or gene stop to consider the prediction of two models '
                                                         'isoforms. Only relevent for isoform prediction. [30]',
                               default=30, type=int, required=False, dest='tolerance')
    optional_args.add_argument('-d', '--delete', help="Delete intermediary files (prediction of the different MGM "
                                                      "models).",
                               action='store_true', dest='delete', required=False, default=False)
    optional_args.add_argument('-AA', '--amino_acids', help='Extract amino acid sequences of predicted proteins.',
                               default=False, action='store_true', required=False, dest='amino_acids')
    optional_args.add_argument('-NT', '--nucleotides', help='Extract nucleotide sequences of predicted proteins.',
                               default=False, action='store_true', required=False, dest='nucleotides')
    optional_args.add_argument('--short_contigs', help='Predict genetic codes of contigs < 5000bp. Prediction may not '
                                                       'be reliable', default=False, action='store_true',
                               required=False)
    optional_args.add_argument('--logfile', required=False, help='Path to log file')
    optional_args.add_argument('--threads', help='Number of processes. [8]', type=int, default=8)
    optional_args.add_argument('--batch_size', help='Batch size for multiprocessing. [256]', type=int, default=256,
                               required=False)
    optional_args.add_argument('--version', action='version', version='%(prog)s 1.0.2')
    parser._action_groups.append(optional_args)
    args = parser.parse_args()
    if not os.path.isdir(args.path_to_output):
        os.makedirs(args.path_to_output)
    if not args.logfile:
        current_time = time.localtime()
        logfile = 'log/mgcod_{}{}{}{}{}{}.log'.format(current_time.tm_year, current_time.tm_mon, current_time.tm_mday,
                                                    current_time.tm_hour, current_time.tm_min, current_time.tm_sec)
    else:
        logfile = args.logfile
    path_to_logfile = '/'.join(logfile.split('/')[:-1])
    if not os.path.isdir(path_to_logfile):
        os.makedirs(path_to_logfile)
    logging.basicConfig(filename=logfile, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M:%S')

    logger = logging.getLogger('Multiprocessed Mgcod')
    install_mp_handler()
    genomes = []
    with open(args.genomes) as in_file:
        for line in in_file:
            genomes.append(line.strip())
    base_command = get_command(args)
    iterator = iter(genomes)
    chunks = [(c, base_command) for c in takewhile(bool, (list(islice(iterator, args.batch_size)) for _ in repeat(None)))]
    pool = mp.Pool(processes=args.threads)
    pool.map(call_mgcod, chunks)
    pool.close()
    pool.join()


if __name__ == '__main__':
    main(sys.argv[1:])
