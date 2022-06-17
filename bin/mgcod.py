#!/usr/bin/env python
# -------------------------------------------------
# Mgcod
# Aaron Pfennig
# Georgia Institute of Technology
# -------------------------------------------------

import os
import sys
import argparse
import re
from datetime import datetime
from pipeline import Pipeline
from visualizations import Visualizations
from results import Results
import warnings
import logging
import time
warnings.simplefilter(action='ignore', category=UserWarning)


def main(argv):
    warnings.simplefilter('always', UserWarning)
    date = datetime.now()
    cwd = '/'.join(sys.path[0].split('/')[:-1])
    parser = argparse.ArgumentParser(description="Mgcod segments contigs based on genetic "
                                                 "code usage and performs genetic-code-informed annotation of coding "
                                                 "regions using MetaGeneMark. Can be run on any prokaryotic sequences "
                                                 "with(-out) stop codon reassignment", prog="mgcod.py")
    optional_args = parser._action_groups.pop()

    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('-i', '--path_to_genome', help="Path to input file, FASTA format",
                               required=True, dest='path_to_genome')
    optional_args.add_argument('-p', '--path_to_mgm_predictions', help="Directory where to save MGM predictions so that"
                                                                       " they can be re-used. If path does not exist, "
                                                                       "it will be created. [./mgm_results/]",
                               required=False, dest='path_to_predictions', default=f'{cwd}/mgm_results/')
    optional_args.add_argument('-o', '--path_to_output', help='Directory where to save final gene annotations and '
                                                              'supporting outputs. If path does not exist, it will be '
                                                              'created. If -AA or -NT flag is set, sequences will be '
                                                              'saved here, too. [./results/]',
                               required=False, dest='output', default=f'{cwd}/results/')
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
    optional_args.add_argument('-st', '--stride', help='Step size in bp, with which window will be moved along '
                                                       'sequence. Only relevant for isoform prediction. [5000]',
                               default=5000, required=False, type=int, dest='stride')
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
    optional_args.add_argument('-v', '--verbose', help='verbose', required=False, dest='verbose', action='store_true',
                               default=False)
    optional_args.add_argument('--version', action='version', version='%(prog)s 1.0.0')
    parser._action_groups.append(optional_args)
    args = parser.parse_args()
    args = vars(args)
    path_to_genome = args['path_to_genome']
    path_to_predictions = args['path_to_predictions']
    plot_dir = args['path_to_plots']
    verbose = args['verbose']
    path_to_output = args['output']
    circular = args['circular']
    isoforms = args['isoforms']
    stride = args['stride']
    tolerance = args['tolerance']
    window_size = args['window_size']
    consecutive_windows = args['consecutive_windows']
    consecutive_gene_labels = args['consecutive_gene_labels']
    nucleotides = args['nucleotides']
    amino_acids = args['amino_acids']
    short_contigs = args['short_contigs']
    delete = args['delete']
    if not args['logfile']:
        current_time = time.localtime()
        logfile = 'log/mgcod_{}{}{}{}{}{}.log'.format(current_time.tm_year, current_time.tm_mon, current_time.tm_mday,
                                           current_time.tm_hour, current_time.tm_min, current_time.tm_sec)
    else:
        logfile = args['logfile']
    path_to_logfile = '/'.join(logfile.split('/')[:-1])
    if not os.path.isdir(path_to_logfile):
        os.makedirs(path_to_logfile)
    logging.basicConfig(filename=logfile, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M:%S')
    logger = logging.getLogger('Mgcod')
    # check if genome exists
    if not os.path.isfile(path_to_genome):
        sys.exit("{} does not exist".format(path_to_genome))
    # create the required directories if they do not exist
    if not os.path.isdir(path_to_predictions):
        os.makedirs(path_to_predictions)
    if not os.path.isdir(path_to_output):
        os.makedirs(path_to_output)
    if plot_dir is not None and not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)
    if not path_to_predictions.endswith('/'):
        path_to_predictions += '/'
    if not path_to_output.endswith('/'):
        path_to_output += '/'
    if not plot_dir is None and not plot_dir.endswith('/'):
        plot_dir += '/'
    suffix = '.' + path_to_genome.split('.')[-1]
    species = path_to_genome.split('/')[-1].split(suffix)[0]
    species = re.sub("[']", "\\'", species)
    species = re.sub("[(]", "\\(", species)
    species = re.sub("[)]", "\\)", species)

    # run pipeline
    pipeline = Pipeline(path_to_genome, path_to_predictions, species, isoforms, circular,
                        consecutive_windows, consecutive_gene_labels, window_size, stride, tolerance, 
                        short_contigs, verbose, cwd)
    pipeline()
    if plot_dir is not None and isoforms:
        logging.info(f'Plotting log-odds ratios for {species}')
        if verbose:
            print(f"Plotting log-odds ratios for {species}")
        for n, seq in enumerate(pipeline.sequence):
            length_sequence = len(seq)
            contig = seq.id
            if pipeline.predicted_switch_regions[contig] is None:
                continue
            else:
                contig_labels, switches = pipeline.predicted_switch_regions[contig]
                compared_gene_set = pipeline.coords_w_isoforms[pipeline.coords_w_isoforms.contig == contig]
                compared_gene_set = compared_gene_set[compared_gene_set.logodd != '.'].reset_index(drop=True)
                gene_label_to_gcode = pipeline.gene_label_to_gcode[contig]
                switch_regions_for_plots = [x for x in switches]
                labels = [x for x in contig_labels]
            pars = {'logodd_scores': pipeline.logodd_scores[n],
                    'length_sequence': length_sequence,
                    'sequence': seq,
                    'window_size': pipeline.window_size,
                    'stride': pipeline.stride,
                    'species': species,
                    'switch_regions': switch_regions_for_plots,
                    'labels': labels,
                    'plot_dir': plot_dir,
                    'gene_label_to_gcode': gene_label_to_gcode,
                    'compared_gene_set': compared_gene_set}
            visualizations = Visualizations(pars, verbose)
            visualizations()
    # write results
    results = Results(pipeline, species, path_to_predictions, path_to_genome, path_to_output, date, delete,
                      amino_acids, nucleotides, verbose, cwd)
    results()
    if verbose:
        sys.exit('----------------------------Done----------------------------')
    else:
        sys.exit()


if __name__ == '__main__':
    main(sys.argv[1:])
