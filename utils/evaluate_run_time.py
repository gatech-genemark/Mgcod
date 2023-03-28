#!/usr/bin/env python
import sys
import argparse
import subprocess
import time
import matplotlib.pyplot as plt
from Bio import SeqIO
import shutil
import os
import numpy as np

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--uniform', help='A file specifying one input genome with a single genetic code (in FASTA format) per line.')
    parser.add_argument('-m', '--multiple', help='A file specifying one input genome with a multiple genetic codes (in FASTA format) per line.')
    args = parser.parse_args()
    mgm_predictions = './mgm_results_tmp'
    results_path = './mgcod_results_tmp'
    if os.path.isdir(mgm_predictions):
        shutil.rmtree(mgm_predictions, ignore_errors=True)
    basic_mgcod_command = 'mgcod.py -i {genome} -p {mgm_pred} -o {results}'
    multiple_gcode_mgcod_command = 'mgcod.py -i {genome} -p {mgm_pred} -o {results} --isoforms'
    input_files = []
    length_of_seqs = []
    with open(args.uniform, 'r') as fp:
        for line in fp:
            input_files.append(line.strip())
            length_of_seqs.append(len(list(SeqIO.parse(line.strip(), 'fasta'))[0].seq))
   
    runtime_basic = []
    for i, genome in enumerate(input_files):
        start = time.time()
        subprocess.check_output(basic_mgcod_command.format(genome=genome, mgm_pred=mgm_predictions, results=results_path).split(' '))
        end = time.time()
        print(genome, length_of_seqs[i], end - start)
        runtime_basic.append(end - start)
    fig, ax = plt.subplots(1, 2, figsize=(9, 4))
    plt.subplots_adjust(wspace=0.3)
    length_of_seqs = np.array(length_of_seqs) / 1000
    sorted_ind =  np.argsort(length_of_seqs)
    length_of_seqs = np.sort(length_of_seqs)
    runtime_basic_sorted = [runtime_basic[i] for i in sorted_ind]
    ax[0].plot(length_of_seqs, runtime_basic_sorted, color='blue', marker='o')
    ax[0].text(-0.15, 1.05, 'A', transform=ax[0].transAxes, fontsize=12, fontweight='bold')
    runtime_multiple_gcodes = []
    shutil.rmtree(mgm_predictions, ignore_errors=True)
    input_files = []
    length_of_seqs = []
    with open(args.multiple, 'r') as fp:
        for line in fp:
            input_files.append(line.strip())
            length_of_seqs.append(len(list(SeqIO.parse(line.strip(), 'fasta'))[0].seq))


    # default parameters
    for i, genome in enumerate(input_files):
        start = time.time()
        subprocess.check_output(multiple_gcode_mgcod_command.format(genome=genome, mgm_pred=mgm_predictions, results=results_path).split(' '))
        end = time.time()
        print(genome, length_of_seqs[i], end - start)
        runtime_multiple_gcodes.append(end - start)
    length_of_seqs = np.array(length_of_seqs) / 1000
    sorted_ind =  np.argsort(length_of_seqs)
    length_of_seqs = np.sort(length_of_seqs)
    runtime_multiple_gcodes_sorted = [runtime_multiple_gcodes[i] for i in sorted_ind]
    ax[1].plot(length_of_seqs, runtime_multiple_gcodes_sorted, color='green', marker='o', label='-st 5000 -w 5000')

    # reduce stride
    runtime_multiple_gcodes = []
    shutil.rmtree(mgm_predictions, ignore_errors=True)
    for i, genome in enumerate(input_files):
        start = time.time()
        subprocess.check_output((multiple_gcode_mgcod_command + ' -st 2500').format(genome=genome, mgm_pred=mgm_predictions, results=results_path).split(' '))
        end = time.time()
        print(genome, length_of_seqs[i], end - start)
        runtime_multiple_gcodes.append(end - start)
    runtime_multiple_gcodes_sorted = [runtime_multiple_gcodes[i] for i in sorted_ind]
    ax[1].plot(length_of_seqs, runtime_multiple_gcodes_sorted, color='blue', marker='o', label='-st 2500 -w 5000')

    # reduce stride and window size
    runtime_multiple_gcodes = []
    shutil.rmtree(mgm_predictions, ignore_errors=True)
    for i, genome in enumerate(input_files):
        start = time.time()
        subprocess.check_output((multiple_gcode_mgcod_command + ' -st 2500 -w 2500').format(genome=genome, mgm_pred=mgm_predictions, results=results_path).split(' '))
        end = time.time()
        print(genome, length_of_seqs[i], end - start)
        runtime_multiple_gcodes.append(end - start)
    runtime_multiple_gcodes_sorted = [runtime_multiple_gcodes[i] for i in sorted_ind]
    ax[1].plot(length_of_seqs, runtime_multiple_gcodes_sorted, color='red', marker='o', label='-st 2500 -w 2500')


    ax[1].legend()
    ax[1].text(-0.15, 1.05, 'B', transform=ax[1].transAxes, fontsize=12, fontweight='bold')
    ax[0].set_xlabel('Sequence length in kb')
    ax[0].set_ylabel('Run time in s')
    ax[1].set_xlabel('Sequence length in kb')

    fig.savefig('runtime.png', dpi=600)


if __name__ == '__main__':
   main(sys.argv[1:])
