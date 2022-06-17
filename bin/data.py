#!/usr/bin/env python
# -------------------------------------------------
# Mgcod
# Aaron Pfennig
# Georgia Institute of Technology
# -------------------------------------------------

import pandas as pd
import subprocess
import os
import numpy as np
from Bio import SeqIO
import logging


class Data:
    def __init__(self, path_to_genome, path_to_predictions, verbose, cwd, species):
        """
        :param path_to_genome: str, path to genome
        :param path_to_predictions: str, path to predictions
        :param verbose: boolean, verbosity
        :param cwd: str, current working_directory
        :param species: str, species
        """
        self.path_to_genome = path_to_genome
        self.path_to_predictions = path_to_predictions
        self.verbose = verbose
        self.cwd = cwd
        self.species = species

    def run_mgm(self):
        """
        Executes all MGM models
        """
        # run all 4 mgm models
        logging.info(f'Running MGM models on {self.species}')
        if self.verbose:
            print(f'Running MGM models on {self.species}')
        if not os.path.isfile(f'{self.path_to_predictions}mgm_11_{self.species}.gff'):
            mgm_11 = [f'{self.cwd}/dependencies/gmhmmp2', '-M', f'{self.cwd}/dependencies/mgm_11.mod', '-s',
                      self.path_to_genome, '-f', 'gff', '-o', f'{self.path_to_predictions}mgm_11_{self.species}.gff']
            subprocess.check_call(mgm_11)
        if not os.path.isfile(f'{self.path_to_predictions}mgm_4_{self.species}.gff'):
            mgm_4 = [f'{self.cwd}/dependencies/gmhmmp2', '-M', f'{self.cwd}/dependencies/mgm_4.mod', '-s',
                     self.path_to_genome, '-f', 'gff', '-o', f'{self.path_to_predictions}mgm_4_{self.species}.gff']
            subprocess.check_call(mgm_4)
        if not os.path.isfile(f'{self.path_to_predictions}mgm_15_{self.species}.gff'):
            mgm_15 = [f'{self.cwd}/dependencies/gmhmmp2', '-M', f'{self.cwd}/dependencies/mgm_15.mod', '-s',
                      self.path_to_genome, '-f', 'gff', '-o', f'{self.path_to_predictions}mgm_15_{self.species}.gff']
            subprocess.check_call(mgm_15)
        if not os.path.isfile(f'{self.path_to_predictions}mgm_101_{self.species}.gff'):
            mgm_101 = [f'{self.cwd}/dependencies/gmhmmp2', '-M', f'{self.cwd}/dependencies/mgm_101.mod', '-s',
                       self.path_to_genome, '-f', 'gff', '-o', f'{self.path_to_predictions}mgm_101_{self.species}.gff']
            subprocess.check_call(mgm_101)

    def extract_gene_coords(self):
        """
        Reads MGM predictions
        :return: pandas DataFrame, for each model
        """
        # extract gene coords from gff files
        logging.info(f"Extracting gene coordinates predicted in {self.species}")
        if self.verbose:
            print(f"Extracting gene coordinates predicted in {self.species}")

        coords_11 = pd.read_csv(f"{self.path_to_predictions}mgm_11_{self.species}.gff", sep='\t',
                                header=None, comment='#',
                                names=['sequence', 'source', 'feature', 'start',
                                       'end', 'logodd', 'strand', 'phase', 'attributes'])
        coords_11['source'] = coords_11['source'] + '_11'
        coords_11['length'] = coords_11['end'] - coords_11['start'] + 1
        coords_11['sequence'] = coords_11['sequence'].astype(str)

        coords_4 = pd.read_csv(f"{self.path_to_predictions}mgm_4_{self.species}.gff", sep='\t',
                               header=None, comment='#',
                               names=['sequence', 'source', 'feature', 'start',
                                      'end', 'logodd', 'strand', 'phase', 'attributes'])
        coords_4['source'] = coords_4['source'] + '_4'
        coords_4['length'] = coords_4['end'] - coords_4['start'] + 1
        coords_4['sequence'] = coords_4['sequence'].astype(str)

        coords_15 = pd.read_csv(f"{self.path_to_predictions}mgm_15_{self.species}.gff", sep='\t',
                                header=None, comment='#',
                                names=['sequence', 'source', 'feature', 'start',
                                       'end', 'logodd', 'strand', 'phase', 'attributes'])
        coords_15['length'] = coords_15['end'] - coords_15['start'] + 1
        coords_15['source'] = coords_15['source'] + '_15'
        coords_15['sequence'] = coords_15['sequence'].astype(str)

        coords_101 = pd.read_csv(f"{self.path_to_predictions}mgm_101_{self.species}.gff", sep='\t',
                                 header=None, comment='#',
                                 names=['sequence', 'source', 'feature', 'start',
                                        'end', 'logodd', 'strand', 'phase', 'attributes'])
        coords_101['length'] = coords_101['end'] - coords_101['start'] + 1
        coords_101['source'] = coords_101['source'] + '_101'
        coords_101['sequence'] = coords_101['sequence'].astype(str)

        # simulating switch in PES
        #annotations = [(float(a), float(b)) for (a, b) in [pd.read_csv('/storage3/w/aaron/data/switching_points/annotations/switch_points.txt',
        #              header=None, index_col=0, sep='\t', names=['switch', 'switch_point']).loc[self.species, 'switch_point'][1:-1].split(', ')]]
        #coords_11['strand'] = np.where(coords_11.start.values < annotations[0][0], '+', '-')
        #coords_4['strand'] = np.where(coords_4.start.values < annotations[0][0], '+', '-')
        #coords_15['strand'] = np.where(coords_15.start.values < annotations[0][0], '+', '-')
        #coords_101['strand'] = np.where(coords_101.start.values < annotations[0][0], '+', '-')

        return coords_11, coords_4, coords_15, coords_101

    def read_dna_sequence(self):
        """
        Reads DNA Sequence from fasta file
        :return: Biopython Seq object
        """
        logging.info(f"Reading {self.path_to_genome}")
        if self.verbose:
            print(f"Reading {self.path_to_genome}")
        # read sequence
        records = list(SeqIO.parse(self.path_to_genome, 'fasta'))
        if len(records) == 1:
            seq = [records[0]]
        else:
            # sort records based on their contig id
            records = [records[x] for x in np.argsort([rec.id for rec in records])]
            seq = [records[i] for i in range(len(records))]
        return seq
