#!/usr/bin/env python
# -------------------------------------------------
# AlternatingGeneticCodeAnnotator
# Aaron Pfennig
# Georgia Institute of Technology
# -------------------------------------------------

import warnings
import pandas as pd
import numpy as np
from Bio.SeqUtils import GC
import logging
warnings.simplefilter(action='ignore', category=UserWarning)


class Pipeline:
    def __init__(self, path_to_genome, path_to_predictions, species, isoforms, circular,
                 consecutive_windows, consecutive_gene_labels, window_size, stride, tolerance,
                 short_contigs, verbosity, working_dir):
        """
        Class to run the analysis
        :param path_to_genome: str, path to genome
        :param path_to_predictions: str, path to predictions
        :param species: str, species
        :param isoforms: boolean, predict isoforms
        :param circular: boolean, ciruclar genome
        :param consecutive_windows: int, number of consecutive windows with same genetic code
        :param consecutive_gene_labels: int, number of consecutive genes with same genetic code
        :param window_size: int, window size in bp
        :param stride: int, stride in bp
        :param tolerance:int, max tolerated difference in start/stop coordinates to consider gene prediction identical
        :param short_contigs: boolean, attempt to predict contigs shorter than 5000bp
        :param verbosity: boolean, verbosity
        :param working_dir: str, current working_directory
        :rtype: object
        """
        # intialize attributes
        self.verbose = verbosity
        self.cwd = working_dir
        self.path_to_genome = path_to_genome
        self.path_to_predictions = path_to_predictions
        self.species = species
        self.circular = circular
        self.tolerance = tolerance
        self.predicted_switch_regions = None
        self.window_size = window_size
        self.stride = stride
        self.isoforms = isoforms
        self.short_contigs = short_contigs
        self.consecutive_windows = consecutive_windows
        self.consecutive_gene_labels = consecutive_gene_labels
        self.genome_label = 11
        self.label_to_gcode = {0: 11, 1: 4, 2: 15, 3: 101}
        self.gcode_to_label = {11: 0, 4: 1, 15: 2, 101: 3}

    def __call__(self):
        from data import Data
        data = Data(self.path_to_genome, self.path_to_predictions, self.verbose, self.cwd, self.species)
        # load sequence
        self.sequence = data.read_dna_sequence()
        # add contig information
        self.nr_contigs = len(self.sequence)
        # run MGM
        data.run_mgm()
        # load MGM predictions
        self.gene_coords_11, self.gene_coords_4, \
        self.gene_coords_15, self.gene_coords_101 = data.extract_gene_coords()
        contig_ids_11 = [x for x in self.gene_coords_11.sequence.values]
        contig_ids_4 = [x for x in self.gene_coords_4.sequence.values]
        contig_ids_15 = [x for x in self.gene_coords_15.sequence.values]
        contig_ids_101 = [x for x in self.gene_coords_101.sequence.values]

        self.gene_coords_11.loc[:, 'contig'] = contig_ids_11
        self.gene_coords_4.loc[:, 'contig'] = contig_ids_4
        self.gene_coords_15.loc[:, 'contig'] = contig_ids_15
        self.gene_coords_101.loc[:, 'contig'] = contig_ids_101

        self.label_to_gene_coords = {0: self.gene_coords_11, 1: self.gene_coords_4,
                                     2: self.gene_coords_15, 3: self.gene_coords_101}
        # classify complete sequence
        if not self.isoforms:
            self.genome_label = self.classify_complete_sequence()
            # naming for compatibility only --> no isoforms are predicted in this mode
            self.coords_w_isoforms = self.label_to_gene_coords[self.gcode_to_label[self.genome_label]]
        # predict isoforms and segmentate genome
        else:
            # identify two most likely genetic codes
            self.genome_label, self.logodd_scores = self.identify_most_common_codes()
            # update window labels
            # add padding for circular genome
            if self.circular and self.nr_contigs == 1:
                beginning = self.genome_label[0][:self.consecutive_windows].copy()
                end = self.genome_label[0][-self.consecutive_windows:].copy()
                if beginning.shape[0] < self.consecutive_windows or end.shape[0] < self.consecutive_windows:
                    self.consecutive_windows = min(beginning.shape[0], end.shape[0])
                    warnings.warn(f"Setting consecutive windows to {self.consecutive_windows}")
                self.genome_label = [np.concatenate([end, self.genome_label[0].copy(), beginning])]
            genome_labels = np.concatenate(self.genome_label)
            updated_labels = []
            for window_labels in self.genome_label:
                if np.unique(window_labels).shape[0] == 1:
                    updated_labels.append(window_labels)
                else:
                    updated_labels.append(self.find_best_sequence_of_labels(window_labels, genome_labels,
                                                                            self.consecutive_windows))
                    # remove padding
            if self.circular and self.nr_contigs == 1:
                updated_labels = [updated_labels[0][self.consecutive_windows:-self.consecutive_windows]]
                self.genome_label = [self.genome_label[0][self.consecutive_windows:-self.consecutive_windows]]
            self.updated_genome_labels = updated_labels
            # unique genetic code predicted
            if np.unique(np.concatenate(self.updated_genome_labels)).shape[0] == 1:
                logging.info(f"{self.species} is genetic code {self.label_to_gcode[self.updated_genome_labels[0][0]]}")
                if self.verbose:
                    print(f"{self.species} is genetic code {self.label_to_gcode[self.updated_genome_labels[0][0]]}")
                self.genome_label = self.label_to_gcode[self.updated_genome_labels[0][0]]
                # naming only for compatibility --> no isoforms predicted in this case
                self.coords_w_isoforms = self.label_to_gene_coords[self.gcode_to_label[self.genome_label]]
            # check for isoforms and segementate genome
            else:
                # compare prediction of the two models --> determine isoforms
                self.coords, self.gene_label_to_gcode, self.coords_w_isoforms = self.compare_sets_of_predicted_genes()
                # self.confident_1, self.confident_2, self.ambiguous, self.segment_labels = \
                # self.segmentate_genome_based_on_gene_labels()
                self.confident_1, self.confident_2, self.segment_labels = self.segment_genome_based_on_window_labels()
                self.predicted_switch_regions = {}
                for n, seq in enumerate(self.sequence):
                    predicted_switch_regions = sorted([seg for segments in [self.confident_1[seq.id],
                                                                            self.confident_2[seq.id]] for seg in
                                                       segments])
                    self.predicted_switch_regions[seq.id] = (
                    np.unique(self.updated_genome_labels[n]), [(predicted_switch_regions[i][1],
                                                                predicted_switch_regions[i + 1][0])
                                                               for i in range(len(predicted_switch_regions) - 1)])

    def classify_complete_sequence(self):
        """
        Classify complete sequence based on highest scoring MGM model
        :return: int, genetic code of sequence
        """
        logging.info(f"Classifying {self.species}")
        if self.verbose:
            print(f"Classifying {self.species}")
        # get model with highest total logodd score
        max_model = np.argmax([self.gene_coords_11['logodd'].sum(),
                               self.gene_coords_4['logodd'].sum(),
                               self.gene_coords_15['logodd'].sum(),
                               self.gene_coords_101['logodd'].sum()])
        logging.info(f"{self.species} is genetic code {self.label_to_gcode[max_model]}")
        if self.verbose:
            print(f"{self.species} is genetic code {self.label_to_gcode[max_model]}")
        return self.label_to_gcode[max_model]

    def identify_most_common_codes(self):
        """
        Applies sliding window approach to find most common genetic codes
        :return: numpy array sequence of labels, numpy array with logodd scores for each MGM model per window
        """
        logging.info(f"Identifying most common genetic codes in {self.species}")
        if self.verbose:
            print(f"Identifying most common genetic codes in {self.species}")
        length_seq = len(self.sequence[0].seq)
        contigs_seen = 0
        contig = self.sequence[contigs_seen].id
        if length_seq < 5000:
            logging.warning(f"Contig {contig} in {self.species} is shorter than 5000bp. "
                            f"Classification results may be not as reliable")
        start = 0
        end = self.window_size
        labels = []
        logodd_scores = []
        contig_labels = []
        contig_logodd_scores = []
        starts = []
        short_contigs = []
        short_contigs_not_classified = []
        # apply sliding window approach
        while start <= length_seq:
            starts.append(start)
            if length_seq <= 5000 and not self.short_contigs:
                short_contigs_not_classified.append(contig)
                contig_labels.append(0)
                contig_logodd_scores.append((0, 0, 0, 0))
                labels.append(np.array([0]))
                logodd_scores.append(np.array([(0, 0, 0, 0)]))
                contigs_seen += 1

                if contigs_seen > self.nr_contigs - 1:
                    break
                contig = self.sequence[contigs_seen].id
                start = 0
                end = self.window_size
                contig_labels = []
                contig_logodd_scores = []
                length_seq = len(self.sequence[contigs_seen].seq)
                continue

            # if contig is short so that less than n windows would fit in there do not apply sliding window approach
            if length_seq <= self.window_size + (self.consecutive_windows - 1) * self.stride:
                # get contig label
                contig_labels.append(
                    np.argmax([self.gene_coords_11[self.gene_coords_11.contig == contig]['logodd'].sum(),
                               self.gene_coords_4[self.gene_coords_4.contig == contig]['logodd'].sum(),
                               self.gene_coords_15[self.gene_coords_15.contig == contig]['logodd'].sum(),
                               self.gene_coords_101[self.gene_coords_101.contig == contig]['logodd'].sum()]))
                # get highest logodd score
                contig_logodd_scores.append((self.gene_coords_11[self.gene_coords_11.contig == contig]['logodd'].sum(),
                                             self.gene_coords_4[self.gene_coords_4.contig == contig]['logodd'].sum(),
                                             self.gene_coords_15[self.gene_coords_15.contig == contig]['logodd'].sum(),
                                             self.gene_coords_101[self.gene_coords_101.contig == contig][
                                                 'logodd'].sum()))
                if len(contig_labels) == 0:
                    contig_labels.append(0)
                    contig_logodd_scores.append((0, 0, 0, 0))
                if length_seq > self.window_size and len(contig_labels) == 1:
                    contig_labels = contig_labels * int(np.ceil((length_seq - self.window_size) / self.stride) + 1)
                    contig_logodd_scores = contig_logodd_scores * int(
                        np.ceil((length_seq - self.window_size) / self.stride) + 1)
                # prepare next iteration
                labels.append(np.array(contig_labels))
                logodd_scores.append(np.array(contig_logodd_scores))
                if length_seq < 5000:
                    short_contigs.append(contig)

                start = 0
                end = self.window_size
                contigs_seen += 1
                if contigs_seen > self.nr_contigs - 1:
                    break
                length_seq = len(self.sequence[contigs_seen].seq)
                contig = self.sequence[contigs_seen].id
                contig_labels = []
                contig_logodd_scores = []
                continue
            # get logodd sum of each model for current window
            window_sum_11 = self.gene_coords_11[(self.gene_coords_11['start'] >= start) &
                                                (self.gene_coords_11['start'] < end) &
                                                (self.gene_coords_11['contig'] == contig)]['logodd'].sum()
            window_sum_4 = self.gene_coords_4[(self.gene_coords_4['start'] >= start) &
                                              (self.gene_coords_4['start'] < end) &
                                              (self.gene_coords_4['contig'] == contig)]['logodd'].sum()
            window_sum_15 = self.gene_coords_15[(self.gene_coords_15['start'] >= start) &
                                                (self.gene_coords_15['start'] < end) &
                                                (self.gene_coords_15['contig'] == contig)]['logodd'].sum()
            window_sum_101 = self.gene_coords_101[(self.gene_coords_101['start'] >= start) &
                                                  (self.gene_coords_101['start'] < end) &
                                                  (self.gene_coords_101['contig'] == contig)]['logodd'].sum()
            # if log odd score is zero but there is gene spanning over the entire window take the logodd of that gene
            try:
                if window_sum_11 == 0.0 and self.gene_coords_11[(self.gene_coords_11['start'] < start) &
                                                                (self.gene_coords_11['contig'] == contig)].iloc[
                    -1, 4] > end:
                    window_sum_11 = self.gene_coords_11[(self.gene_coords_11['start'] < start) &
                                                        (self.gene_coords_11['contig'] == contig)].iloc[-1, 5]
            except IndexError:
                pass
            try:
                if window_sum_4 == 0.0 and self.gene_coords_4[(self.gene_coords_4['start'] < start) &
                                                              (self.gene_coords_4['contig'] == contig)].iloc[
                    -1, 4] > end:
                    window_sum_4 = self.gene_coords_4[(self.gene_coords_4['start'] < start) &
                                                      (self.gene_coords_4['contig'] == contig)].iloc[-1, 5]
            except IndexError:
                pass
            try:
                if window_sum_15 == 0.0 and self.gene_coords_15[(self.gene_coords_15['start'] < start) &
                                                                (self.gene_coords_15['contig'] == contig)].iloc[
                    -1, 4] > end:
                    window_sum_15 = self.gene_coords_15[(self.gene_coords_15['start'] < start) &
                                                        (self.gene_coords_15['contig'] == contig)].iloc[-1, 5]
            except IndexError:
                pass
            try:
                if window_sum_101 == 0.0 and self.gene_coords_101[(self.gene_coords_101['start'] < start) &
                                                                  (self.gene_coords_101['contig'] == contig)].iloc[
                    -1, 4] > end:
                    window_sum_101 = self.gene_coords_101[(self.gene_coords_101['start'] < start) &
                                                          (self.gene_coords_101['contig'] == contig)].iloc[-1, 5]
            except IndexError:
                pass
            # determine highest scoring model
            window_label = np.argmax([window_sum_11, window_sum_4, window_sum_15, window_sum_101])
            # keep logodd scores for each model and wsndow
            contig_logodd_scores.append((window_sum_11, window_sum_4, window_sum_15, window_sum_101))
            # save window label
            contig_labels.append(window_label)
            # update parameters
            start += self.stride
            end += self.stride
            # go to next contig
            if start > length_seq and contigs_seen < self.nr_contigs - 1:
                if len(contig_labels) == 0:
                    contig_labels.append(0)
                    contig_logodd_scores.append((0, 0, 0, 0))
                if length_seq > self.window_size and len(contig_labels) == 1:
                    contig_labels = contig_labels * int(np.ceil((length_seq - self.window_size) / self.stride) + 1)
                    contig_logodd_scores = contig_logodd_scores * int(
                        np.ceil((length_seq - self.window_size) / self.stride) + 1)
                labels.append(np.array(contig_labels))
                logodd_scores.append(np.array(contig_logodd_scores))
                if length_seq < 5000:
                    short_contigs.append(contig)
                contigs_seen += 1
                if contigs_seen > self.nr_contigs - 1:
                    break
                contig = self.sequence[contigs_seen].id
                length_seq = len(self.sequence[contigs_seen].seq)

                start = 0
                end = self.window_size
                contig_labels = []
                contig_logodd_scores = []

            elif start > length_seq and contigs_seen >= self.nr_contigs - 1:
                if len(contig_labels) == 0:
                    contig_labels.append(0)
                    contig_logodd_scores.append((0, 0, 0, 0))
                if length_seq > self.window_size and len(contig_labels) == 1:
                    contig_labels = contig_labels * int(np.ceil((length_seq - self.window_size) / self.stride) + 1)
                    contig_logodd_scores = contig_logodd_scores * int(
                        np.ceil((length_seq - self.window_size) / self.stride) + 1)
                labels.append(np.array(contig_labels))
                logodd_scores.append(np.array(contig_logodd_scores))
        if len(short_contigs) >= 1 and self.short_contigs:
            warnings.warn(f"Contig {contig} in {self.species} is shorter than 5000bp. "
                          f"Classification results may be not as reliable")
        if len(short_contigs_not_classified) >= 1:
            warnings.warn(
                f"Contig(s) {', '.join(short_contigs_not_classified)} in {self.species} are shorter than 5000bp. "
                f"Genetic code has been set to 11. Use --short_contigs to predict genetic code of these contigs")

        return np.array(labels), np.array(logodd_scores)

    def segment_genome_based_on_window_labels(self):
        """
        Determines coordinates of switch points
        :return: sequence of label
        """
        logging.info(f"Segmenting genome of {self.species}")
        if self.verbose:
            print(f"Segmenting genome of {self.species}")
        # iterate over all contigs
        confident_1 = {}
        confident_2 = {}
        segment_labels = {}
        n = 0
        # non_canonical_contigs = {1: [], 2: [], 3: []}
        # contigs_with_switch = {}
        while n <= len(self.sequence) - 1:
            current_contig = self.sequence[n].id
            current_genome_labels = self.updated_genome_labels[n]
            reference = current_genome_labels[0]
            start = 1
            segment_start = 1
            conf_1 = []
            conf_2 = []
            labels = []
            label_1 = np.unique(current_genome_labels)[0]
            try:
                label_2 = np.unique(current_genome_labels)[1]
            except IndexError:
                pass
            i = 1
            while i <= len(current_genome_labels) - 1:
                # uniform genetic code
                if reference == current_genome_labels[i]:
                    i += 1
                else:
                    start = (i - 1) * self.stride
                    end = i * self.stride + self.window_size
                    # previous genetic code
                    original_coords = self.label_to_gene_coords[current_genome_labels[i - 1]][
                        self.label_to_gene_coords[current_genome_labels[i - 1]].contig == current_contig]
                    # previous genetic code
                    new_coords = self.label_to_gene_coords[current_genome_labels[i]][
                        self.label_to_gene_coords[current_genome_labels[i]].contig == current_contig]
                    # strand values in previous and current window
                    strand_values = original_coords[(original_coords.start >= start)
                                                    & (original_coords.start < end)].strand.values
                    # stop coordinates in previous and current window
                    end_values = original_coords[(original_coords.start >= start)
                                                 & (original_coords.start < end)].end.values
                    if end_values.shape[0] == 0:
                        # center of two windows between switch occurred
                        end_values = np.array([(start + end) / 2])
                    # switch of encoding strand
                    if np.unique(strand_values).shape[0] > 1:
                        switch1 = np.where(np.diff(np.where(strand_values == '+', 1, 0)) != 0)[0][0]
                        switch_point1 = end_values[switch1]
                    # no siwtch of encoding strand
                    else:
                        try:
                            switch_point1 = \
                                original_coords.end.values[original_coords.start.values <= end_values[0]][-1]
                        except IndexError:
                            switch_point1 = end_values[0]
                    try:
                        switch_point2 = new_coords.start.values[new_coords.start.values >= switch_point1][0]
                    except IndexError:
                        switch_point2 = switch_point1 + 2
                    segment_end = int((switch_point1 + switch_point2) / 2)
                    if reference == label_1:
                        conf_1.append((segment_start, segment_end))
                        labels.append(self.label_to_gcode[label_1])
                    else:
                        conf_2.append((segment_start, segment_end))
                        labels.append(self.label_to_gcode[label_2])
                    logging.info(f"Contig {current_contig} in {self.species} has genetic code "
                                 f"{self.label_to_gcode[reference]} from {segment_start} to {segment_end}")
                    if self.verbose:
                        print(
                            f"Contig {current_contig} in {self.species} has genetic code "
                            f"{self.label_to_gcode[reference]} from {segment_start} to {segment_end}")
                    segment_start = segment_end + 1
                    reference = current_genome_labels[i]
                    i += 1
            segment_end = len(self.sequence[n])
            if reference == label_1:
                conf_1.append((segment_start, segment_end))
                labels.append(self.label_to_gcode[label_1])
            else:
                conf_2.append((segment_start, segment_end))
                labels.append(self.label_to_gcode[label_2])
            confident_1[current_contig] = conf_1
            confident_2[current_contig] = conf_2
            segment_labels[current_contig] = labels
            logging.info(f"Contig {current_contig} is {self.label_to_gcode[reference]} from {segment_start} to {segment_end}")
            if self.verbose:
                print(
                    f"Contig {current_contig} is {self.label_to_gcode[reference]} from {segment_start} to {segment_end}")
            n += 1
        return confident_1, confident_2, segment_labels

    def get_child_1(self, coords, start, end, gene_id, attributes, label):
        """
        Extract isoform 1 and format attributes
        :param coords: df, gene set 1
        :param start: int, start coordinates
        :param end: int, stop coordinates
        :param gene_id: int, current gene id
        :param attributes: str, current gene attributes
        :param label: int, label of corresponding genetic code
        :return: pd.Series, format isoform 1
        """
        child_1 = coords[(coords.start == start) & (coords.end == end)].copy()
        child_1.reset_index(drop=True, inplace=True)
        child_1.loc[0, 'attributes'] = 'gene_id {}.{}; Parent={};'.format(gene_id, self.label_to_gcode[label],
                                                                          gene_id) + ';'.join(attributes.split(';')[1:])
        child_1.loc[0, 'gene_label'] = 1
        return child_1

    def format_child_2(self, child_2, gene_id, label):
        """
        Format isoform 2
        :param child_2: pd.Series, isoform 2
        :param gene_id: int, current gene id
        :param label: int, label corresponding to genetic code of gene set 2
        :return: pd.Series, formatted isoform 2
        """
        child_2.reset_index(drop=True, inplace=True)

        child_2.loc[0, 'attributes'] = 'gene_id {}.{}; Parent={};'.format(gene_id,
                                                                          self.label_to_gcode[label],
                                                                          gene_id) + \
                                       ';'.join(child_2.attributes.values[0].split(';')[1:])
        if 'seen' in child_2.columns.values:
            child_2.drop('seen', axis=1, inplace=True)
        child_2.loc[0, 'gene_label'] = 2
        return child_2

    def format_parent(self, parent, gene_id, attributes, label_1, label_2, length=False, seq=None):
        """
        Format parent
        :param parent: pd.Series, parent gene
        :param gene_id: int, current gene id
        :param attributes: str, gene attributes
        :param label_1: int, label to genetic code 1
        :param label_2: int, label to genetic code 2
        :param length: boolean, whether to update length
        :param seq: Bio.Seq, required if length to be updated --> needed for updating GC% too
        :return: pd.Series, formatted parent
        """
        parent.reset_index(drop=True, inplace=True)
        parent.loc[0, 'source'] = parent.source.values[0].split('_')[0] + '_' + \
                                  '_'.join(sorted(map(str, [self.label_to_gcode[label_1],
                                                            self.label_to_gcode[label_2]])))
        if not length:
            parent.loc[0, 'attributes'] = 'gene_id {};'.format(gene_id) + ';'.join(attributes.split(';')[1:])
        else:
            length = parent.end.values[0] - parent.start.values[0] + 1
            parent.loc[0, 'length'] = length
            gc = int(round(GC(seq[parent.start.values[0] - 1: parent.end.values[0]]), 0))
            parent.loc[0, 'attributes'] = 'gene_id {}; gene_type bacteria; gc {}; length {};'.format(gene_id,
                                                                                                     gc, length)
        parent.loc[0, 'gene_label'] = 3

        return parent

    def compare_sets_of_predicted_genes(self):
        """
        Compare gene sets of two most common genetic codes to identify isoforms
        :return: pd.DataFrame  plain annotations with most likely genetic code
                 dict mapping of labels to genetic code per contig
                 pd.DataFrame annotations with isoforms
        """
        logging.info(f"Comparing sets of predicted genes in {self.species}")
        if self.verbose:
            print(f"Comparing sets of predicted genes in {self.species}")
        n = 0
        gene_labels_to_gcode = {}
        compared_gene_set = pd.DataFrame(columns=self.label_to_gene_coords[0].columns.values, dtype=object)
        compared_gene_set_w_isoforms = pd.DataFrame(columns=self.label_to_gene_coords[0].columns.values, dtype=object)
        while n <= len(self.updated_genome_labels) - 1:
            if np.unique(self.updated_genome_labels[n]).shape[0] == 1:
                coords = self.label_to_gene_coords[self.updated_genome_labels[n][0]].copy()
                coords.loc[:, 'updated_gene_label'] = 1
                # compared_gene_set_w_isoforms = compared_gene_set_w_isoforms.append(coords,
                #                                                                   sort=True)
                compared_gene_set = compared_gene_set.append(coords, sort=True)
                gene_labels_to_gcode[self.sequence[n].id] = {1: self.label_to_gcode[self.updated_genome_labels[n][0]]}
                n += 1
                continue
            contig = self.sequence[n].id
            current_gene_set = pd.DataFrame(columns=compared_gene_set.columns.values, dtype=object)
            current_gene_set_w_isoforms = pd.DataFrame(columns=compared_gene_set_w_isoforms.columns.values, dtype=object)
            label_1 = np.unique(self.updated_genome_labels[n])[0]
            coords_1 = self.label_to_gene_coords[label_1].copy()
            coords_1 = coords_1[coords_1.contig == contig]
            label_2 = np.unique(self.updated_genome_labels[n])[1]
            coords_2 = self.label_to_gene_coords[label_2].copy()
            coords_2 = coords_2[coords_2.contig == contig]
            logging.info("Contig {} in {} carries genetic code {} and {}".format(contig, self.species,
                                                                                 self.label_to_gcode[label_1],
                                                                                 self.label_to_gcode[label_2]))
            if self.verbose:
                print("Contig {} in {} carries genetic code {} and {}".format(contig, self.species,
                                                                              self.label_to_gcode[label_1],
                                                                              self.label_to_gcode[label_2]))
            gene_labels_to_gcode[contig] = {1: str(self.label_to_gcode[label_1]), 2: str(self.label_to_gcode[label_2]),
                                            3: str(self.label_to_gcode[label_1]) + '_' + str(
                                                self.label_to_gcode[label_2])}

            coords_2.loc[:, 'seen'] = 0
            previous_end = 0
            gene_labels = []
            sources = []
            gene_id = 1
            for start, end, strand, logodd, attributes in zip(coords_1.start.values, coords_1.end.values,
                                                              coords_1.strand.values,
                                                              coords_1.logodd.values, coords_1.attributes.values):
                if end <= previous_end:
                    continue
                # identical gene
                if coords_2[(coords_2.start == start) & (coords_2.end == end) & (coords_2.strand == strand)].shape[
                    0] == 1:
                    gene_labels.append(3)
                    sources.append('GeneMark.hmm2_' + gene_labels_to_gcode[contig][3])
                    logodd_2 = coords_2.loc[(coords_2.start == start) & (coords_2.end == end) & (
                                coords_2.strand == strand), 'logodd'].values[0]
                    coords_2.loc[
                        (coords_2.start == start) & (coords_2.end == end) & (coords_2.strand == strand), 'seen'] = 1

                    # append gene with greater logodd
                    if logodd >= logodd_2:
                        current_gene_set = current_gene_set.append(coords_1[(coords_1.start == start) &
                                                                            (coords_1.end == end)], ignore_index=True,
                                                                   sort=True)
                    elif logodd < logodd_2:
                        current_gene_set = current_gene_set.append(coords_2[(coords_2.start == start) &
                                                                            (coords_2.end == end) & (
                                                                                        coords_2.strand == strand)].drop(
                            'seen', axis=1),
                                                                   ignore_index=True, sort=True)

                    # get isoforms
                    child_1 = self.get_child_1(coords_1, start, end, gene_id, attributes, label_1)
                    child_2 = coords_2[
                        (coords_2.start == start) & (coords_2.end == end) & (coords_2.strand == strand)].copy()
                    child_2 = self.format_child_2(child_2, gene_id, label_2)
                    parent = coords_1[(coords_1.start == start) & (coords_1.end == end)].copy()
                    parent = self.format_parent(parent, gene_id, attributes, label_1, label_2)
                    parent.loc[0, 'logodd'] = max([child_1.logodd.values[0], child_2.logodd.values[0]])
                    current_gene_set_w_isoforms = current_gene_set_w_isoforms.append(parent,
                                                                                     ignore_index=True, sort=True)

                # identical stop different start --> treat as identical
                elif coords_2[(coords_2.end == end) & (coords_2.strand == strand)].shape[0] == 1:
                    gene_labels.append(3)
                    sources.append('GeneMark.hmm2_' + gene_labels_to_gcode[contig][3])
                    logodd_2 = coords_2.loc[(coords_2.end == end) & (coords_2.strand == strand), 'logodd'].values[0]
                    coords_2.loc[(coords_2.end == end) & (coords_2.strand == strand), 'seen'] = 1

                    # append gene with greater logodd
                    if logodd >= logodd_2:
                        current_gene_set = current_gene_set.append(coords_1[(coords_1.start == start) &
                                                                            (coords_1.end == end)], ignore_index=True,
                                                                   sort=True)
                    elif logodd < logodd_2:
                        current_gene_set = current_gene_set.append(
                            coords_2[(coords_2.end == end) & (coords_2.strand == strand)].drop('seen', axis=1),
                            ignore_index=True, sort=True)

                    # get isoforms
                    child_1 = self.get_child_1(coords_1, start, end, gene_id, attributes, label_1)

                    child_2 = coords_2[(coords_2.end == end) & (coords_2.strand == strand)].copy()
                    child_2 = self.format_child_2(child_2, gene_id, label_2)

                    parent = coords_1[(coords_1.start == start) & (coords_1.end == end)].copy()
                    parent.reset_index(drop=True, inplace=True)
                    parent.loc[0, 'start'] = min([child_1.start.values[0], child_2.start.values[0]])
                    parent = self.format_parent(parent, gene_id, attributes, label_1, label_2, length=True,
                                                seq=self.sequence[n].seq)
                    parent.loc[0, 'logodd'] = '.'

                    current_gene_set_w_isoforms = current_gene_set_w_isoforms.append([parent, child_1, child_2],
                                                                                     ignore_index=True, sort=True)

                # identical start different stop and stop coordinates differ by less than thresh --> treat as isoforms
                elif coords_2[(coords_2.start == start) & (coords_2.strand == strand)].shape[0] == 1 and np.absolute(
                        coords_2.loc[(coords_2.start == start) & (coords_2.strand == strand), 'end'].values[
                            0] - end) <= self.tolerance:
                    gene_labels.append(3)
                    sources.append('GeneMark.hmm2_' + gene_labels_to_gcode[contig][3])
                    logodd_2 = coords_2.loc[(coords_2.start == start) & (coords_2.strand == strand), 'logodd'].values[0]
                    coords_2.loc[(coords_2.start == start) & (coords_2.strand == strand), 'seen'] = 1

                    # append gene with greater logodd
                    if logodd >= logodd_2:
                        current_gene_set = current_gene_set.append(coords_1[(coords_1.start == start) &
                                                                            (coords_1.end == end)], ignore_index=True,
                                                                   sort=True)
                    elif logodd < logodd_2:
                        current_gene_set = current_gene_set.append(
                            coords_2[(coords_2.start == start) & (coords_2.strand == strand)].drop('seen', axis=1),
                            ignore_index=True, sort=True)
                        previous_end = \
                        coords_2.loc[(coords_2.start == start) & (coords_2.strand == strand), 'end'].values[0]

                    # get isoforms
                    child_1 = self.get_child_1(coords_1, start, end, gene_id, attributes, label_1)

                    child_2 = coords_2[(coords_2.start == start) & (coords_2.strand == strand)].copy()
                    child_2 = self.format_child_2(child_2, gene_id, label_2)

                    parent = coords_1[(coords_1.start == start) & (coords_1.end == end)].copy()
                    parent.reset_index(drop=True, inplace=True)
                    parent.loc[0, 'end'] = max([child_1.end.values[0], child_2.end.values[0]])
                    parent = self.format_parent(parent, gene_id, attributes, label_1, label_2, length=True,
                                                seq=self.sequence[n].seq)
                    parent.loc[0, 'logodd'] = '.'
                    current_gene_set_w_isoforms = current_gene_set_w_isoforms.append([parent, child_1, child_2],
                                                                                     ignore_index=True, sort=True)

                # gene from second set overspanning area
                elif coords_2[(coords_2.start < start) & (coords_2.end > end) & (coords_2.strand == strand)].shape[
                    0] == 1:
                    overspanning_gene = coords_2.loc[(coords_2.start < start) &
                                                     (coords_2.end > end) & (coords_2.strand == strand), ['start',
                                                                                                          'end',
                                                                                                          'logodd']].copy()
                    coords_2.loc[
                        (coords_2.start < start) & (coords_2.end > end) & (coords_2.strand == strand), 'seen'] = 1
                    # compare logodd scores
                    logodd_1 = coords_1.loc[(coords_1.start >= overspanning_gene.start.values[0]) &
                                            (coords_1.end < overspanning_gene.end.values[0]), 'logodd'].sum()

                    # decide if isoform
                    isoforms = False
                    if start - overspanning_gene.start.values[0] < self.tolerance and overspanning_gene.end.values[
                        0] - end < self.tolerance:
                        parent = coords_2[
                            (coords_2.start < start) & (coords_2.end > end) & (coords_2.strand == strand)].copy()
                        parent = self.format_parent(parent, gene_id, parent.attributes.values[0], label_1, label_2)
                        parent.loc[0, 'logodd'] = '.'
                        child_1 = self.get_child_1(coords_1, start, end, gene_id, attributes, label_1)

                        child_2 = coords_2[
                            (coords_2.start < start) & (coords_2.end > end) & (coords_2.strand == strand)].copy()
                        child_2 = self.format_child_2(child_2, gene_id, label_2)

                        current_gene_set_w_isoforms = current_gene_set_w_isoforms.append([parent, child_1, child_2],
                                                                                         ignore_index=True, sort=True)
                        isoforms = True

                    # gene set 2 wins
                    if overspanning_gene.logodd.values[0] > logodd_1:
                        gene_labels.append(2)
                        sources.append('GeneMark.hmm2_' + gene_labels_to_gcode[contig][2])
                        previous_end = overspanning_gene.end.values[0]
                        current_gene_set = current_gene_set.append(coords_2[(coords_2.start < start) &
                                                                            (coords_2.end > end) & (
                                                                                        coords_2.strand == strand)].drop(
                            'seen', axis=1),
                                                                   ignore_index=True, sort=True)
                        if not isoforms:
                            gene = coords_2[
                                (coords_2.start < start) & (coords_2.end > end) & (coords_2.strand == strand)].copy()
                            gene.reset_index(drop=True, inplace=True)
                            gene.drop('seen', axis=1, inplace=True)
                            gene.loc[0, 'attributes'] = 'gene_id {};'.format(gene_id) + ';'.join(
                                gene.attributes.values[0].split(
                                    ';')[1:])
                            gene.loc[0, 'gene_label'] = 2
                            current_gene_set_w_isoforms = current_gene_set_w_isoforms.append(gene,
                                                                                             ignore_index=True,
                                                                                             sort=True)

                    # gene set 1 wins
                    elif overspanning_gene.logodd.values[0] < logodd_1:
                        gene_labels.extend([1] * coords_1.loc[(coords_1.start >= overspanning_gene.start.values[0]) &
                                                              (coords_1.end < overspanning_gene.end.values[0])].shape[
                            0])
                        sources.extend(['GeneMark.hmm2_' + gene_labels_to_gcode[contig][1]] *
                                       coords_1.loc[(coords_1.start >= overspanning_gene.start.values[0]) &
                                                    (coords_1.end < overspanning_gene.end.values[0])].shape[0])
                        previous_end = coords_1.loc[(coords_1.start >= overspanning_gene.start.values[0]) &
                                                    (coords_1.end < overspanning_gene.end.values[0]), 'end'].max()
                        current_gene_set = current_gene_set.append(
                            coords_1.loc[(coords_1.start >= overspanning_gene.start.values[0]) &
                                         (coords_1.end < overspanning_gene.end.values[0])], ignore_index=True,
                            sort=True)
                        if not isoforms:
                            gene = coords_1.loc[(coords_1.start >= overspanning_gene.start.values[0]) &
                                                (coords_1.end < overspanning_gene.end.values[0])].copy()
                            gene.reset_index(drop=True, inplace=True)
                            gene.loc[0, 'attributes'] = 'gene_id {};'.format(gene_id) + ';'.join(
                                gene.attributes.values[0].split(
                                    ';')[1:])
                            gene.loc[0, 'gene_label'] = 1
                            current_gene_set_w_isoforms = current_gene_set_w_isoforms.append(gene, ignore_index=True,
                                                                                             sort=True)
                    # undecided
                    elif overspanning_gene.logodd.values[0] == logodd_1:
                        gene_labels.append(3)
                        sources.append('GeneMark.hmm2_' + gene_labels_to_gcode[contig][3])

                        current_gene_set = current_gene_set.append(coords_2[(coords_2.start < start) &
                                                                            (coords_2.end > end) & (
                                                                                        coords_2.strand == strand)].drop(
                            'seen', axis=1),
                                                                   ignore_index=True, sort=True)
                        if not isoforms:
                            gene = coords_2[
                                (coords_2.start < start) & (coords_2.end > end) & (coords_2.strand == strand)].copy()
                            gene.reset_index(drop=True, inplace=True)
                            gene.drop('seen', axis=1, inplace=True)
                            gene.loc[0, 'attributes'] = 'gene_id {};'.format(gene_id) + ';'.join(
                                gene.attributes.values[0].split(
                                    ';')[1:])
                            gene.loc[0, 'gene_label'] = 2
                            current_gene_set_w_isoforms = current_gene_set_w_isoforms.append(gene,
                                                                                             ignore_index=True,
                                                                                             sort=True)

                # gene from set 1 overspans genes from set 2
                elif coords_2[(coords_2.start >= start) & (coords_2.end < end) & (coords_2.strand == strand)].shape[
                    0] >= 1:
                    logodd_2 = coords_2.loc[(coords_2.start >= start) & \
                                            (coords_2.end < end) & (coords_2.strand == strand), 'logodd'].sum()
                    coords_2.loc[
                        (coords_2.start >= start) & (coords_2.end < end) & (coords_2.strand == strand), 'seen'] = 1
                    overspanned_genes = coords_2[
                        (coords_2.start >= start) & (coords_2.end < end) & (coords_2.strand == strand)]

                    # decide if isoforms
                    isoforms = False
                    if overspanned_genes.shape[0] == 1 and overspanned_genes.start.values - start < self.tolerance and \
                            end - overspanned_genes.end.values[0] < self.tolerance:
                        parent = coords_1[(coords_1.start == start) & (coords_1.end == end)].copy()
                        parent.reset_index(drop=True, inplace=True)
                        parent.loc[0, 'source'] = parent.source.values[0].split('_')[0] + '_' + \
                                                  '_'.join(sorted(map(str, [self.label_to_gcode[label_1],
                                                                            self.label_to_gcode[label_2]])))
                        parent.loc[0, 'logodd'] = '.'
                        child_1 = self.get_child_1(coords_1, start, end, gene_id, attributes, label_1)

                        child_2 = overspanned_genes.copy()
                        child_2 = self.format_child_2(child_2, gene_id, label_2)

                        current_gene_set_w_isoforms = current_gene_set_w_isoforms.append([parent, child_1, child_2],
                                                                                         ignore_index=True, sort=True)
                        isoforms = True

                    # gene set 1 wins
                    if logodd > logodd_2:
                        gene_labels.append(1)
                        sources.append('GeneMark.hmm2_' + gene_labels_to_gcode[contig][1])
                        current_gene_set = current_gene_set.append(coords_1[(coords_1.start == start) &
                                                                            (coords_1.end == end)], ignore_index=True,
                                                                   sort=True)
                        if not isoforms:
                            gene = coords_1[(coords_1.start == start) & (coords_1.end == end)].copy()
                            gene.reset_index(drop=True, inplace=True)
                            gene.loc[0, 'attributes'] = 'gene_id {};'.format(gene_id) + ';'.join(
                                gene.attributes.values[0].split(
                                    ';')[1:])
                            gene.loc[0, 'gene_label'] = 1
                            current_gene_set_w_isoforms = current_gene_set_w_isoforms.append(gene,
                                                                                             ignore_index=True,
                                                                                             sort=True)
                    # gene set 2 wins
                    elif logodd < logodd_2:
                        gene_labels.extend([2] * coords_2.loc[(coords_2.start >= start) & \
                                                              (coords_2.end < end) & (coords_2.strand == strand)].shape[
                            0])
                        sources.extend(['GeneMark.hmm2_' + gene_labels_to_gcode[contig][2]] *
                                       coords_2.loc[(coords_2.start >= start) & \
                                                    (coords_2.end < end) & (coords_2.strand == strand)].shape[0])
                        current_gene_set = current_gene_set.append(coords_2.loc[(coords_2.start >= start) &
                                                                                (coords_2.end < end) & (
                                                                                            coords_2.strand == strand)].drop(
                            'seen',
                            axis=1),
                                                                   ignore_index=True, sort=True)
                        if not isoforms:
                            gene = coords_2.loc[
                                (coords_2.start >= start) & (coords_2.end < end) & (coords_2.strand == strand)].copy()
                            gene.reset_index(drop=True, inplace=True)
                            gene.drop('seen', axis=1, inplace=True)
                            gene.loc[0, 'gene_label'] = 2
                            gene.loc[0, 'attributes'] = 'gene_id {};'.format(gene_id) + ';'.join(
                                gene.attributes.values[0].split(
                                    ';')[1:])
                            current_gene_set_w_isoforms = current_gene_set_w_isoforms.append(gene,
                                                                                             ignore_index=True,
                                                                                             sort=True)
                    # undecided
                    elif logodd == logodd_2:
                        gene_labels.append(3)
                        sources.append('GeneMark.hmm2_' + gene_labels_to_gcode[contig][3])
                        current_gene_set = current_gene_set.append(coords_1[(coords_1.start == start) &
                                                                            (coords_1.end == end)], ignore_index=True,
                                                                   sort=True)
                        if not isoforms:
                            gene = coords_1[(coords_1.start == start) & (coords_1.end == end)].copy()
                            gene.reset_index(drop=True, inplace=True)
                            gene.loc[0, 'gene_label'] = 1
                            gene.loc[0, 'attributes'] = 'gene_id {};'.format(gene_id) + ';'.join(
                                gene.attributes.values[0].split(
                                    ';')[1:])
                            current_gene_set_w_isoforms = current_gene_set_w_isoforms.append(gene,
                                                                                             ignore_index=True,
                                                                                             sort=True)
                else:
                    gene_labels.append(1)
                    sources.append('GeneMark.hmm2_' + gene_labels_to_gcode[contig][1])
                    current_gene_set = current_gene_set.append(coords_1[(coords_1.start == start) &
                                                                        (coords_1.end == end)], ignore_index=True,
                                                               sort=True)
                    gene = coords_1[(coords_1.start == start) & (coords_1.end == end)].copy()
                    gene.reset_index(drop=True, inplace=True)
                    gene.loc[0, 'attributes'] = 'gene_id {};'.format(gene_id) + ';'.join(
                        gene.attributes.values[0].split(';')[1:])
                    gene.loc[0, 'gene_label'] = 1
                    current_gene_set_w_isoforms = current_gene_set_w_isoforms.append(gene,
                                                                                     ignore_index=True,
                                                                                     sort=True)
                gene_id += 1

            current_gene_set.loc[current_gene_set.contig == contig, 'gene_label'] = gene_labels
            current_gene_set.loc[current_gene_set.contig == contig, 'source'] = sources

            if np.any(coords_2.seen.values == 0):
                unseen_genes = coords_2.loc[coords_2.seen == 0].drop('seen', axis=1)
                unseen_genes.loc[:, 'gene_label'] = 2
                current_gene_set = current_gene_set.append(unseen_genes, ignore_index=True, sort=True)
                current_gene_set_w_isoforms = self.insert_unseen_genes(current_gene_set_w_isoforms, unseen_genes)

            current_gene_set.sort_values('start', inplace=True)
            current_gene_set.reset_index(inplace=True, drop=True)
            updated_gene_labels = self.find_best_sequence_of_labels(current_gene_set.gene_label.values,
                                                                    genome_label=None,
                                                                    consecutive_labels=self.consecutive_gene_labels,
                                                                    allow_multiple_labels=True)
            # updated_gene_labels = [gene_labels_to_gcode[x] for x in updated_gene_labels]
            current_gene_set.loc[:, 'updated_gene_label'] = updated_gene_labels
            compared_gene_set = compared_gene_set.append(current_gene_set, sort=True)
            compared_gene_set_w_isoforms = compared_gene_set_w_isoforms.append(current_gene_set_w_isoforms, sort=True)
            n += 1
            logging.info("Identified {} possible isoforms in {} in {}".format(
                    current_gene_set_w_isoforms[current_gene_set_w_isoforms.logodd == '.'].shape[0], contig,
                self.species))
            if self.verbose:
                print("Identified {} possible isoforms in {} in {}".format(
                    current_gene_set_w_isoforms[current_gene_set_w_isoforms.logodd == '.'].shape[0], contig,
                    self.species))
        compared_gene_set_w_isoforms.reset_index(drop=True, inplace=True)
        # compared_gene_set_w_isoforms.updated_gene_label.astype(int)
        compared_gene_set.updated_gene_label.astype(int)
        return compared_gene_set, gene_labels_to_gcode, compared_gene_set_w_isoforms

    def search_for_change_in_coding_strand_in_vicinity(self, coords, i):
        """
        Search for change of protein encoding strand in vicinity to switch region of identified switch
        :param coords: pd.DataFrame, current gene prediction
        :param i: int, index at which switch happens
        """
        # look at strand of 5 surrounding genes
        strand_values = coords.strand[i - 5: i + 5]
        if np.unique(strand_values).shape[0] == 1:
            ind = i
        else:
            switch = np.where(np.diff(np.where(strand_values == '+', 1, 0)) != 0)[0][0]
            ind = switch - 5 + i
        return ind

    # def segmentate_genome_based_on_gene_labels(self):
    #     """
    #     Segmentate genome based on gene labels
    #     :return: confident_1 -> list of tuples with (start, end) coordinates of segments with genetic code label 1
    #              confident_2 -> list of tuples with (start, end) coordinates of segments with genetic code label 2
    #              ambiguous -> list of tuples with (start, end) coordinates of segments with genetic code label 1 & 2
    #     """
    #     confident_1 = {}
    #     confident_2 = {}
    #     ambiguous = {}
    #     segment_labels = {}
    #     n = 0
    #     while n <= len(self.sequence) - 1:
    #         current_contig = self.sequence[n].id
    #         current_coords = self.coords[self.coords.contig == current_contig].copy()
    #         current_coords.reset_index(inplace=True, drop=True)
    #         # set reference
    #         reference = current_coords.loc[0, 'updated_gene_label']
    #         start = 1
    #         conf_1 = []
    #         conf_2 = []
    #         amb = []
    #         labels = []
    #         # iterate over all genes
    #         i = 1
    #         while i < current_coords.shape[0]:
    #             # same as reference --> continue
    #             if reference == current_coords.loc[i, 'updated_gene_label']:
    #                 i += 1
    #             # end of segment --> save boundaries
    #             else:
    #                 ind = self.search_for_change_in_coding_strand_in_vicinity(current_coords, i)
    #                 if current_coords.loc[ind - 1, 'start'] == current_coords.loc[ind, 'start']:
    #                     end = self.coords_w_isoforms.loc[(self.coords_w_isoforms.contig == current_contig) &
    #                                                      (self.coords_w_isoforms.start == current_coords.loc[
    #                                                          ind, 'start']),
    #                                                      'end'].max()
    #
    #                 else:
    #                     end = self.coords_w_isoforms.loc[(self.coords_w_isoforms.contig == current_contig) &
    #                                                      (self.coords_w_isoforms.start == current_coords.loc[
    #                                                          ind - 1, 'start']),
    #                                                      'end'].max()
    #
    #                 if reference == 1:
    #                     conf_1.append((start, end))
    #                     labels.append(self.gene_label_to_gcode[current_contig][1])
    #                 elif reference == 2:
    #                     conf_2.append((start, end))
    #                     labels.append(self.gene_label_to_gcode[current_contig][2])
    #                 elif reference == 3:
    #                     amb.append((start, end))
    #                     labels.append(self.gene_label_to_gcode[current_contig][3])
    #                 logging.info('Contig {} in {} has genetic code {} from {} to {}'.format(current_contig,
    #                                                                                         self.species,
    #                                                                                         self.gene_label_to_gcode[
    #                                                                                         current_contig][reference],
    #                                                                                         start, end))
    #                 if self.verbose:
    #                     print('Contig {} in {} has genetic code {} from {} to {}'.format(current_contig, self.species,
    #                                                                                      self.gene_label_to_gcode[
    #                                                                                      current_contig][reference],
    #                                                                                      start, end))
    #                 while current_coords.loc[ind, 'start'] < end:
    #                     ind += 1
    #                 start = current_coords.loc[ind, 'start']
    #                 reference = current_coords.loc[i, 'updated_gene_label']
    #                 if ind > i + 1:
    #                     i = ind
    #                 else:
    #                     i += 1
    #         end = len(self.sequence[n])
    #
    #         if reference == 1:
    #             conf_1.append((start, end))
    #             labels.append(self.gene_label_to_gcode[current_contig][1])
    #         elif reference == 2:
    #             conf_2.append((start, end))
    #             labels.append(self.gene_label_to_gcode[current_contig][2])
    #         elif reference == 3:
    #             amb.append((start, end))
    #             labels.append(self.gene_label_to_gcode[current_contig][3])
    #         confident_1[current_contig] = conf_1
    #         confident_2[current_contig] = conf_2
    #         ambiguous[current_contig] = amb
    #         segment_labels[current_contig] = labels
    #         logging.info('Contig {} in {} has genetic code {} from {} to {}'.format(current_contig, self.species,
    #                                                                                 self.gene_label_to_gcode[current_contig][reference],
    #                                                                                 start, end))
    #         if self.verbose:
    #             print('Contig {} in {} has genetic code {} from {} to {}'.format(current_contig, self.species,
    #                                                                              self.gene_label_to_gcode[current_contig][reference],
    #                                                                              start, end))
    #         n += 1
    #     return confident_1, confident_2, ambiguous, segment_labels

    @staticmethod
    def find_best_sequence_of_labels(labels, genome_label, consecutive_labels,
                                     allow_multiple_labels=False):
        """
        Merges segments with one genetic code which are shorter than n labels.
        Extends the longest segment in each direction until it reaches the start/end or finds n consecutive labels with
        different label. This is done iteratively until all window labels have been checked
        :param labels: array-like, sequence of labels
        :param genome_label: array-like, window labels
        :param consecutive_labels: int, minimum number of consecutive genes with same label to start new segment
        :param allow_multiple_labels: boolean, whether two allow more than two genetic codes (needed for gene labels)
        :return: numpy array, updated labels
        """
        # intialize scoring matrix
        scoring_matrix = np.zeros_like(labels).astype(int)
        scoring_matrix[0] = 1
        updated_labels = labels.copy()
        updated = np.zeros_like(updated_labels, dtype=bool)
        # if more than two 2 genetic codes have been predicted, consider the least like ones as noise
        # --> identify which genetic code(s)
        if not allow_multiple_labels:
            likely_false_predictions = np.unique(genome_label)[np.argsort(np.unique(genome_label,
                                                                                    return_counts=True)[1])[::-1]][2:]
        else:
            likely_false_predictions = []
        i = 1
        # determine length of segments with same genetic code
        while i <= len(labels) - 1:
            # create scoring based on length of sequence of windows with same genetic code
            if labels[i] == labels[i - 1] or labels[i] in likely_false_predictions:
                scoring_matrix[i] = scoring_matrix[i - 1] + 1
            else:
                scoring_matrix[i] = 1
            i += 1
        # merge segments until all it has been looked at all windows once
        while not np.all(scoring_matrix == 0):
            # find longest unseen sequence of windows with same genetic code
            highest_scoring_seed = np.argmax(scoring_matrix)
            # determine borders of current segment
            start = highest_scoring_seed - (scoring_matrix[highest_scoring_seed])
            end = highest_scoring_seed + 1
            # determine label of current segment
            reference = updated_labels[highest_scoring_seed]
            # go to previous window if labels is among genetic codes previously identified as noise
            while reference in likely_false_predictions:
                highest_scoring_seed -= 1
                reference = updated_labels[highest_scoring_seed]
            # extend start until n mismatches in a row
            mismatch = 0
            while start >= 0 and mismatch < consecutive_labels:
                # mismatch
                if updated_labels[start] != reference:
                    mismatch += 1
                # no mismatch --> reset counter
                else:
                    mismatch = 0
                start -= 1
            # extend end  until n consecutive mismateches
            mismatch = 0
            while end <= updated_labels.shape[0] - 1 and mismatch < consecutive_labels:
                # mismatch
                if updated_labels[end] != reference:
                    mismatch += 1
                # no mismatch --> reset counter
                else:
                    mismatch = 0
                end += 1
            if start < 0:
                start = 0
            else:
                # plus n+1 because start is included and n mismatches
                start += consecutive_labels + 1
            if end > updated_labels.shape[0] - 1:
                end = updated_labels.shape[0]
            else:
                # minus n because end is excluded and n mismatches
                end -= consecutive_labels
            # update labels and mark as updated
            if end - start < consecutive_labels:
                length_left = 0
                i = start - 1
                ref_left = updated_labels[i]
                length_right = 0
                j = end
                try:
                    ref_right = updated_labels[j]
                except IndexError:
                    ref_right = ref_left
                if ref_left == ref_right:
                    reference = ref_left
                else:
                    while i >= 0 and updated_labels[i] == ref_left:
                        i -= 1
                        length_left += 1
                    while j < updated_labels.shape[0] and updated_labels[j] == ref_right:
                        j += 1
                        length_right += 1
                    if length_left > length_right:
                        reference = ref_left
                    elif length_right >= length_left:
                        reference = ref_right
            scoring_matrix[start: end] = 0
            updated_labels[start: end] = reference
            updated[start: end] = True
        return updated_labels

    def insert_unseen_genes(self, coords, unseen_genes):
        for i in range(unseen_genes.shape[0]):
            gene = unseen_genes.iloc[i].copy()
            coords_1 = coords[coords.start <= gene.start]
            try:
                ind = coords_1.index.values[-1]
                ind += 1
            except IndexError:
                gene.attributes = 'gene_id {};'.format(0) + \
                                  ';'.join(gene.attributes.split(';')[1:])
                coords_1 = self.update_gene_ids(coords.copy())
                coords = pd.DataFrame(gene).T.append(coords_1, ignore_index=True, sort=True)
                coords.reset_index(drop=True, inplace=True)
                continue
            # if ind >= coords.shape[0]:
            #    gene.attributes = 'gene_id {};'.format(int(coords.loc[coords.shape[0] - 1, 'attributes'].\
            #    split(';')[0].split(' ')[1].split('.')[0]) + 1) + \
            #                      ';'.join(gene.attributes.split(';')[1:])
            #    coords = coords.append(gene, ignore_index=True, sort=True)
            #    continue
            while ind <= coords.shape[0] - 1 and 'Parent' in coords.loc[ind, 'attributes']:
                ind += 1
            if ind >= coords.shape[0]:
                gene.attributes = 'gene_id {};'.format(
                    int(coords.loc[coords.shape[0] - 1, 'attributes'].split(';')[0].split(' ')[1].split('.')[0]) + 1) +\
                                  ';'.join(gene.attributes.split(';')[1:])
                coords = coords.append(gene, ignore_index=True, sort=True)
                continue
            else:
                gene.attributes = 'gene_id {};'.format(
                    int(coords.loc[ind, 'attributes'].split(';')[0].split(' ')[1].split('.')[0])) + \
                                  ';'.join(gene.attributes.split(';')[1:])
                coords_1 = coords.loc[:ind - 1].copy()
                coords_2 = self.update_gene_ids(coords.loc[ind:].copy())
                coords = coords_1.append([pd.DataFrame(gene).T, coords_2], ignore_index=True,
                                         sort=True)
                coords.reset_index(drop=True, inplace=True)

        return coords

    @staticmethod
    def update_gene_ids(coords):
        """
        Update gene ids --> increment each id by 1
        """
        coords.reset_index(drop=True, inplace=True)
        current_gene_id = int(coords.loc[0, 'attributes'].split(';')[0].split(' ')[1]) + 1
        reference = coords.loc[0, 'attributes'].split(';')[0].split(' ')[1]
        coords.loc[0, 'attributes'] = coords.loc[0, 'attributes'].replace('gene_id {}'.format(reference),
                                                                          'gene_id {}'.format(current_gene_id)).replace(
            'Parent={}'.format(reference),
            'Parent={}'.format(current_gene_id))
        for i in range(1, coords.shape[0]):
            c_attribute = coords.loc[i, 'attributes'].split(';')[0]
            if c_attribute.startswith('gene_id {}'.format(reference)):
                coords.loc[i, 'attributes'] = coords.loc[i, 'attributes'].replace('gene_id {}'.format(reference),
                                                                                  'gene_id {}'.format(
                                                                                      current_gene_id)).replace(
                    'Parent={}'.format(reference),
                    'Parent={}'.format(current_gene_id))
            else:
                current_gene_id += 1
                reference = coords.loc[i, 'attributes'].split(';')[0].split(' ')[1]
                coords.loc[i, 'attributes'] = coords.loc[i, 'attributes'].replace('gene_id {}'.format(reference),
                                                                                  'gene_id {}'.format(
                                                                                      current_gene_id)).replace(
                    'Parent={}'.format(reference),
                    'Parent={}'.format(current_gene_id))
        return coords
