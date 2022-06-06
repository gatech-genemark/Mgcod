#!/usr/bin/env python
# -------------------------------------------------
# AlternatingGeneticCodeAnnotator
# Aaron Pfennig
# Georgia Institute of Technology
# -------------------------------------------------
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Data.CodonTable import CodonTable
import subprocess
import pandas as pd
import numpy as np
import os
import warnings
import logging
warnings.simplefilter(action='ignore', category=FutureWarning)


class Results:
    def __init__(self, pipeline, species, path_to_predictions, path_to_genome, path_to_output, date,
                 delete, amino_acids, nucleotides, verbosity, working_dir):
        """
        Write results to file
        :param pipeline: Pipeline object
        :param species: str, species
        :param path_to_predictions: str, path to mgm predictions
        :param path_to_genome: str, path to input genome
        :param path_to_output: str, path to ouput
        :param date: datetime object, start time of pipeline
        :param delete: boolean, delete intermediary files
        :param amino_acids: boolean, extract amino acid sequences of predicted genes
        :param nucleotides: boolean, extract nucleotide sequences of predicted genes
        :param verbose: boolean, verbosity
        :param working_dir: str, current working_directory
        """
        self.verbose = verbosity
        self.cwd = working_dir
        self.results = pipeline
        self.species = species
        self.path_to_predictions = path_to_predictions
        self.path_to_final_annotations = path_to_output + species + '.gff'
        self.path_to_nt_seqs = '{}proteins_nt_{}.fasta'.format(path_to_output, species)
        self.path_to_aa_seqs = '{}proteins_aa_{}.fasta'.format(path_to_output, species)
        self.path_to_genome = path_to_genome
        self.path_to_output = path_to_output + species + ".tab"
        self.date = date
        self.delete = delete
        self.amino_acids = amino_acids
        self.nucleotides = nucleotides
        # extract sequences
        self.stop_recodings = self.extract_sequences()
        logging.info(f'Writing results for {self.species}')
        if self.verbose:
            print(f"Writing results for {self.species}")
        # save results
        self.save_results()
        # add number of recodings to gff
        attributes = self.results.coords_w_isoforms.attributes.values
        gene_ids = [float(attr.split(';')[0].split(' ')[1]) for attr in attributes]
        self.results.coords_w_isoforms['gene_id'] =  gene_ids
        for i in range(len(attributes)):
            if attributes[i].endswith(';'):
                pass
            else:
                attributes[i] += ';'
            if self.results.coords_w_isoforms.loc[i, 'logodd'] == '.':
                self.stop_recodings.insert(i, {'TAA': [], 'TGA': [], 'TAG': []})
                continue
            attributes[i] = attributes[i] + ' number of recodings: ' +\
                            f"total {len(self.stop_recodings[i]['TAA']) + len(self.stop_recodings[i]['TGA']) + len(self.stop_recodings[i]['TAG'])}, " +\
                            f"TAA {len(self.stop_recodings[i]['TAA'])}, " +\
                            f"TGA {len(self.stop_recodings[i]['TGA'])}, " +\
                            f"TAG {len(self.stop_recodings[i]['TAG'])};"
        self.results.coords_w_isoforms.attributes = attributes
        #self.results.coords_w_isoforms.sort_values(['sequence', 'gene_id', 'start'], ascending=[True, True, True], inplace=True)

    def __call__(self):
        self.write_stats()
        self.write_final_results()
        if self.delete:
            os.remove(f'{self.path_to_predictions}mgm_11_{self.species}.gff')
            os.remove(f'{self.path_to_predictions}mgm_4_{self.species}.gff')
            os.remove(f'{self.path_to_predictions}mgm_15_{self.species}.gff')
            os.remove(f'{self.path_to_predictions}mgm_101_{self.species}.gff')

    def save_results(self):
        """
        save classification results
        """
        if type(self.results.genome_label) == int:
            index = [(self.results.sequence[0].id, 'Major genetic code'),
                     (self.results.sequence[0].id, 'Total gene count')]
            results = pd.DataFrame(index=pd.MultiIndex.from_tuples(index))
            results.loc[(self.results.sequence[0].id, 'Major genetic code'), 'Sequence'] = self.results.genome_label
            results.loc[(self.results.sequence[0].id, 'Total gene count'), 'Sequence'] = self.results.coords_w_isoforms.shape[0]
        # multiple contigs or switch point
        else:
            results = pd.DataFrame()
            level_1 = [seq.id for seq in self.results.sequence]
            for i, contig in enumerate(level_1):
                gene_label_to_gcode = self.results.gene_label_to_gcode[contig]
                # only need to initialize index once
                if i == 0:
                    
                    level_2 = ['Coordinates in bp', 'Major genetic code', 'Nr genes with {}'.format(gene_label_to_gcode[1]),
                               'Nr identical genes', 'Nr isoforms',
                               'Total gene count']
                    index = pd.MultiIndex.from_product([level_1, level_2])
                    results = results.append(pd.DataFrame(index=index), sort=True)
                #segments = sorted([seg for seg_type in [self.results.ambiguous[contig],
                #                                        self.results.confident_1[contig],
                #                                        self.results.confident_2[contig]] for seg in seg_type])
                segments = sorted([seg for seg_type in [self.results.confident_1[contig],
                                                        self.results.confident_2[contig]] for seg in seg_type]) 
                coords = self.results.coords_w_isoforms[self.results.coords_w_isoforms.contig == contig]
                if len(gene_label_to_gcode.keys()) == 1:
                    results.loc[(contig, 'Coordinates in bp'), "Sequence"] = str(segments[0])
                    results.loc[(contig, 'Major genetic code'), "Sequence"] = gene_label_to_gcode[1]
                    results.loc[(contig, 'Total gene count'), "Sequence"] = coords.shape[0]
                    continue
                for n, segment in enumerate(segments, start=1):
                    results.loc[(contig, 'Coordinates in bp'), "Segment {}".format(n)] = str(segment)
                    #if segment in self.results.ambiguous[contig]:
                    #    results.loc[(contig, 'Major genetic code'), "Segment {}".format(n)] = gene_label_to_gcode[3]
                    if segment in self.results.confident_1[contig]:
                        results.loc[(contig, 'Major genetic code'), "Segment {}".format(n)] = gene_label_to_gcode[1]
                    elif segment in self.results.confident_2[contig]:
                        results.loc[(contig, 'Major genetic code'), "Segment {}".format(n)] = gene_label_to_gcode[2]
                    current_coords = coords[(coords.start >= segment[0]) & (coords.start <= segment[1])]
                    nr_isoforms = current_coords[current_coords.logodd == '.'].shape[0]
                    nr_confident_1 = current_coords[current_coords.source ==
                                                    'GeneMark.hmm2_{}'.format(gene_label_to_gcode[1])].shape[0]
                    nr_confident_2 = current_coords[current_coords.source ==
                                                    'GeneMark.hmm2_{}'.format(gene_label_to_gcode[2])].shape[0]
                    nr_identical = current_coords[current_coords.source ==
                                                  'GeneMark.hmm2_{}_{}'.format(*sorted([gene_label_to_gcode[1],
                                                                               gene_label_to_gcode[2]]))].shape[0]
                    results.loc[(contig, 'Nr genes with {}'.format(gene_label_to_gcode[1])),
                                "Segment {}".format(n)] = nr_confident_1 - nr_isoforms
                    results.loc[(contig, 'Nr genes with {}'.format(gene_label_to_gcode[2])),
                                "Segment {}".format(n)] = nr_confident_2 - nr_isoforms
                    results.loc[(contig, 'Nr identical genes'), "Segment {}".format(n)] = nr_identical - nr_isoforms
                    results.loc[(contig, 'Nr isoforms'), "Segment {}".format(n)] = nr_isoforms
                    results.loc[(contig, 'Total gene count'), "Segment {}".format(n)] = (nr_identical - nr_isoforms) +\
                                                                                        (nr_confident_1 - nr_isoforms) +\
                                                                                        (nr_confident_2 - nr_isoforms) +\
                                                                                        nr_isoforms
                results.loc[(contig, 'Nr genes with {}'.format(gene_label_to_gcode[1])),
                            "Sequence"] = results.loc[(contig, 'Nr genes with {}'.format(gene_label_to_gcode[1]))].sum()
                results.loc[(contig, 'Nr genes with {}'.format(gene_label_to_gcode[2])),
                            "Sequence"] = results.loc[(contig, 'Nr genes with {}'.format(gene_label_to_gcode[2]))].sum()
                results.loc[(contig, 'Nr identical genes'), "Sequence"] = results.loc[(contig, 'Nr identical genes')].sum()
                results.loc[(contig, 'Nr isoforms'), "Sequence"] = results.loc[(contig, 'Nr isoforms')].sum()
                results.loc[(contig, 'Total gene count'), "Sequence"] = results.loc[(contig, 'Total gene count')].sum()
        results.fillna('-', inplace=True)
        results.to_csv(self.path_to_output, sep='\t', index=True)

    def write_stats(self):
        """
        Write stats file: gene id, number recoded stops, positions of recoded stops + stop codon
        """
        with open(self.path_to_output.replace('.tab', '_stats.tab'), 'w+') as stat:
            for n, recodings in enumerate(self.stop_recodings):
                recoded_positions = []
                [recoded_positions.extend(x)for x in recodings.values()]
                recoded_stops = []
                [recoded_stops.extend([x] * len(recodings[x])) for x in recodings.keys()]
                if len(set(recoded_stops)) == 1:
                    pass
                else:
                    sorted_inds = np.argsort(recoded_positions)
                    recoded_positions = np.array(recoded_positions)[sorted_inds]
                    recoded_stops = np.array(recoded_stops)[sorted_inds]
                string_to_write = f"{self.results.coords_w_isoforms.loc[n, 'attributes'].split(';')[0]} recodings {len(recoded_positions)} codon_index_values"
                for pos, stop in zip(recoded_positions, recoded_stops):
                    string_to_write += f" {pos} {stop}"
                string_to_write = string_to_write.replace(' ', '\t')
                string_to_write += '\n'
                stat.writelines(string_to_write)

    def write_final_results(self):
        """
        Write gff file
        """
        mgm_version = str(subprocess.run([f"{self.cwd}/dependencies/gmhmmp2"], capture_output=True).stderr).split("\\n")[1]
        if isinstance(self.results.genome_label, int):
            translation_tables = [str(self.results.genome_label)]
        else:
            translation_tables = [str(self.results.label_to_gcode[x]) for y in self.results.genome_label for x in y]
        translation_tables = sorted(set(translation_tables))
        parameters = ", ".join([f"{self.cwd}/dependencies/mgm_{x}.mod" for x in translation_tables])
        translation_tables = ", ".join(translation_tables)
        size = len(self.results.sequence[0].seq)
        try:
            accession = self.results.coords_w_isoforms.iloc[0, 0]
        except IndexError:
            warnings.warn('No genes were predicted, automatically assign genetic code 11')
        # restore correct column order
        original_column_order= self.results.label_to_gene_coords[0].columns.values
        with open(self.path_to_final_annotations, 'w+') as output:
            #write gff header
            output.writelines('##gff-version 2\n')
            output.writelines('#' + mgm_version + "\n")
            output.writelines(f"# File with sequence: {self.path_to_genome}\n")
            output.writelines(f"# File with MetaGeneMark parameters: {parameters}\n")
            output.writelines(f"# translation table: {translation_tables}\n")
            output.writelines(f"# output date start: {self.date.strftime('%a %b %d %Y %H:%M:%S')}\n")
            for n in range(self.results.nr_contigs):
                size = len(self.results.sequence[n].seq)
                current_coords = self.results.coords_w_isoforms[self.results.coords_w_isoforms.contig ==
                                                                self.results.sequence[n].id].copy()
                current_coords = current_coords.loc[:, original_column_order]
                #current_coords = self.results.coords[self.results.coords.contig == self.results.sequence[n].id].copy()
                if current_coords.shape[0] == 0:
                    continue
                current_coords.reset_index(drop=True, inplace=True)
                accession = current_coords.iloc[0, 0]
                output.writelines("\n")
                output.writelines(f"##sequence-region {accession} 1 {size}\n")
                # filter out isoforms for statistics
                inds_of_isoform_parents = current_coords[current_coords.logodd == '.'].index.values.tolist()
                inds_of_isoforms = []
                [inds_of_isoforms.append(i + 1) for i in inds_of_isoform_parents]
                [inds_of_isoforms.append(i + 2) for i in inds_of_isoform_parents]
                inds_of_isoforms = sorted(inds_of_isoforms)
                current_coords_without_isoforms = current_coords[(~current_coords.index.isin(inds_of_isoforms)) &
                                                                 (~current_coords.index.isin(inds_of_isoform_parents))]
                mean_length_isoforms = []
                mean_logodd_isoforms = []
                for i in range(0, len(inds_of_isoforms), 2):
                    mean_length_isoforms.append(current_coords.loc[[inds_of_isoforms[i], inds_of_isoforms[i + 1]],
                                                                   'length'].mean())
                    mean_logodd_isoforms.append(current_coords.loc[[inds_of_isoforms[i], inds_of_isoforms[i + 1]],
                                                                   'logodd'].mean())
                # get average gene length
                average_gene_length = int(np.round(np.concatenate([current_coords_without_isoforms.length.values,
                                                                   mean_length_isoforms]).mean()))
                current_coords.drop(["length", "contig"], axis=1, inplace=True)
                # get total logodd score
                # TODO: how to correctly handle logodd scores of isoforms?
                total_logodd_score = np.round(current_coords_without_isoforms.logodd.sum() + sum(mean_logodd_isoforms),
                                              decimals=1)
                # get average gene density
                average_gene_density = np.round((current_coords.shape[0] -
                                                 current_coords[current_coords.logodd == '.'].shape[0] * 2) /
                                                (size / 1000), decimals=2)
                # write final predictions
                current_coords.to_csv(output, header=False, sep='\t', index=False)
                # write summary
                output.writelines(f"# {accession}\ttotal_logodd\t{total_logodd_score}\taverage_length\t"
                                  f"{average_gene_length}\taverage_density\t{average_gene_density}\n")

    @staticmethod
    def find_recoded_stops(seq, stop):
        total_counts = seq.count(stop)
        i = 0
        n = 0
        start_recodings = []
        while i <= len(seq) - 1 and n < total_counts:
            start = seq[i:].find(stop)
            # if in-frame must be multiple of 3 and not last codon
            if (i + start) % 3 == 0 and (start + i) != len(seq) - 3:
                seq = seq[: start + i] + seq[start + i: start + i + 3].upper() + seq[start + i + 3:]
                start_recodings.append(int((start + i) / 3))
            i += (start + 1)
            n += 1
        return start_recodings, seq

    @staticmethod
    def make_seqrecord(prediction, sequence, n):
        """
        Extract sequences of predicted genes and convert to Biopython SeqRecord obkect
        :param prediction: pandas Series, prediction
        :param sequence: Biopython Seq object, DNA sequence
        :param n: int, gene index
        :return: Biopython SeqRecord, SeqRecord of gene
                 dict, codon positions of recoded stops
        """
        if prediction.strand == '-':
            seq = sequence.seq[prediction.start -1: prediction.end].reverse_complement()
        elif prediction.strand == '+':
            seq = sequence.seq[prediction.start - 1: prediction.end]
        # find recoded in frame stop codons
        seq = seq.lower()
        stop_recodings = {}
        stop_recodings['TAG'], seq = Results.find_recoded_stops(seq, 'tag')
        stop_recodings['TGA'], seq = Results.find_recoded_stops(seq, 'tga')
        stop_recodings['TAA'], seq = Results.find_recoded_stops(seq, 'taa')
        record = SeqRecord(seq, id=f"{n + 1} {sequence.id}", description=f"{prediction.start}..{prediction.end},"
                                                                         f"{prediction.source}, {prediction.attributes} "
                                                                         f"number recodings: TAA {len(stop_recodings['TAA'])},"
                                                                         f"TGA {len(stop_recodings['TGA'])},"
                                                                         f"TAG {len(stop_recodings['TAG'])}")
        return record, stop_recodings

    @staticmethod
    def translate_seq(prediction, nucleotide_seq, start_recodings):
        """
        Translates nucleotide sequences of predicted genes.
        :param prediction: pandas Series, prediction
        :param nucleotide_seq: Biopython SeqRecord, nucleotide sequence of corresponding prediction
        :param start_recodings: dict, codon positions of recoded stops
        :return: BioPython SeqRecord, amino acid sequence of corresponding gene
        """
        #initialize trans table
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
        stop_codons = ['TAA', 'TAG', 'TGA']
        start_codons = ['TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG']
        # adjust translation table based on genetic code
        if '4' in prediction.source:
            current_stop_codons = stop_codons.copy()
            current_stop_codons.remove('TGA')
            current_table = table.copy()
            current_table['TGA'] = 'W'
        elif '15' in prediction.source:
            current_stop_codons = stop_codons.copy()
            current_stop_codons.remove('TAG')
            current_table = table.copy()
            current_table['TAG'] = 'Q'
        elif '101' in prediction.source:
            current_stop_codons = stop_codons.copy()
            current_stop_codons.remove('TAA')
            current_table = table.copy()
            current_table['TAA'] = 'Y'
        # genetic code 11
        else:
            current_stop_codons = stop_codons
            current_table = table
        codon_table = CodonTable(forward_table=current_table, start_codons=start_codons,
                                 stop_codons=current_stop_codons)
        # translate to amino acids
        aa_seq = nucleotide_seq.translate(table=codon_table, stop_symbol='')
        # aa from recoded stops will be capital case letter, everything else lower case letter
        aa_seq = aa_seq.lower()
        recoded_aas = []
        [recoded_aas.extend(x) for x in start_recodings.values()]
        for position in recoded_aas:
            aa_seq = aa_seq[:position] + aa_seq[position].upper() + aa_seq[position + 1:]
        aa_seq.id = nucleotide_seq.id
        aa_seq.description = nucleotide_seq.description
        return aa_seq

    def extract_sequences(self):
        """
        Extract nucleotide, amino acid sequences
        """
        nucleotide_seqs = []
        stop_recodings = []
        # get nucleotide seqs and positions of recoded stops
        for n in range(self.results.nr_contigs):
            current_coords = self.results.coords_w_isoforms[self.results.coords_w_isoforms.contig == self.results.sequence[n].id].copy()
            current_coords = current_coords[current_coords.logodd != '.']
            current_coords.reset_index(drop=True, inplace=True)
            for i in range(current_coords.shape[0]):
                nucleotide_seq, stop_recoding = self.make_seqrecord(current_coords.iloc[i], self.results.sequence[n], i)
                nucleotide_seqs.append(nucleotide_seq)
                stop_recodings.append(stop_recoding)
        # write nucleotide seqs
        if self.nucleotides:
            SeqIO.write(nucleotide_seqs, self.path_to_nt_seqs, "fasta")
        # translate nucleotide seqs and write aa seqs
        if self.amino_acids:
            amino_acid_seqs = []
            for i in range(current_coords.shape[0]):
                aa_seq = self.translate_seq(current_coords.iloc[i], nucleotide_seqs[i],
                                            stop_recodings[i])
                if not aa_seq is None:
                    amino_acid_seqs.append(aa_seq)
            SeqIO.write(amino_acid_seqs, self.path_to_aa_seqs, "fasta")
        return stop_recodings
