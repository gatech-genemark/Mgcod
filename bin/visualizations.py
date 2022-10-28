#!/usr/bin/env python
# -------------------------------------------------
# Mgcod
# Aaron Pfennig
# Georgia Institute of Technology
# -------------------------------------------------

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import os


class Visualizations:
    def __init__(self, pars, verbosity):
        """
        Plots logodd ratios for each model. Each model is represented in one panel. The middle of each window is
        indicated by a dot on the bar. The highest scoring model is indicated by a dot at the top of the corresponding
        panel
        :param pars: dict with parameters
        :param verbosity: boolean, verbosity
        """
        self.verbose = verbosity
        self.logodd_scores = np.repeat(pars['logodd_scores'], repeats=2, axis=0)[:, pars['labels']]
        self.length_sequence = pars['length_sequence']
        self.sequence = pars['sequence']
        self.window_size = pars['window_size']
        self.stride = pars['stride']
        self.plot_dir = pars['plot_dir']
        self.species = pars['species']
        # self.strands = pars['strand']
        # self.starts = pars['gene_starts']
        self.switch_regions = pars['switch_regions']
        self.labels = pars['labels']
        # self.confident_1 = pars['confident_1']
        # self.confident_2 = pars['confident_2']
        # `self.ambiguous = pars['ambiguous']
        self.compared_gene_set = pars['compared_gene_set'].reset_index(drop=True)
        self.gene_label_to_gcode = pars['gene_label_to_gcode']
        # f len(self.gene_label_to_gcode) == 1:
        #   self.gene_label_to_gcode = self.gene_label_to_gcode[0]
        self.label_to_gcode = {0: 11, 1: 4, 2: 15, 3: 101}
        self.gcode_to_label = {11: 0, 4: 1, 15: 2, 101: 3}

    def __call__(self):
        self.logodd_ratios = self.get_logodd_ratios()
        self.plot_logodd_ratios()

    def get_logodd_ratios(self):
        """
        Compute logodd ratios
        :return: numpy array [number mgm models, number windows], with log odd ratios
        """
        np.seterr(all='ignore')
        try:
            logodd_ratios = np.array(self.logodd_scores) / np.array(self.logodd_scores).sum(axis=1)[:, np.newaxis]
        except np.AxisError:
            logodd_ratios = np.array(self.logodd_scores) / np.array(self.logodd_scores)
        return np.nan_to_num(logodd_ratios)

    def plot_logodd_ratios(self):
        """
        Plot logodd ratios in panels
        """
        legend_elements = {}
        colors = ['blue', 'green', 'red', 'orange']
        x_values_1 = np.arange(0, len(self.sequence), self.window_size)
        x_values_2 = np.arange(self.window_size / 2, len(self.sequence) + self.window_size / 2, self.window_size)
        # center of windows
        annotations_x = np.arange(self.window_size / 2, len(self.sequence) + self.window_size / 2, self.window_size)
        x_values = [None] * (len(x_values_2) + len(x_values_2))
        annotations_y = self.logodd_ratios[::2]
        # highest scoring model in window
        try:
            highest_scoring_models = np.argmax(self.logodd_ratios[::2], axis=1)
        except np.AxisError:
            highest_scoring_models = np.argmax(self.logodd_ratios[::2])
        highest_scoring_models_matrix = np.zeros_like(self.logodd_ratios[::2])
        for i in range(self.logodd_ratios.shape[1]):
            winning_windows = np.unique(np.where(highest_scoring_models == i)[0])
            highest_scoring_models_matrix[winning_windows, i] = annotations_x[winning_windows]
        x_values[::2] = x_values_1
        x_values[1::2] = x_values_2
        fig = plt.figure()
        gs = fig.add_gridspec(14, 14, hspace=0)

        # plot super ylabel
        ax2 = fig.add_subplot(gs[7:9, :1])
        ax2.set_ylabel('Normalized log-odds scores', labelpad=19, fontdict={'fontsize': 10})
        ax2.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False,
                        labelright=True)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax2.spines['left'].set_visible(False)

        ax = list()
        ax.append(fig.add_subplot(gs[2:8, 1:]))
        ax.append(fig.add_subplot(gs[8:, 1:], sharex=ax[0]))

        # plot encoding strand
        ax1 = fig.add_subplot(gs[0:2, 1:], sharex=ax[0])

        starts_pos = self.compared_gene_set.start.values[np.where(self.compared_gene_set.strand.values == '+')[0]]
        starts_neg = self.compared_gene_set.start.values[np.where(self.compared_gene_set.strand.values == '-')[0]]
        stops_pos = self.compared_gene_set.end.values[np.where(self.compared_gene_set.strand.values == '+')[0]]
        stops_neg = self.compared_gene_set.end.values[np.where(self.compared_gene_set.strand.values == '-')[0]]
        for start, stop in zip(starts_pos, stops_pos):
            ax1.axvspan(start, stop, 0.75, 1, color='black')
        for start, stop in zip(starts_neg, stops_neg):
            ax1.axvspan(start, stop, 0.25, 0.5, color='black')
        ax1.set_yticks([5, 15])
        ax1.set_yticklabels(['Genes predicted in negative strand', 'Genes predicted in positive strand'], fontsize=10)
        ax1.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False,
                        labelright=True)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        # underlay genes predicted by specific models in corresponding color
        for i in range(self.compared_gene_set.shape[0]):
            gene_label = int(self.compared_gene_set.loc[i, 'gene_label'])
            gcode = self.gene_label_to_gcode[gene_label]
            if gene_label == 1 or gene_label == 2:
                gcode_label = self.gcode_to_label[int(gcode)]
                ax_ind = self.labels.index(self.gcode_to_label[int(gcode)])
                ax[ax_ind].axvspan(self.compared_gene_set.loc[i, 'start'],
                                        self.compared_gene_set.loc[i, 'end'], facecolor=colors[gcode_label],
                                        alpha=0.6, zorder=0)
            # underlay dual coded segment in gray
            elif gene_label == 3:
                ax[0].axvspan(self.compared_gene_set.loc[i, 'start'],
                              self.compared_gene_set.loc[i, 'end'], facecolor='lightgrey', zorder=0, alpha=0.5)
                ax[1].axvspan(self.compared_gene_set.loc[i, 'start'],
                              self.compared_gene_set.loc[i, 'end'], facecolor='lightgrey', zorder=0, alpha=0.5)

        for n in range(self.logodd_ratios.shape[1]):
            ax[n].plot(x_values, self.logodd_ratios[:, n], lw=2, ls='-', c=colors[self.labels[n]], zorder=0)
            # mark center of window
            ax[n].scatter(annotations_x, annotations_y[:, n], marker='o', s=6, color=colors[self.labels[n]],)
            # mark highest scoring model
            x_coords = np.unique(highest_scoring_models_matrix[:, n])[1:]
            if x_coords.shape[0] >= 1:
                y_coords = np.repeat(1.0, repeats=x_coords.shape[0])
                ax[n].scatter(x_coords, y_coords, marker='o', color='black', s=6,)

        for gene_label in self.compared_gene_set.gene_label.unique():
            gcode = self.gene_label_to_gcode[gene_label]
            if gene_label == 1 or gene_label == 2:
                gcode_label = self.gcode_to_label[int(gcode)]
                if not str(gcode) in legend_elements:
                    legend_elements[str(gcode)] = Patch(facecolor=colors[gcode_label], alpha=0.6,
                                                        label="Genes predicted in code " + str(gcode))
            elif gene_label == 3:
                if not str(gcode.replace('_', '&')) in legend_elements:
                    legend_elements[str(gcode.replace('_', '&'))] = Patch(facecolor='lightgrey', alpha=0.5,
                                                                          label="Genes predicted in codes " +\
                                                                                str(gcode.replace('_', '&')))
        legend_elements['SP'] = Line2D([0], [0], color='black', ls='--', label='Positions of predicted code switches')

        for sp in self.switch_regions:
            if sp[1] < sp[0]:
                continue
            ax[0].axes.axvline((sp[0] + sp[1]) / 2, color='black', ls='--', alpha=0.75)
            ax[1].axes.axvline((sp[0] + sp[1]) / 2, color='black', ls='--', alpha=0.75)
        legend_elements = [v for v in legend_elements.values()]

        # set limits
        [ax[i].set_xlim([0.0, self.length_sequence]) for i in range(self.logodd_ratios.shape[1])]
        [ax[i].set_ylim([-0.05, 1.05]) for i in range(self.logodd_ratios.shape[1])]
        # set labels
        [ax[i].set_ylabel('Code {}'.format(self.label_to_gcode[self.labels[i]], fontdict={'fontsize': 8})) for i in
         range(self.logodd_ratios.shape[1])]
        [ax[i].get_xaxis().set_visible(False) for i in range(self.logodd_ratios.shape[1])]
        ax[-1].set_xlabel('Position in genome in bp', fontdict={'fontsize': 10})
        ax[1].legend(handles=legend_elements, bbox_to_anchor=(1.01, 1.3),
                     frameon=False, fontsize=10)  # bbox_to_anchor=(0.5, -.14), loc='upper center', ncol=3)

        ax[-1].get_xaxis().set_visible(True)
        ax[-1].tick_params(axis='x', labelrotation=45, labelsize=9)
        ax[-1].tick_params(axis='y', labelsize=9)
        ax[0].tick_params(axis='y', labelsize=9)

        # check if path exists otherwise create
        if not os.path.exists(self.plot_dir):
            os.mkdir(self.plot_dir)
        # save
        fig.savefig(
            self.plot_dir + 'logodds_scores_in_{}_{}.pdf'.format(self.species, self.sequence.id), bbox_inches='tight',
            dpi=600)
        plt.close('all')

