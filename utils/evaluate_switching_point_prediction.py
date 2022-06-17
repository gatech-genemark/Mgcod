#!/usr/bin/env python
"""
Script to evaluate the accuracy with which Mgcod predicts the switch of genetic code in simulated
genomes with multiple genetic codes
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import sys
import argparse


def get_predictions(input_dir):
    """
    Function to parse Mgcod output files (*.tab)
    :param input_dir: str, directory containing Mgcod output files
    :return: pd.DataFrame, dataframe of segmentations of target genomes
    """
    predictions = pd.DataFrame(columns=['genetic_code', 'coordinates'])
    for i in range(1000):
        prediction_path = f"{input_dir}switch_point_{i}.tab"
        try:
            c_prediction = pd.read_csv(prediction_path, header=0, index_col=0, sep='\t')
        except FileNotFoundError:
            continue
        if not 'Segment 1' in c_prediction.columns:
            coordinates = [(1, 200000)]
            genetic_codes = [c_prediction.iloc[0].values[1]]
        else:
            coordinates = [(int(x[1:-1].split(', ')[0]), int(x[1:-1].split(', ')[1]))
                           for x in c_prediction.iloc[0].values[1:-1]]
            genetic_codes = [x for x in c_prediction.iloc[1].values[1:-1]] 
        predictions.loc[f"switching_point_{i}"] = [genetic_codes, coordinates]
    return predictions


def filter_annotation_predictions(annotations, predictions):
    """
    Filter out Mgcod predictions that predicted the wrong sequence of genetic codes
    :param annotations: pd.DataFrame, annotated segmentation of genomes
    :param predictions: pd.DataFrame, predicted segmentation of genomes
    :return: pd.DataFrame, predicted segmentation of genomes with correct sequence of multiple genetic codes
    """
    annotations = annotations[[len(x) == 2 for x in predictions.genetic_code.values]]
    predictions = predictions[[len(x) == 2 for x in predictions.genetic_code.values]]
    to_keep = []
    for anno, pred in zip(annotations.switch.values, predictions.genetic_code.values):
        if anno == "[11, 4]" and (pred == ["11", "4"] or pred == ["11_4", "4"] or pred == ["11", "11_4"]):
            to_keep.append(True)
        elif anno == "[4, 11]" and (pred == ["4", "11"] or pred == ["11_4", "11"] or pred == ["4", "11_4"]):
            to_keep.append(True)
        else:
            to_keep.append(False)
    annotations = annotations[to_keep]
    predictions = predictions[to_keep]
    return annotations, predictions


def get_prediction_error(annotations, predictions, switch_type=1):
    """
    Calculate haw far the predicted genetic code switch is in bp
    :param annotations: pd.DataFrame, annotated genetic code switch coordinates
    :param predictions: pd.DataFrame, predicted genetic code switch coordinates
    :param switch_type: int, code indicating type of switch
    :return: np.array, error in predicted genetic code switch coordinates in bp for all target genomes
    """
    if switch_type == 1:
        predictions = predictions[(annotations.switch == "[11, 4]").values]
        annotations = annotations[(annotations.switch == "[11, 4]").values]
    elif switch_type == 2:
        predictions = predictions[(annotations.switch == "[4, 11]").values]
        annotations = annotations[(annotations.switch == "[4, 11]").values]
    annotated_switch_points = np.array([((int(x[1:-1].split(', ')[0]) + int(x[1:-1].split(', ')[1])) / 2)
                                        for x in annotations.switch_point.values])
    predicted_switch_points = np.array([int((x[0][1] + x[1][0]) / 2) for x in predictions.coordinates.values])
    return predicted_switch_points - annotated_switch_points


def plot_error(error_switch_type_1, error_switch_type_2, output_dir):
    """
    Plot histograms of errors, distinguishing the two types of error
    :param error_switch_type_1: np.array, errors in bp for type 1 switches
    :param error_switch_type_2: np.array, errors in bp for type 2 switches
    :param output_dir: str, directory where to save output figure error_boundary_prediction.pdf
    """
    fig, ax = plt.subplots()
    handle1 = Line2D([], [], c='blue')
    handle2 = Line2D([], [], c='green')
    bins = np.arange(-37500, 42500, 5000) 
    ax.hist(error_switch_type_1, alpha=0.7, bins=bins, align='right', weights=np.ones(error_switch_type_1.shape[0]) /
                                                                              error_switch_type_1.shape[0],
            color='blue')
    ax.hist(error_switch_type_2, alpha=0.7, bins=bins, align='left', weights=np.ones(error_switch_type_2.shape[0]) /
                                                                             error_switch_type_2.shape[0],
            color='green')
    ax.set_xlabel("Error in predicting block boundary in bp")
    ax.set_ylabel("Density")
    ax.legend(handles=[handle1, handle2], labels=[r'$11 \rightarrow 4$', r'$4 \rightarrow 11$',], ncol=2,
              loc='upper center', bbox_to_anchor=(0.5, -0.13))
    fig.savefig(f"{output_dir}error_boundary_prediction.pdf", bbox_inches='tight', dpi=600)


def main(argv):
    parser = argparse.ArgumentParser(description="Evaluate the accuracy with which Mgcod predicts the switch of "
                                                 "genetic code in simulated genomes with multiple genetic codes")
    parser.add_argument('-a', '--annotations', help='Path to switch point annotations of simulated genomes')
    parser.add_argument('-r', '--results', help='Directory to Mgcod results')
    parser.add_argument('-o', '--output_dir', help="Directory where to save figure")
    args = parser.parse_args()

    annotations = pd.read_csv(args.annotations, header=None, index_col=0, sep='\t', names=['switch', 'switch_point'])
    results_path = args.results
    if not results_path.endswith('/'):
        results_path += '/'
    output_dir = args.output_dir
    if not output_dir.endswith('/'):
        output_dir += '/'

    predictions = get_predictions(results_path)
    #annotations = annotations.loc[predictions.index.values] 
    switch_type_counts = annotations.switch.value_counts()
    print("Number of genomes with simulated switch:")
    print(switch_type_counts)
    correct_annotations, correct_predictions = filter_annotation_predictions(annotations, predictions)
    correct_switch_type_counts = correct_annotations.switch.value_counts()
    print("Number of genomes with correctly predicted switch:")
    print(correct_switch_type_counts)
    # 11 --> 4
    error_switch_type_1 = get_prediction_error(correct_annotations, correct_predictions, switch_type=1)
    # 4 --> 11
    error_switch_type_2 = get_prediction_error(correct_annotations, correct_predictions, switch_type=2) 
    print("Error switch type I: {} =/- {}".format(np.mean(error_switch_type_1), np.std(error_switch_type_1)))
    print("Error switch type II: {} =/- {}".format(np.mean(error_switch_type_2), np.std(error_switch_type_2)))
    print("Overall Error: {} +/- {}".format(np.mean(np.concatenate([error_switch_type_1, error_switch_type_2])),
                                            np.std(np.concatenate([error_switch_type_1, error_switch_type_2]))))
    print('Correct switches within 5000bp:')
    print('11 -> 4: {}'.format((np.abs(error_switch_type_1) <= 5000).sum()))
    print('4 -> 11: {}'.format((np.abs(error_switch_type_2) <= 5000).sum()))
    plot_error(error_switch_type_1, error_switch_type_2, output_dir)


if __name__ == '__main__':
    main(sys.argv[1:])
