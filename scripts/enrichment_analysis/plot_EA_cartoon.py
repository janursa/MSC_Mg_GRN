"""
Map different annotations to each other.
String data should be available beforehand.
# geneID (entrenz gene id). https://www.uniprot.org/id-mapping/
"""
import sys
import os
import typing
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from common_tools import enrich_analysis
from imports import ENRICH_DIR, F_selected_models, F_model_name_2_method_and_DE_type, F_protnames_to_genenames

def load_enrichment_data(de_type: str) -> pd.DataFrame:
    """Load enrichment data for the given DE type from a CSV file."""
    filename = os.path.join(ENRICH_DIR, f'enrichment_all_{de_type}.csv')
    with open(filename, 'r') as f:
        df = pd.read_csv(f, index_col=False)
    return df


def shorten_term_names(df: pd.DataFrame, threshold: int = 50) -> pd.DataFrame:
    """Shorten the description field of a DataFrame if it exceeds a given threshold."""
    short_descript = []
    for item in df.loc[:, 'Description']:
        if len(item) > threshold:
            desc = item[0:threshold] + '...'
        else:
            desc = item
        short_descript.append(desc)
    df['Description'] = short_descript
    return df


def select_top_enriched_terms(df: pd.DataFrame, category: str, top_n: int = 10, length_limit:int = 30, filter_tag: str = 'Strength') -> pd.DataFrame:
    """Select the top n enriched terms for the given category, sorted by the given filter tag."""
    df = df.loc[df['category'] == category, :]
    descriptions = df['Description']
    descriptions = [name[0:length_limit] for name in descriptions]
    df['Description'] = descriptions
    df = df.sort_values(filter_tag, ascending=False).head(top_n)
    return df


def plot_enrichment_data(df_targets: list, categories: list, DE_type: str):
    """Plot the enrichment data for the given terms using the specified parameters."""
    # - assignments
    size_tag, color_tag, xlabel = 'Strength', 'FDR', 'ProteinCount'
    markers = ['o', '>', '*', 's', 'p']
    legend_marker = False
    legend_size = False
    legend_color = False
    figsize = (3, 4)
    scale_factor = 1.3

    fig = enrich_analysis.plot_enrich(df_targets, categories, size_tag, color_tag, xlabel, markers, figsize=figsize,
                                 legend_color=legend_color, legend_size=legend_size, legend_marker=legend_marker,
                                 title='', scale_factor=scale_factor)
    ENRICH_PLOTS_DIR = Path(ENRICH_DIR)/'plots'
    if not os.path.isdir(ENRICH_PLOTS_DIR):
        os.makedirs(ENRICH_PLOTS_DIR)
    fig.savefig(ENRICH_PLOTS_DIR/ f'{DE_type}_cartoon.png', dpi=300, transparent=True, bbox_inches='tight')
    fig.savefig(ENRICH_PLOTS_DIR/ f'{DE_type}_cartoon.pdf', bbox_inches='tight')



if __name__ == '__main__':
    # - parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--top_n', type=int, default=4, help='Top number of enriched terms to show for each category')
    parser.add_argument('--length_limit', type=int, default=20, help='Limits term length to be fit into the graph')
    args, remaining_args = parser.parse_known_args()
    top_n = args.top_n
    length_limit = args.length_limit
    #- EA for each model
    model_name = 'early_MinProb'
    _, DE_type = F_model_name_2_method_and_DE_type(model_name)
    #- retrieve the enriched terms and reformat them
    df = load_enrichment_data(DE_type)
    # - shorten term names for visualization purposes
    df = shorten_term_names(df)
    # - divide df based on each term and keep only top n enriched items
    categories = ['GO Biological Process', 'GO Biological Function', 'GO Cellular Component']
    df_targets = []
    for category in categories:
        df_target = select_top_enriched_terms(df, category, top_n=top_n, length_limit=length_limit, filter_tag= 'Strength')
        df_targets.append(df_target)

    #- plot
    plot_enrichment_data(df_targets, categories, DE_type)

