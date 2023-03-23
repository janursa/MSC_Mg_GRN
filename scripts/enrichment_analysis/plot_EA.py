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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from utils import enrich_analysis
from imports import ENRICH_DIR, F_DE_proteins

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


def select_top_enriched_terms(df: pd.DataFrame, term: str, top_n: int = 10, filter_tag: str = 'Strength') -> pd.DataFrame:
    """Select the top n enriched terms for the given category, sorted by the given filter tag."""
    df = df.loc[df['category'] == term, :]
    df = df.sort_values(filter_tag, ascending=False).head(top_n)
    return df


def plot_enrichment_data(df_targets: list, terms: list):
    """Plot the enrichment data for the given terms using the specified parameters."""
    # - assignments
    size_tag, color_tag, xlabel = 'Strength', 'FDR', 'ProteinCount'
    markers = ['o', '>', '*', 's', 'p']
    if 'day1_11' in DE_type:
        legend_marker = False
        legend_size = False
        legend_color = False
        title = '(A) Early phase'
        figsize = (3, 8)
        scale_factor = .4
    else:
        legend_marker = True
        legend_size = True
        legend_color = True
        title = '(B) Late phase'
        figsize = (3.5, 8)
        scale_factor = 1.8
    fig = enrich_analysis.plot_enrich(df_targets, terms, size_tag, color_tag, xlabel, markers, figsize=figsize,
                                 legend_color=legend_color, legend_size=legend_size, legend_marker=legend_marker,
                                 title=title, scale_factor=scale_factor)
    fig.savefig(os.path.join(ENRICH_DIR, f'{title}.png'), dpi=300, transparent=True, bbox_inches='tight')
    fig.savefig(os.path.join(ENRICH_DIR, f'{title}.pdf'), bbox_inches='tight')

def write_term_enrichment_data_to_file(df_stack: List[pd.DataFrame], terms: List[str], DE_type: str) -> None:
    """Writes enrichment data to files for each tag (for R plots)."""
    FOLDER = Path(ENRICH_DIR)/'For_R'
    if not os.path.isdir(FOLDER):
        os.makedirs(FOLDER)
    for term, df_term in zip(terms, df_stack):
        tag = term.split(" ")[-1].lower()
        file_name = Path(FOLDER) / f'enrichment_{DE_type}_{tag}.csv'
        with open(file_name, 'w') as f:
            df_term.to_csv(file_name, index=False)

if __name__ == '__main__':
    selected_models = ['day1_11_KNN', 'day1_21_KNN']
    for DE_type, DE_proteins in F_DE_proteins().items():
        if DE_type not in selected_models:
            continue
        #- retrieve the enriched terms and reformat them
        df = load_enrichment_data(DE_type)
        # - shorten term names for visualization purposes
        df = shorten_term_names(df)
        # - divide df based on each term and keep only top n enriched items
        terms = ['GO Biological Process', 'GO Biological Function', 'GO Cellular Component']
        df_targets = []
        for term in terms:
            df_targets.append(select_top_enriched_terms(df, term, top_n= 10, filter_tag= 'Strength'))
        #- output each term's df (for R plot)
        write_term_enrichment_data_to_file(df_targets, terms, DE_type)
        #- plot
        plot_enrichment_data(df_targets, terms)

