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
    if 'early' in DE_type:
        legend_marker = False
        legend_size = False
        legend_color = False
        title = '(A) Short term'
        figsize = (3, 8)
        scale_factor = .4
    else:
        legend_marker = True
        legend_size = True
        legend_color = True
        title = '(B) Long term'
        figsize = (3.5, 8)
        scale_factor = 1.8
    fig = enrich_analysis.plot_enrich(df_targets, categories, size_tag, color_tag, xlabel, markers, figsize=figsize,
                                 legend_color=legend_color, legend_size=legend_size, legend_marker=legend_marker,
                                 title=title, scale_factor=scale_factor)
    ENRICH_PLOTS_DIR = Path(ENRICH_DIR)/'plots'
    if not os.path.isdir(ENRICH_PLOTS_DIR):
        os.makedirs(ENRICH_PLOTS_DIR)
    fig.savefig(ENRICH_PLOTS_DIR/ f'{DE_type}.png', dpi=300, transparent=True, bbox_inches='tight')
    fig.savefig(ENRICH_PLOTS_DIR/ f'{DE_type}.pdf', bbox_inches='tight')

def write_term_enrichment_data_to_file(df_stack: List[pd.DataFrame], categories: List[str], DE_type: str) -> None:
    """Writes enrichment data to files for each tag (for R plots)."""
    FOLDER = Path(ENRICH_DIR)/'For_R'
    if not os.path.isdir(FOLDER):
        os.makedirs(FOLDER)
    for term, df_term in zip(categories, df_stack):
        tag = term.split(" ")[-1].lower()
        file_name = Path(FOLDER) / f'enrichment_{DE_type}_{tag}.csv'
        with open(file_name, 'w') as f:
            df_term.to_csv(file_name, index=False)
def change_protnames_2_genenames(df):
    protnames_stack = df['ProteinNames']
    genenames_stack = []
    for protnames in protnames_stack:
        for protname, genename in F_protnames_to_genenames().items():
            protnames = protnames.replace(protname, genename)
        genenames_stack.append(protnames)
    df['ProteinNames'] = genenames_stack
    return df

if __name__ == '__main__':
    # - parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--top_n', type=int, default=10, help='Top number of enriched terms to show for each category')
    parser.add_argument('--length_limit', type=int, default=40, help='Limits term length to be fit into the graph')
    args, remaining_args = parser.parse_known_args()
    top_n = args.top_n
    length_limit = args.length_limit
    #- EA for each model
    for model_name in F_selected_models()[1:]: #TODO: fix
        from local_utils import * #TODO: togo
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
            df_target = change_protnames_2_genenames(df_target)
            df_targets.append(df_target)
        dff = pd.concat(df_targets).reset_index()
        dff.drop(columns=['index','ProteinNames', 'BackgroundGeneCount'], inplace=True)
        print_df(dff)
        aa
        #- output each term's df (for R plot)
        write_term_enrichment_data_to_file(df_targets, categories, DE_type)
        #- plot
        plot_enrichment_data(df_targets, categories, DE_type)

