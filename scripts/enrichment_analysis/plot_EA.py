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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.utils import enrichplot
from scripts.imports import ENRICH_DIR, F_DE_protiens


if __name__ == '__main__':
    selected_models = ['day1_11_KNN', 'day1_21_KNN']
    for DE_type, DE_proteins in F_DE_protiens().items():
        if DE_type not in selected_models:
            continue
        # - assignments
        size_tag, color_tag, xlabel = 'Strength', 'FDR', 'ProteinCount'
        markers = ['o', '>', '*', 's', 'p']
        #- retrieve the enriched terms and reformat them
        df = pd.read_csv(os.path.join(ENRICH_DIR, f'enrichment_all_{DE_type}.csv'), index_col=False)
        # - shorten term names for visualization purposes
        term_length_cut_threshold = 50
        short_descript = []
        for item in df.loc[:, 'Description']:
            if len(item) > term_length_cut_threshold:
                aa = item[0:term_length_cut_threshold] + '...'
            else:
                aa = item
            short_descript.append(aa)
        df['Description'] = short_descript
        #-----
        # terms = ['GO Biological Process', 'GO Biological Function', 'GO Cellular Component', 'Reactome Pathways', 'Uniprot Keywords']
        terms = ['GO Biological Process', 'GO Biological Function', 'GO Cellular Component']

        df_targets = []
        for term in terms:
            df_targets.append(df.loc[df['category']==term,:])
        #- keep only top n enriched items
        top_n = 10 #based on Strength
        filter_tag = 'Strength'
        df_targets = [df.sort_values(filter_tag, ascending=False) for df in df_targets]
        df_targets = [df.head(top_n) for df in df_targets]
        #- plot
        if 'day1_11' in DE_type:
            legend_marker=False
            legend_size = False
            legend_color = False
            title = '(A) Early phase'
            figsize = (3, 8)
            scale_factor = .4
        else:
            # stren = df['Strength'].to_numpy(str)
            # print(min(stren),max(stren))
            # aa
            legend_marker = True
            legend_size = True
            legend_color = True
            title = '(B) Late phase'
            figsize = (3.5, 8)
            scale_factor = 1.8
        fig = enrichplot.plot_enrich(df_targets, terms, size_tag, color_tag, xlabel, markers, figsize=figsize,
                                   legend_color=legend_color, legend_size=legend_size, legend_marker=legend_marker,
                                   title=title, scale_factor=scale_factor)
        fig.savefig(os.path.join(ENRICH_DIR, f'{title}.png'), dpi=300, transparent=True, bbox_inches='tight')
        fig.savefig(os.path.join(ENRICH_DIR, f'{title}.pdf'), bbox_inches='tight')
    # plt.show()
