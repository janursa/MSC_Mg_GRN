"""
Map different annotations to each other.
String data should be available beforehand.
# geneID (entrenz gene id). https://www.uniprot.org/id-mapping/
"""
import sys
import os
import typing
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.utils import enrichplot
from scripts.imports import *


if __name__ == '__main__':
    #- retrieve the enriched terms
    read_EA = lambda term: pd.read_csv(os.path.join(OUTPUT_DIR, 'enrichment_analysis/For_R', f'enrichment_{term}.csv'), index_col=False)
    data_function = read_EA('function')
    data_function = data_function.iloc[0:5, :]
    data_uniprot = read_EA('uniprot')
    data_uniprot = data_uniprot.iloc[0:5, :]
    data_reactome = read_EA('reactome')
    data_reactome = data_reactome.iloc[0:3,:]
    # - assignments
    size_tag = 'Strength'
    color_tag = 'FDR'
    xlabel = 'ProteinCount'
    # - data and their names
    marker_types_2 = ['o']
    datas_3 = [data_function, data_uniprot, data_reactome]
    tags_3 = ['Gene Ontology', 'UniProt Keyword', 'Reactome Pathways']  # also shape tag
    marker_types_3 = ['*', 's', 'p']


    fig = enrichplot.plot_enrich(datas_3, tags_3, size_tag, color_tag, xlabel, marker_types_3, figsize=(3.5, 3),
                               legend_color=True, legend_size=True, legend_marker=True,
                               title='')
    fig.savefig(os.path.join(OUTPUT_DIR, 'enrichment_analysis', 'cartoon.png'), bbox_inches='tight', dpi=300, transparent=True)
    fig.savefig(os.path.join(OUTPUT_DIR, 'enrichment_analysis', 'cartoon.pdf'), bbox_inches='tight')
    # plt.show()

