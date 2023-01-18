"""
Map different annotations to each other.
String data should be available beforehand.
# geneID (entrenz gene id). https://www.uniprot.org/id-mapping/
"""
import sys
import os
import typing

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import *


if __name__ == '__main__':
    #- retrieve the enriched terms
    read_EA = lambda term: pd.read_csv(os.path.join(OUTPUT_DIR, 'enrichment_analysis', f'enrichment_{term}.csv'), index_col=False)
    data_process = read_EA('process')
    data_component = read_EA('component')
    data_function = read_EA('function')
    data_uniprot = read_EA('uniprot')
    data_reactome = read_EA('reactome')
    # - assignments
    size_tag = 'Strength'
    color_tag = 'FDR'
    xlabel = 'ProteinCount'
    # - data and their names
    datas_1 = [data_process]
    tags_1 = ['GO: Process']  # also shape tag
    marker_types_1 = ['o']
    datas_2 = [data_component]
    tags_2 = ['GO: Component']  # also shape tag
    marker_types_2 = ['o']
    datas_3 = [data_function, data_uniprot, data_reactome]
    tags_3 = ['GO: Biological Function', 'UniProt Keyword', 'Reactome Pathways']  # also shape tag
    marker_types_3 = ['*', 's', 'p']

    # fig = utils.enrichplot.plot(datas_1, tags_1, size_tag, color_tag, xlabel, marker_types_1, figsize=(6.2, 5),
    #                            legend_color=True, legend_size=True, legend_marker=False, title='GO: Biological Process')
    # fig.savefig(os.path.join(OUTPUT_DIR, 'enrichment_analysis', 'biological_process.png'), dpi=300, transparent=True)

    fig = utils.enrichplot.plot(datas_2, tags_2, size_tag, color_tag, xlabel, marker_types_2, figsize=(6, 5),
                               legend_color=True, legend_size=True, legend_marker=False,
                               title='GO: Biological Component')
    fig.savefig(os.path.join(OUTPUT_DIR, 'enrichment_analysis', 'biological_component.png'), dpi=300, transparent=True)
    fig.savefig(os.path.join(OUTPUT_DIR, 'enrichment_analysis', 'biological_component.pdf'))


    fig = utils.enrichplot.plot(datas_3, tags_3, size_tag, color_tag, xlabel, marker_types_3, figsize=(6, 3.5),
                               legend_color=True, legend_size=True, legend_marker=True, title='')
    fig.savefig(os.path.join(OUTPUT_DIR, 'enrichment_analysis', 'multiple.png'), dpi=300, transparent=True,
                bbox_inches='tight')
    fig.savefig(os.path.join(OUTPUT_DIR, 'enrichment_analysis', 'multiple.pdf'),
                bbox_inches='tight')