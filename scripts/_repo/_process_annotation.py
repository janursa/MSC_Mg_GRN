"""
Map different annotations to each other.
String data should be available beforehand.
# geneID (entrenz gene id). https://www.uniprot.org/id-mapping/
"""
import sys
import os
import typing
import json
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import ENRICH_DIR, DATA_DIR


def main(DE_protnames: list, phase) -> typing.Tuple[dict,dict,dict,dict]:
    #- mao protname to geneID
    map_protname_geneID = annotaion_map.protname_geneID(DE_protnames)
    #- map protname to genename according to model_selection: because it is different using the website
    string_map = pd.read_csv(os.path.join(ENRICH_DIR, f'string_mapping_{phase}.tsv'), sep='\t',
                             index_col=False)
    map_protname_genename = {key: value for key, value in zip(string_map['queryItem'], string_map['preferredName'])}
    #- check if model_selection imported all protnames successfully
    try:
        assert (len(map_protname_genename) == len(DE_protnames))
    except:
        for prot in DE_protnames:
            if prot not in string_map['queryItem'].values:
                print(f'prot {prot} was not mapped in vs_string')
    #- gene name to protname
    map_genename_protname = {value: key for value, key in
                             zip(map_protname_genename.values(), map_protname_genename.keys())}
    #- genename to geneID
    map_genename_geneID = {key: map_protname_geneID[value] for key, value in map_genename_protname.items()}

    return map_protname_geneID, map_protname_genename, map_genename_protname, map_genename_geneID
if __name__ == '__main__':
    for phase in ['early','late']:
        #- retreive DE proteins
        with open(os.path.join(DATA_DIR, 'DE_protnames.txt')) as f:
            DE_protnames = eval(f.read())[phase]
        #- main function
        map_protname_geneID, map_protname_genename, map_genename_protname, map_genename_geneID = main(DE_protnames, phase)
        #- output
        # with open(os.path.join(OUTPUT_DIR,'data', 'map_protname_geneID.json'), 'w') as f:
        #     f.write(json.dumps({'map': map_protname_geneID}))

        # with open(os.path.join(OUTPUT_DIR, 'data', 'map_protname_genename.json'), 'w') as f:
        #     f.write(json.dumps({'map': map_protname_genename}))

        with open(os.path.join(ENRICH_DIR, f'map_genename_protname_{phase}.json'), 'w') as f:
            f.write(json.dumps({'map': map_genename_protname}))
        #
        # with open(os.path.join(OUTPUT_DIR, 'data', 'map_genename_geneID.json'), 'w') as f:
        #     f.write(json.dumps({'map': map_genename_geneID}))

        # write genenames to a file
        # DE_genenames = list(map_protname_genename.values())
        # DE_genenames_str = ""
        # for gene in DE_genenames:
        #     DE_genenames_str += gene + ","
        # with open(os.path.join(MAIN_DIR, 'results/postprocess', 'DE_genenames.json'), 'w') as f:
        #     json.dump({"DE_genenames": DE_genenames}, f)