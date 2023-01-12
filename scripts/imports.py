from utils import *
import utils
# import utils
OUTPUT_DIR = os.path.join(MAIN_DIR, 'results')
CALIBRATION_DIR = os.path.join(OUTPUT_DIR, 'calibration')

#- import training data
df_target = pd.read_csv(os.path.join(MAIN_DIR,'results/postprocess/DE_data.csv'))
#- DE protnames
protnames = list(df_target['Protein'].values)
n_features = len(protnames)
#- funcs
def retreive_grn(study, method):
    if method == 'portia' or method=='ridge':
        links = utils.links.read_write_links(study=study, mode='read_links',method=method,OUTPUT_DIR=OUTPUT_DIR)
    else:
        links = pd.read_pickle(os.path.join(OUTPUT_DIR,'GRN', f'links_{study}_{method}.csv'))
    return links

param_grid_RF = dict(
    decay_coeff=np.arange(0,1,.05),
    # min_samples_leaf=np.arange(1,5,1),
    # max_depth=np.arange(5,34,1),
    max_features=np.arange(int(np.sqrt(n_features)),n_features,1),
)

param_grid_ridge = get_estimator_wrapper('ridge').get_grid_parameters()
param_grid_ridge = {**param_grid_ridge,'decay_coeff':np.arange(0,1,.05)}

