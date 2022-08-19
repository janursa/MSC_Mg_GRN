"""This is the example module.

This module does stuff.
"""
__all__ = ['a', 'b', 'c']
__version__ = '0.1'
__author__ = 'Jalil Nourisa'

import time
import numpy as np
import itertools
import operator 
import matplotlib.pyplot as plt
from pathos.pools import ParallelPool as Pool

from sklearn import tree
from sklearn import ensemble
from sklearn import base
from sklearn import utils
from sklearn import model_selection


#TODO: lag h: can be more than 1. can be in the format of hour/day 
#TODO: add stready state data
#TODO: https://peps.python.org/pep-0008/
#TODO: https://peps.python.org/pep-0257/
#TODO: https://github.com/numpy/numpydoc
#TODO: https://docs.python.org/3/library/doctest.html
#TODO: standalone examples to illustrate the usage, model visualisation, and benefits/benchmarks of particular algorithms;
#TODO: https://sklearn-template.readthedocs.io/en/latest/quick_start.html#develop-your-own-scikit-learn-estimators
#TODO: scripts to manage continuous integration (testing on Linux and Windows)

# TODO: sensitivity of the model to the parameters. 
#       This can be done both using sensitivity analysis or by evaluating the results of grid search 
#       (by checking the variations in the best chosen parameters)

# TODOC: pathos is used instead of multiprocessing, which is an external dependency. 
#        This is because multiprocessing uses pickle that has problem with lambda function.

# def SA(Xs, ys, param, grid_range, specs):
#     run_fun = # create the model by getting a grid_range, 
#     def run_fun()
#     specs['']
#     obj = barneySA.tools.SA(free_params = grid_range,settings = specs)
#     obj.sample()

#     obj.run()

#     obj.postprocessing()

def get_link_list(links, gene_names=None, regulators='all', maxcount='all', file_name=None):
    
    """Gets the ranked list of (directed) regulatory links.
    
    Parameters
    ----------
    
    
    gene_names: list of strings, optional
        List of length p, where p is the number of rows/columns in VIM, containing the names of the genes. The i-th item of gene_names must correspond to the i-th row/column of VIM. When the gene names are not provided, the i-th gene is named Gi.
        default: None
        
    regulators: list of strings, optional
        List containing the names of the candidate regulators. When a list of regulators is provided, the names of all the genes must be provided (in gene_names), and the returned list contains only edges directed from the candidate regulators. When regulators is set to 'all', any gene can be a candidate regulator.
        default: 'all'
        
    maxcount: 'all' or positive integer, optional
        Writes only the first maxcount regulatory links of the ranked list. When maxcount is set to 'all', all the regulatory links are written.
        default: 'all'
        
    file_name: string, optional
        Writes the ranked list of regulatory links in the file file_name.
        default: None
        
        
    
    Returns
    -------
    
    The list of regulatory links, ordered according to the edge score. Auto-regulations do not appear in the list. Regulatory links with a score equal to zero are randomly permuted. In the ranked list of edges, each line has format:
        
        regulator   target gene     score of edge
    """
    
    # Check input arguments     
    VIM =  np.array(links)
    # if not isinstance(VIM,ndarray):
    #     raise ValueError('VIM must be a square array')
    # elif VIM.shape[0] != VIM.shape[1]:
    #     raise ValueError('VIM must be a square array')
        
    ngenes = VIM.shape[0]
        
    if gene_names is not None:
        if not isinstance(gene_names,(list,tuple)):
            raise ValueError('input argument gene_names must be a list of gene names')
        elif len(gene_names) != ngenes:
            raise ValueError('input argument gene_names must be a list of length p, where p is the number of columns/genes in the expression data')
        
    if regulators != 'all':
        if not isinstance(regulators,(list,tuple)):
            raise ValueError('input argument regulators must be a list of gene names')

        if gene_names is None:
            raise ValueError('the gene names must be specified (in input argument gene_names)')
        else:
            sIntersection = set(gene_names).intersection(set(regulators))
            if not sIntersection:
                raise ValueError('The genes must contain at least one candidate regulator')
        
    if maxcount != 'all' and not isinstance(maxcount,int):
        raise ValueError('input argument maxcount must be "all" or a positive integer')
        
    if file_name is not None and not isinstance(file_name,str):
        raise ValueError('input argument file_name must be a string')
    
    

    # Get the indices of the candidate regulators
    if regulators == 'all':
        input_idx = list(range(ngenes))
    else:
        input_idx = [i for i, gene in enumerate(gene_names) if gene in regulators]
               
    nTFs = len(input_idx)
    
    # Get the non-ranked list of regulatory links
    vInter = [(i,j,score) for (i,j),score in ndenumerate(VIM) if i in input_idx and i!=j]
    
    # Rank the list according to the weights of the edges        
    vInter_sort = sorted(vInter,key=operator.itemgetter(2),reverse=True)
    nInter = len(vInter_sort)
    
    # Random permutation of edges with score equal to 0
    flag = 1
    i = 0
    while flag and i < nInter:
        (TF_idx,target_idx,score) = vInter_sort[i]
        if score == 0:
            flag = 0
        else:
            i += 1
            
    if not flag:
        items_perm = vInter_sort[i:]
        items_perm = np.random.permutation(items_perm)
        vInter_sort[i:] = items_perm
        
    # Write the ranked list of edges
    nToWrite = nInter
    if isinstance(maxcount,int) and maxcount >= 0 and maxcount < nInter:
        nToWrite = maxcount
        
    if file_name:
    
        outfile = open(file_name,'w')
    
        if gene_names is not None:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                outfile.write('%s\t%s\t%.6f\n' % (gene_names[TF_idx],gene_names[target_idx],score))
        else:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                outfile.write('G%d\tG%d\t%.6f\n' % (TF_idx+1,target_idx+1,score))
            
        
        outfile.close()
        
    else:
        
        if gene_names is not None:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                print('%s\t%s\t%.6f' % (gene_names[TF_idx],gene_names[target_idx],score))
        else:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                print('G%d\tG%d\t%.6f' % (TF_idx+1,target_idx+1,score))
def compute_feature_importances(estimator):
    
    """Computes variable importances from a trained tree-based model."""
    
    if isinstance(estimator, tree.BaseDecisionTree):
        return estimator.tree_.compute_feature_importances(normalize=False)
    else:
        importances = [e.tree_.compute_feature_importances(normalize=False)
                       for e in estimator.estimators_]
        importances = np.array(importances)
        return np.sum(importances,axis=0) / len(estimator)

def post_grid_search(best_params,best_scores):
    """plots the results of grid search"""
    sorted_best_params = {key: np.array([item[key] for item in best_params]) for key in best_params[0].keys()}
    n_sorted_best_params = {key: values/max(values) for key,values in sorted_best_params.items()}
 
    fig, axs = plt.subplots(1,2, tight_layout = True, figsize = (10,5),  gridspec_kw={'width_ratios': [3, 1]})
    axs[0].boxplot(n_sorted_best_params.values(), showfliers= True, labels=n_sorted_best_params.keys())
    axs[0].set_ylabel('Normalized quantity')
    axs[0].set_title('Estimated params stats')

    axs[1].boxplot(best_scores, showfliers= True)
    axs[1].set_ylabel('Score')
    axs[1].set_title('Best scores distribution')
    axs[1].set_xticklabels([])
def process_data(TS_data, time_points, regulators = 'all'): 
    """ Reformats the raw data for sklearn applications

    Receives data in the raw format (n_exp*n_time*n_genes) and returns a standard format for sklearn (n_samples*n_genes)
    for each gene.

    
    Arguments:
    TS_data --

    Return: 
    Xs -- A list of training inputs for each target gene, containing data in (n_exp*n_time - n_exp)*n_regulators format.
    y -- A list of training outputs for each target gene, containing data in (n_exp*n_time - n_exp) format
    """
    # Check input arguments
    if not isinstance(TS_data,(list,tuple)):
        raise ValueError('TS_data must be a list of arrays, where each row of an array corresponds to a time point/sample and each column corresponds to a gene')
    
    ngenes = TS_data[0].shape[1]

    if len(TS_data) > 1:
        for expr_data in TS_data[1:]:
            if expr_data.shape[1] != ngenes:
                raise ValueError('The number of columns/genes must be the same in every array of TS_data.')
    if not isinstance(time_points,(list,tuple)):
        raise ValueError('time_points must be a list of n one-dimensional arrays, where n is the number of time series experiments in TS_data')
    
    if len(time_points) != len(TS_data):
        raise ValueError('time_points must be a list of n one-dimensional arrays, where n is the number of time series experiments in TS_data')
    
    # for tp in time_points:
    #     if (not isinstance(tp,(list,tuple,ndarray))) or (isinstance(tp,ndarray) and tp.ndim > 1):
    #         raise ValueError('time_points must be a list of n one-dimensional arrays, where n is the number of time series in TS_data')
        
    for (i,expr_data) in enumerate(TS_data):
        if len(time_points[i]) != expr_data.shape[0]:
            raise ValueError('The length of the i-th vector of time_points must be equal to the number of rows in the i-th array of TS_data')

    # Re-order time points in increasing order
    for (i,tp) in enumerate(time_points):
        tp = np.array(tp, np.float32)
        indices = np.argsort(tp)
        time_points[i] = tp[indices]
        expr_data = TS_data[i]
        TS_data[i] = expr_data[indices,:]
    # obtain X and y for each target in a n_sample * n_feature, where n_sample is (n_exp*n_time - n_exp)
    Xs = [] 
    ys = []
    h = 1 # lag used for the finite approximation of the derivative of the target gene expression

    for i_gene in range(TS_data[0].shape[1]):
        if regulators == 'all':
            input_idx = list(range(TS_data[0].shape[1])) #TODO: needs to be based on the regulators
        else:
            input_idx = regulators[i_gene]
        ngenes = TS_data[0].shape[1]
        nexp = len(TS_data)
        nsamples_time = sum([expr_data.shape[0] for expr_data in TS_data]) 
        ninputs = len(input_idx)

        # Time-series data
        input_matrix_time = np.zeros((nsamples_time-h*nexp,ninputs))
        output_vect_time = []

        for (i,exp_timeseries) in enumerate(TS_data):
            exp_time_points = time_points[i]
            n_time = exp_timeseries.shape[0]
            exp_time_diff = exp_time_points[h:] - exp_time_points[:n_time-h]
            exp_timeseries_x = exp_timeseries[:n_time-h,input_idx]
            # current_timeseries_output = (exp_timeseries[h:,i_gene] - exp_timeseries[:n_time-h,i_gene]) / exp_time_diff + alphas[i_gene]*exp_timeseries[:n_time-h,i_gene]
            for ii in range(len(exp_time_diff)):
                f_dy_dt = lambda alpha_i,i=i, ii=ii,i_gene=i_gene: float((TS_data[i][ii+1:ii+2,i_gene] - TS_data[i][ii:ii+1,i_gene])/exp_time_diff[ii] + alpha_i*TS_data[i][ii:ii+1,i_gene])
                output_vect_time.append(f_dy_dt)
            
            exp_n_samples = exp_timeseries_x.shape[0]
            input_matrix_time[i*exp_n_samples:(i+1)*exp_n_samples,:] = exp_timeseries_x
            

        # Steady-state data (if any)
        if False: 
            pass
        
        else:
      
            input_all = input_matrix_time
            output_all = output_vect_time
            
        Xs.append(input_all)
        ys.append(output_all)

    return Xs,ys
def rand_search(n_sample = 100):
    """samples among all possible scenarios, and returns the best params, also the stats of each params.
    This can be then used to evaluate which params are more packt (more important), for further search by grid search. 
    """
    # TODO: before implementing this, be sure about how alpha will be included in the rate of prediction. 
    # read the relevant paper
def evaluate(Xs, ys, param, cv = 4):
    # TODO: This is for SA. for now, its not practical
    """ evaluetes the whole network and returns a test score which is the average score of all estimator"""
    scores = []
    for X, y in zip(Xs, ys):
        _, score = evaluate_single(X, y, param, cv)
        scores.append(score)
    return np.mean(scores)

def evaluate_single(X, y, param, cv = 4):
    """ evalutes the target estimator and returns a test score 
    for RF, the test score is oob_score. For the rest, cv is applied.
    """
    if param['estimator_t'] == 'RF' :
        use_oob_flag = True
    else:
        use_oob_flag = False
    
    est = TargetEstimator(**param)
    if use_oob_flag:
        fit = est.fit(X,y)
        score = est.test()
    else:
        cv_a = model_selection.ShuffleSplit(n_splits=cv, test_size=(1/cv))
        scores = model_selection.cross_val_score(est, X, y, cv=cv_a)
        score = np.mean(scores)

    return est, score
def grid_search_single_gene(X, y, param, permts, cv = 4):
    """ evalute one gene for the given params, and returns the best fit """
    # evalute each permutation
    fits = []
    scores = []
    for permt in permts:  
        param_a = {**param,**permt}
        fit, score = evaluate_single(X, y, param_a,cv)
        fits.append(fit)
        scores.append(score)
    # find the best candidate. Max score is considered best score. 
    best_score = max(scores)
    index = scores.index(best_score)
    best_param = permts[index]
    best_est = fits[index]

    return best_score, best_param, best_est
def grid_search_single_permut(Xs, ys, param, permt, cv = 4):
    """ evalute all genes for the one permutation of param, and returns the best fit """
    param_a = {**param,**permt}
    fits = []
    scores = []
    for X, y in zip(Xs, ys):
        fit, score = evaluate_single(X, y, param_a,cv)
        fits.append(fit)
        scores.append(score)
    return scores, fits

def map_gene(args):
    """ maps the args to the grid search function for single target, used for multi threating """
    i = args['i'] # index of each target
    args_rest = {key:value for key,value in args.items() if key != 'i'}
    return (i, grid_search_single_gene(**args_rest))

def map_permut(args):
    """ maps the args to the grid search function for single permutation of param, used for multi threating """
    i = args['i'] # index of each target
    args_rest = {key:value for key,value in args.items() if key != 'i'}
    return (i, grid_search_single_permut(**args_rest))

def grid_search(Xs,ys, param, param_grid, n_jobs = 1, **specs):
    """  """
    # TODO: verify the inputs
    n_genes = Xs[0].shape[1]
    time_start = time.time()
    # permuation of param grid
    print('Grid params:', param_grid)
    tags = list(param_grid.keys())
    values = list(param_grid.values())
    permts = []
    for value in list(itertools.product(*values)):
        permts.append({tag:i for tag,i in zip(tags,value)})
    print('%d genes %d permts -> %d jobs on %d threads' % (n_genes, len(permts), n_genes*len(permts), n_jobs))
    best_scores = [None for i in range(n_genes)]
    best_params = [None for i in range(n_genes)]
    best_ests = [None for i in range(n_genes)]
    if n_jobs == 1: # serial
        map_input =  [{'i':i, 'X': Xs[i], 'y': ys[i], 'param': param, 'permts': permts, **specs} for i in range(n_genes)]
        all_output = list(map(map_gene, map_input))
        all_output.sort(key = lambda x: x[0])
        for i_gene, output in enumerate(all_output): # each output is for a gene
            best_score, best_param, best_est = output[1]
            best_scores[i_gene] = best_score
            best_params[i_gene] = best_param
            best_ests[i_gene] = best_est
    else: # parallel
        # multithreading happens either on gene_n or permuts, depending which one is bigger
        pool = Pool(n_jobs)
        if n_genes >= len(permts):
            print('Gene-based multi threading')
            map_input =  [{'i':i, 'X': Xs[i], 'y': ys[i], 'param': param, 'permts': permts, **specs} for i in range(n_genes)]
            all_output = pool.map(map_gene, map_input)
            all_output.sort(key = lambda x: x[0])
            for i_gene,output in enumerate(all_output): # each output is for a gene
                best_score, best_param, best_est = output[1]

                best_scores[i_gene] = best_score
                best_params[i_gene] = best_param
                best_ests[i_gene] = best_est
        else: # when there is more permuts
            print('Permutation-based multi threading')
            input_data =  [{'i':i, 'Xs': Xs, 'ys': ys, 'param': param, 'permt': permts[i], **specs} for i in range(len(permts))]
            all_output = pool.map(map_permut, input_data)
            all_output.sort(key = lambda x: x[0])

            scores = np.empty([n_genes, 0])
            ests = np.empty([n_genes, 0])
            for output in all_output: # each output is for a permut
                scores_all_genes, ests_all_genes = output[1] # for all genes 
                scores = np.c_[scores, scores_all_genes]
                ests = np.c_[ests, ests_all_genes]
            best_scores = scores.max(1)
            best_indices = np.array(scores.argmax(1))
            best_ests = ests[range(n_genes),best_indices]
            best_params = [permts[i] for i in best_indices]


    time_end = time.time()
    print('Param search is completed in %.3f seconds' % ((time_end-time_start)))
    return best_scores, best_params, best_ests
def sklearn_search(Xs, ys, param, param_grid, specs): # should use oob score when specified, for RF
    '''  '''
    n_genes = len(ys)
    ests = [TargetEstimator(**param) for i in range(n_genes)]    

    res = []
    for est, X, y in zip(ests,Xs, ys):
        # search = model_selection.GridSearchCV(est, param_grid, **specs)
        search = model_selection.RandomSearchCV(est, param_grid, **specs)
        res.append(search.fit(X,y))
    best_scores = [r.best_score_ for r in res]
    best_params = [r.best_params_ for r in res] 
    best_fitted_est = [r.best_estimator_ for r in res]
    return best_scores, best_params, best_fitted_est  
def network_inference(Xs, ys, param  , param_unique = None, ests = None):
    """ Determines links of network inference
    If the ests are given, use them instead of creating new ones.

    Xs -- 
    """
    n_genes = len(ys)
    if ests == None:
        if param_unique == None:
            ests = [TargetEstimator(**param) for i in range(n_genes)]
        else: 
            ests = [TargetEstimator(**{**param,**param_unique[i]}) for i in range(n_genes)]  
        fits = [ests[i].fit(X,y) for i, (X, y) in enumerate(zip(Xs,ys))]
    else: # fitted estimators are given 
        assert(len(ests) == n_genes)
    scores_train = [ests[i].score(X,y) for i, (X, y) in enumerate(zip(Xs,ys))]
    if param['estimator_t'] == 'RF':
        oob_scores = [est.test() for est in ests]
        
    else:
        oob_scores = None

    # wights of links from all features to i target
    links = [ests[i].weight_links() for i in range(n_genes)]
    return ests, scores_train, links, oob_scores
class TargetEstimator(base.BaseEstimator,base.RegressorMixin):
    """The docstring for a class should summarize its behavior and list the public methods and instance variables """
    def __init__(self,estimator_t, alpha = 0, **params):
        '''args should all be keyword arguments with a default value -> kwargs should be all the keyword params of all regressors with values'''
        '''they should not be documented under the “Attributes” section, but rather under the “Parameters” section for that estimator.'''
        '''every keyword argument accepted by __init__ should correspond to an attribute on the instance'''
        '''There should be no logic, not even input validation, and the parameters should not be changed. The corresponding logic should be put where the parameters are used, typically in fit'''
        '''algorithm-specific unit tests,'''
        # self.alpha = alpha
        self.params = params
        self.estimator_t = estimator_t
        self.alpha = alpha
        # self._required_parameters = () #estimators also need to declare any non-optional parameters to __init__ in the
    def fit(self,X,y):
        """ fit X to y
        X -- Array-like of shape (n_samples, n_features)
        y -- Array-like of shape (n_samples,)
        kwargs -- Optional data-dependent parameters
        """
        '''Attributes that have been estimated from the data must always have a name ending with trailing underscore'''
        '''The estimated attributes are expected to be overridden when you call fit a second time.'''
        
        # apply alpha to y
        y = [y_i(self.alpha) for y_i in y]

        utils.check_array(X)
        utils.check_X_y(X,y)
        utils.indexable(X)
        utils.indexable(y)

        if self.estimator_t != 'HGB': #check this. https://scikit-learn.org/stable/developers/utilities.html#developers-utils
            utils.assert_all_finite(X)
            utils.assert_all_finite(y)
        self.X_ = X
        self.y_ = y

        if self.estimator_t == 'RF':
            self.est = ensemble.RandomForestRegressor(oob_score = True,**self.params)
        elif self.estimator_t == 'HGB':
            self.est = ensemble.HistGradientBoostingRegressor(**self.params)
        else:
            raise ValueError('Define estimator_t')
        self.est.fit(X,y)
        return self
    def predict(self,X):
        """ """
        # apply alpha to y
        y = [y_i(self.alpha) for y_i in y]
        utils.validation.check_is_fitted(self)
        utils.check_X_y(X,y)
        return self.est.predict(X)

    def score(self,X,y): 
        """ """
        # apply alpha to y
        y = [y_i(self.alpha) for y_i in y]
        utils.validation.check_is_fitted(self)
        utils.check_array(X)
        utils.check_X_y(X,y)
        utils.indexable(X)
        utils.indexable(y)
        return self.est.score(X,y)
    def test(self,X = None, y = None):
        """ evaluate the fit
        if RF, it returns back oob_score
        """
        utils.validation.check_is_fitted(self)
        if self.estimator_t == 'RF':
            return np.mean(self.est.oob_score_)
        else:
            return self.score(X,y)


    def weight_links(self):
        """ """
        vi = compute_feature_importances(self.est)
        # Normalize importance scores
        vi_sum = sum(vi)
        if vi_sum > 0:
            vi = vi / vi_sum
        return vi
    # def transform(self):
    #     return
    # def get_params(self,deep=True):
    #     #The get_params function takes no arguments and returns a dict of the __init__ parameters of the estimator, together with their values.
    #     #get_params mechanism is not essential
    #     return {'regressor'}

    def set_params(self, **parameters):
        """ """
        for parameter, value in parameters.items():
            setattr(self, parameter, value)
        return self
    def _more_tags(self):
        """ """
        if self.estimator_t == 'HGB':
            allow_nan = True 
        else:
            allow_nan = False
        return {'requires_fit': True, 'allow_nan': allow_nan, 'multioutput': True, 
            'requires_y': True,}

    
   