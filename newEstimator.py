"""This is the example module.

This module does stuff.
"""
__all__ = ['a', 'b', 'c']
__version__ = '0.1'
__author__ = 'Jalil Nourisa'

import time
import numpy as np
import itertools

from numpy.random import permutation, uniform
from operator import itemgetter
from multiprocessing import Pool
import matplotlib.pyplot as plt

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
    vInter_sort = sorted(vInter,key=itemgetter(2),reverse=True)
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
        items_perm = random.permutation(items_perm)
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
def estimate_degradation_rates(TS_data,time_points):
    
    """ Estimation of degradation rates from data
    For each gene, the degradation rate is estimated by assuming that the gene expression x(t) follows:
    x(t) =  A exp(-alpha * t) + C_min,
    between the highest and lowest expression values.
    C_min is set to the minimum expression value over all genes and all samples.

    keyword arguments:
    TS_data -- 
    time_points -- 
    """
    
    ngenes = TS_data[0].shape[1]
    nexp = len(TS_data)
    
    C_min = TS_data[0].min()
    if nexp > 1:
        for current_timeseries in TS_data[1:]:
            C_min = min(C_min,current_timeseries.min())
    
    alphas = np.zeros((nexp,ngenes))
    
    for (i,current_timeseries) in enumerate(TS_data):
        current_time_points = time_points[i]
        
        for j in range(ngenes):
            
            idx_min = np.argmin(current_timeseries[:,j])
            idx_max = np.argmax(current_timeseries[:,j])
            
            xmin = current_timeseries[idx_min,j]
            xmax = current_timeseries[idx_max,j]
            
            tmin = current_time_points[idx_min]
            tmax = current_time_points[idx_max]
            
            xmin = max(xmin-C_min,1e-6)
            xmax = max(xmax-C_min,1e-6)
                
            xmin = np.log(xmin)
            xmax = np.log(xmax)
            
            alphas[i,j] = (xmax - xmin) / abs(tmin - tmax)
                
    alphas = alphas.max(axis=0)
    return alphas

def process_data(TS_data, time_points, regulators = 'all',alpha='from_data'): 
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

    if alpha == 'from_data':
        alphas = estimate_degradation_rates(TS_data,time_points)
    elif isinstance(alpha,(int,float)):
        alphas = zeros(ngenes) + float(alpha)    
    else:
        alphas = [float(a) for a in alpha]

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
        output_vect_time = np.zeros(nsamples_time-h*nexp)

        nsamples_count = 0
        for (i,current_timeseries) in enumerate(TS_data):
            current_time_points = time_points[i]
            npoints = current_timeseries.shape[0]
            time_diff_current = current_time_points[h:] - current_time_points[:npoints-h]
            current_timeseries_input = current_timeseries[:npoints-h,input_idx]
            current_timeseries_output = (current_timeseries[h:,i_gene] - current_timeseries[:npoints-h,i_gene]) / time_diff_current + alphas[i_gene]*current_timeseries[:npoints-h,i_gene]
            nsamples_current = current_timeseries_input.shape[0]
            input_matrix_time[nsamples_count:nsamples_count+nsamples_current,:] = current_timeseries_input
            output_vect_time[nsamples_count:nsamples_count+nsamples_current] = current_timeseries_output
            nsamples_count += nsamples_current
        
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
def grid_search_single(X, y, param, param_grid, cv = 4, use_oob_score = False):
    """ """
    # TODO: check that estimator_t is given
    # TODO: check that param and param_grid dont share same params
    # TODO: if use_oob_score is true, param should have oob_score given
    use_oob_flag = False
    if param['estimator_t'] == 'RF' and use_oob_score == True:
        use_oob_flag = True

    # permuation of param grid
    tags = list(param_grid.keys())
    values = list(param_grid.values())
    permts = []
    for value in list(itertools.product(*values)):
        permts.append({tag:i for tag,i in zip(tags,value)})

    # evalute each permutation
    fits = []
    scores = []
    for permt in permts:  
        param_a = {**param,**permt}
        est = TargetEstimator(**param_a)
        if use_oob_flag:
            fit = est.fit(X,y)
            p_scores = fit.est.oob_score_
        else:
            # TODO: use CV
            raise ValueError('Define first')
        score_m = np.mean(p_scores)

        fits.append(fit)
        scores.append(score_m)
    # find the best candidate. Max score is considered best score. 
    best_score = max(scores)
    index = scores.index(best_score)
    best_param = permts[index]
    best_est = fits[index]

    return best_score, best_param, best_est
def pool_map_function(args):
    """ maps the args to the grid search function, used for multi threating """
    i = args['i'] # index of each target
    args_rest = {key:value for key,value in args.items() if key != 'i'}
    return (i, grid_search_single(**args_rest))
def grid_search(Xs,ys, param, param_grid, n_jobs = 1, **specs):
    """  """
    # TODO: verify the inputs
    n_genes = Xs[0].shape[1]
    print('Running param search on %d threads' % n_jobs)
    best_scores = [None for i in range(n_genes)]
    best_params = [None for i in range(n_genes)]
    best_ests = [None for i in range(n_genes)]
    if n_jobs == 1: # serial
        for i_gene in range(n_genes):
            X = Xs[i_gene]
            y = ys[i_gene]
            best_score, best_param, best_est = grid_search_single(X, y, param, param_grid, **specs)
            best_scores[i_gene] = best_score
            best_params[i_gene] = best_param
            best_ests[i_gene] = best_est
    else: # parallel
        pool = Pool(n_jobs)
        input_data =  [{'i':i, 'X': Xs[i], 'y': ys[i], 'param': param, 'param_grid': param_grid, **specs} for i in range(n_genes)]
        all_output = pool.map(pool_map_function, input_data)
        for output in all_output:
            i_gene = output[0]
            best_score, best_param, best_est = output[1]

            best_scores[i_gene] = best_score
            best_params[i_gene] = best_param
            best_ests[i_gene] = best_est

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
def network_inference(Xs, ys, param , param_unique = None, ests = None):
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
        if 'oob_score' in param:
            oob_scores = [est.est.oob_score_ for est in ests]
        else:
            print('oob_score is not specified in RF')
            oob_scores = None 
    else:
        oob_scores = None

    # wights of links from all features to i target
    links = [ests[i].weight_links() for i in range(n_genes)]
    return ests, scores_train, links, oob_scores
class TargetEstimator(base.BaseEstimator,base.RegressorMixin):
    """The docstring for a class should summarize its behavior and list the public methods and instance variables """
    def __init__(self,estimator_t,**params):
        '''args should all be keyword arguments with a default value -> kwargs should be all the keyword params of all regressors with values'''
        '''they should not be documented under the “Attributes” section, but rather under the “Parameters” section for that estimator.'''
        '''every keyword argument accepted by __init__ should correspond to an attribute on the instance'''
        '''There should be no logic, not even input validation, and the parameters should not be changed. The corresponding logic should be put where the parameters are used, typically in fit'''
        '''algorithm-specific unit tests,'''
        # self.alpha = alpha
        self.params = params
        self.estimator_t = estimator_t
        # self._required_parameters = () #estimators also need to declare any non-optional parameters to __init__ in the
    def fit(self,X,y, alpha = 0):
        """ fit X to y
        X -- Array-like of shape (n_samples, n_features)
        y -- Array-like of shape (n_samples,)
        kwargs -- Optional data-dependent parameters
        """
        '''Attributes that have been estimated from the data must always have a name ending with trailing underscore'''
        '''The estimated attributes are expected to be overridden when you call fit a second time.'''
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
            self.est = ensemble.RandomForestRegressor(**self.params)
        else:
            raise ValueError('Define estimator_t')
        self.est.fit(X,y)
        return self
    def predict(self,X):
        """ """
        utils.validation.check_is_fitted(self)
        utils.check_X_y(X,y)
        return self.est.predict(X)

    def score(self,X,y): 
        """ """
        utils.validation.check_is_fitted(self)
        utils.check_array(X)
        utils.check_X_y(X,y)
        utils.indexable(X)
        utils.indexable(y)
        return self.est.score(X,y)
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


