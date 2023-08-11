from pathlib import Path
import os
import json
import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, TypeAlias, Any, Callable
from abc import ABC, abstractmethod
import pickle
from scipy.spatial.distance import jensenshannon
from tqdm import tqdm

from proteomics_MSC.VSA.vsa_aux import retrieve_VSA_results
from proteomics_MSC.common_tools.links import run_portia
from proteomics_MSC.common_tools import flatten, read_write_data
from proteomics_MSC.imports import DIVERGENCE_ANALYSIS_DIR, DIVERGENCE_ANALYSIS_ROBUSTNESS_DIR
from proteomics_MSC.common_tools import F_model_name_2_method_and_DE_type, F_DE_genenames, F_selected_models
ASPS_series_type: TypeAlias = List[Tuple[float, float]]  # data type for list of (Q,P)

def add_noise(data_org, noise_type, noise_intensity):

    if noise_type == 'mp':
        noisy_data = _add_mpnoise(data_org, std=noise_intensity)
    else:
        noisy_data = _add_adnoise(data_org, std=noise_intensity)

    return noisy_data
def create_noisy_data(n_repeat, **kwargs) -> List[pd.DataFrame]:
    """Create n number of random noise"""
    noisy_sets = []
    for i in range(n_repeat):
        noisy_data = add_noise(**kwargs)
        noisy_sets.append(noisy_data)
    return noisy_sets
def create_noisy_links(data_org:pd.DataFrame, n_repeat:int, noise_type:str,
                       noise_intensity:float,  gene_names: List[str],
                       grn_function: Callable[[pd.DataFrame, List[str]], pd.DataFrame]) -> List[pd.DataFrame]:
    """
    add noise to the data, run grn, and save it. If file already exist, skip it unless forced.
    Collect the links and retunt a list.
    """

    # - create n noisy dataset
    noisy_datasets = create_noisy_data(data_org=data_org, noise_type=noise_type, noise_intensity=noise_intensity, n_repeat=n_repeat)
    # - run grn for each noisy link
    noisy_links = []
    with tqdm(total=n_repeat, desc=f'Run GRN for noisy data: {noise_type} {noise_intensity}') as pbar:
        for data_i, data in enumerate(noisy_datasets):
            # - run GRN
            links = grn_function(data, gene_names)
            noisy_links.append(links)
            pbar.update(1)
    return noisy_links



def _add_mpnoise(data: np.ndarray, std=.3) -> np.array:
    """ Multiplicative noise
         Creates n_noisy_datasets data
    """
    data_std = np.std(data)
    applied_std = std * np.array(data_std)

    frac_noise = np.random.rand()
    noise_mask = np.random.choice([0, 1], size=data.shape, p=[1 - frac_noise, frac_noise])
    rand_values = np.random.normal(loc=1, scale=applied_std, size=data.shape)
    noisy_data = data * (1 + (rand_values - 1) * noise_mask)


    return noisy_data

def _add_adnoise(data: np.ndarray, std=.05) -> np.array:
    """ Additive noise
             Creates n_relica noised links
    """
    data_std = np.std(data)
    applied_std = data_std * std
    rand_values = np.random.normal(loc=1, scale=applied_std, size=data.shape)
    noisy_data = rand_values + data

    return noisy_data


class NoiseAnalysis(ABC):
    """Main base class to conduct noise analysis"""
    def __init__(self, model_name, studies, base_dir, force, n_repeat=500):
        self.model_name = model_name
        self.force = force
        self.studies = studies
        self.base_dir = base_dir
        self.noise_types = ['mp', 'ad']
        self.noise_intensities = {'mp': [0.05], 'ad': [.1]}
        self.robust_thresholds = {'mp': 0.05, 'ad': .1}

        self.n_repeat = n_repeat
        # - create the base dir
        if not os.path.isdir(Path(self.base_dir)):
            os.makedirs(Path(self.base_dir))
    @abstractmethod
    def _run_analysis(self, model_name, noisy_links_stack_ctr, noisy_links_stack_sample, save_file):
        pass
    @abstractmethod
    def _run_post_analysis(self, model_name: str, noise_type: str, noise_intensity: float, save_file: str):
        pass
    @abstractmethod
    def _save_file(self, noise_type, noise_intensity):
        """The file name that the analysis results is saved to"""
        pass
    def _save_dir(self):
        """
            Assigning a unique folder to the model
        """
        results_save_dir = Path(self.base_dir) / self.model_name
        if not os.path.isdir(results_save_dir):
            os.makedirs(results_save_dir)
        return results_save_dir
    def _save_plot_dir(self):
        """
            Assigns a unique folder for the plots
        """
        results_save_dir = self._save_dir()
        plots_save_dir = Path(results_save_dir) / 'plots'
        if not os.path.isdir(plots_save_dir):
            os.makedirs(plots_save_dir)
        return plots_save_dir
    @staticmethod
    def _create_noisy_links_studies(model_name, studies, **kwargs_noise):
        """
            Gets the data and runs create_noisy_links for studies
        """
        method, DE_type = F_model_name_2_method_and_DE_type(model_name)
        genenames = F_DE_genenames()[DE_type]
        # - get the data
        data_studies = [read_write_data(mode='read', tag=f'{DE_type}_{study}') for study in studies]
        # - create noisy links for ctr and sample
        kwargs_grn = dict(gene_names=genenames, grn_function=run_portia)
        kwargs = {**kwargs_grn, **kwargs_noise}
        noisy_links_stack_ctr, noisy_links_stack_sample = [create_noisy_links(data_org=data, **kwargs) for data in
                                                           data_studies]
        return noisy_links_stack_ctr, noisy_links_stack_sample
    def run_decorator(self) -> None:
        """
            Runs the main analysis functions for each noise type and noise intensity
        """
        # - run for each noise type
        for noise_type in self.noise_types:
            # - run for each noise intensity
            for noise_intensity in self.noise_intensities[noise_type]:
                # - run the analysis
                self._run_analysis(noise_type, noise_intensity)
                self._run_post_analysis(noise_type, noise_intensity)



class NetworkNoiseAnalysis(NoiseAnalysis):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.links_template = None
        self.divergence_scores = {}

    def _save_file(self, noise_type, noise_intensity) -> str:
        save_dir = self._save_dir()
        return f'{save_dir}/links_{noise_type}_{noise_intensity}.pkl'  # file to save the results

    @classmethod
    def __pool_links(cls, links_stack: List[pd.DataFrame]) -> List[List[float]]:
        """Gets list of links and pool the weights"""

        # - create pool of weights
        pool_weights = [links['Weight'].values.tolist() for links in links_stack]
        pool_weights: List[List[float]] = [list(row) for row in zip(*pool_weights)] # transpose
        return pool_weights

    def _run_analysis(self, noise_type: str, noise_intensity: float) -> None:
        save_file = self._save_file(noise_type, noise_intensity)
        studies = self.studies
        # - run the analysis if the output file doesnt exist
        if not os.path.exists(save_file) or self.force:
            print(f'Noise analysis for {self.model_name} {noise_type} {noise_intensity}')
            # - run grn to get the links
            noisy_links_stack_ctr, noisy_links_stack_sample = self._create_noisy_links_studies(
                model_name=self.model_name,
                studies=studies,
                noise_type=noise_type,
                noise_intensity=noise_intensity,
                n_repeat=self.n_repeat)
            results = {'ctr': noisy_links_stack_ctr, 'sample': noisy_links_stack_sample}

            with open(save_file, 'wb') as ff:
                pickle.dump(results, ff)

    def _run_post_analysis(self, noise_type: str, noise_intensity: float) -> None:
        """The analysis of the results and plots.
            Get the links, pool them for control and sample, calculate distance
        """
        # - retrieve the noisy links
        save_file = self._save_file(noise_type, noise_intensity)
        with open(save_file, 'rb') as ff:
            results = pickle.load(ff)
        noisy_links_stack_ctr, noisy_links_stack_sample = results['ctr'], results['sample']
        # - extract the links as template
        self.links_template = noisy_links_stack_ctr[0].copy()
        self.links_template.drop('Weight', axis=1, inplace=True)
        # - pool the weights
        pool_ctr = self.__pool_links(noisy_links_stack_ctr)
        pool_sample = self.__pool_links(noisy_links_stack_sample)
        # - calculate divergence
        distances = []  # between ctr and sample
        for row_ctr, row_sample in zip(pool_ctr, pool_sample):
            # std_ctr, std_ = np.std(row_ctr)
            js_divergence = jensenshannon(row_ctr, row_sample)
            js_distance = np.sqrt(js_divergence)  # between 0 and 1
            distances.append(js_distance)
        # - store it
        if noise_type not in self.divergence_scores:
            self.divergence_scores[noise_type] = {}
        self.divergence_scores[noise_type][noise_intensity] = distances

    def determine_overall_divergence_scores(self) -> pd.DataFrame:
        """ Determine overall distance score between control and sample for a given model
            It takes into account the distances calculated for each noise type and internsity.
            For each link, a distance score is calculated to quantify the degree of change from ctr to sample.
            It prints the results to a file
        """
        assert (self.divergence_scores != {})
        assert (self.links_template is not None)

        _distances = []
        for noise_type, noise_type_data in self.divergence_scores.items():
            threshold = self.robust_thresholds[noise_type]  # the noise intensity that we consider for robustness
            distances_noisetype = noise_type_data[threshold]
            _distances.append(distances_noisetype)
        # - final distance is root sum of both noise types
        _distances = np.asarray(_distances)
        _distances = _distances ** 2
        divergence_scores = np.sum(_distances, axis=0)
        # - put the score into a template
        divergence_scores_df = self.links_template.copy()
        divergence_scores_df['divergence_score'] = divergence_scores
        return divergence_scores_df
