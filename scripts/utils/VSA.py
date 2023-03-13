"""
    Sets of functions useful vester's sensitivity analysis
"""
import sys
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import math
from typing import List, Tuple
from scipy.stats import ttest_ind
import scipy
from tqdm import tqdm


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from scripts.utils import serif_font, comic_font
from scripts.utils.links import choose_top_quantile


def __Q_P_thresholds(AS, PS) -> Tuple[float,float]: #TODO:deprecated
    '''
        Calculate the threshold to divide AS and PS into 2
    '''
    # def pool(df_1, df_2, key):
    #     return list(df_1[key])+list(df_2[key])+
    # PS = pool(df_ctr, df_sample, 'PS')
    # AS = pool(df_ctr, df_sample, 'AS')
    # PS_t, AS_t = np.mean(PS), np.mean(AS)
    PS_t, AS_t = np.quantile(PS,.75), np.quantile(AS,.75)
    return AS_t, PS_t
def role_analysis(links, gene_names) -> pd.DataFrame:
    '''
        Vester's sensitivity analysis. 
            - AS: active sum. Sum along rows of the influence matrix and it indicates how much does a variable influence all the others.
            - PS: passive sum. Its is the sum along columns of the influence matrix and it indicates how sensitive a variable is, how does it react to the influence of others
            - Q: AS/PS -> how dominant
            - P: AS.PS -> how participative a variable is
            - Active: +Q
            - Passive: -Q, -P
            - Critical: +Q, +P
            - Buffering: -Q, -P
    ------------------------------------------------------
    inputs: links (DataFrame)
    output: VSA (DataFrame) -> protname: AS, PS, Role
    '''
    #- active sum and passive sum
    links = choose_top_quantile(links, 0.5)
    AS = np.asarray([sum(links.query(f"Regulator == '{gene}'")['Weight']) for gene in gene_names])
    PS = np.asarray([sum(links.query(f"Target == '{gene}'")['Weight']) for gene in gene_names])

    #- Quotient and Product
    Q, P = AS/PS, PS*PS

    #- standardize
    def standardize_array(arr):
        # Calculate the 0.75 quantile of the array
        q75 = np.quantile(arr, 0.75)

        # Subtract the 0.75 quantile from each value in the array
        arr = arr - q75

        # Calculate half of the absolute maximum value of the resulting array
        max_abs_val = np.max(np.abs(arr)) / 2

        # Divide each value in the resulting array by half of the absolute maximum value
        arr = arr / max_abs_val

        # Multiply each value in the resulting array by 0.5 to shift the range from [-1,1] to [-0.5,0.5]
        arr = arr * 0.5

        return arr


    #- define the roles: ['Buffering', 'Passive', 'Active', 'Critical'] -> [0, 1, 2, 3] 
    # Q_t, P_t = Q_P_thresholds(Q, P)
    # Q = Q - Q_t
    # P = P - P_t
    Q, P = standardize_array(Q), standardize_array(P)
    Q_flags = Q>0
    P_flags = P>0
    roles = [2*Q_flag+P_flag for Q_flag, P_flag in zip(Q_flags,P_flags)]
    vsa_results = pd.DataFrame(data={'Entry': gene_names, 'Q': Q, 'P': P, 'Role':roles })
    return vsa_results
def determine_critical_role_change(df_ctr, df_sample, target_role:int=3):
    '''
        Finds those genes with a change in the role from ctr to sample.

    Inputs: target_role -> critical:3
    '''
    # - find those genes in both ctr and sample that have target role, i.e. critical
    indices_target = list(map(lambda df: df.index[df['Role']==target_role].values, [df_ctr, df_sample]))
    indices_target = [i for j in indices_target for i in j] #flatten
    indices_target = list(set(indices_target))

    df_ctr_c = df_ctr.loc[indices_target,:].reset_index(drop=True)
    df_sample_c = df_sample.loc[indices_target,:].reset_index(drop=True)

    df_role_change = pd.DataFrame()
    df_role_change['Entry'] = df_ctr_c['Entry']
    df_role_change[['Q1','P1']]=df_ctr_c[['Q','P']]
    df_role_change[['Q2','P2']]=df_sample_c[['Q','P']]
    return df_role_change
def determine_top_role_change(control: pd.DataFrame, sample: pd.DataFrame, top_quantile: float = 0.95) -> pd.DataFrame:

    '''
    Find the distance in the role from control to sample and return the top_quantile of the distances.

    Args:
        control (pd.DataFrame): The control DataFrame with columns 'Entry', 'P', and 'Q'.
        sample (pd.DataFrame): The sample DataFrame with columns 'Entry', 'P', and 'Q'.
        top_quantile (float): The quantile to use for selecting the top distances.

    Returns:
        pd.DataFrame: A DataFrame with columns 'Entry' and 'Distance', containing the entries with the top distances.

    '''
    # - calculate distance
    distance = (sample['Q'] - control['Q'])**2 + (sample['P'] - control['P'])**2
    df_distance = pd.DataFrame({'Entry': control['Entry'], 'Distance':distance})
    df_distance = df_distance.assign(P1=control['P'], Q1=control['Q'])
    df_distance = df_distance.assign(P2=sample['P'], Q2=sample['Q'])

    #-
    cut_off = np.quantile(df_distance['Distance'].values.tolist(), q=top_quantile)
    short_list = df_distance.loc[df_distance['Distance'] >= cut_off, :].reset_index(drop=True)
    return short_list
class RolePlot:
    """
    Plots VSA results
    """
    linewidth = 2
    markers = ['o', 'o', 'o', '*']
    roles_colors = ['lightgreen', 'lightblue', 'yellow', 'red']
    ctr_sample_colors = ['cyan', 'purple']
    sizes = [40 * i for i in [1, 1.5, 2, 5]]
    roles = ['Buffering', 'Passive', 'Active', 'Critical']
    @staticmethod
    def create_role_legends(ax):
        handles = []
        for i, color in enumerate(RolePlot.roles_colors):
            handles.append(ax.scatter([], [], marker='o', label=RolePlot.roles[i], s=100,
                                                  edgecolor='black', color=color, linewidth=.5))
        return handles
    def plot_role_change(df, ax=None):
        '''
            Plots Q/P for role change betwen ctr and sample
        '''
        xlim = ax.get_xlim()[1] - ax.get_xlim()[0]
        ylim = ax.get_ylim()[1] - ax.get_ylim()[0]
        # an arrow from ctr to mg only for any change from or to critical role
        # names of ctr or mg on the symbols
        # - first ctr
        for i in [1,2]:
            ax.scatter(df[f'P{i}'], df[f'Q{i}'],
                            color=RolePlot.ctr_sample_colors[i-1],
                            alpha=1,
                            linewidths=.2,
                            edgecolors='black',
                            s=100
                            )

        for prot, x0, y0, x1, y1 in zip(df['Entry'].values, df['P1'].values, df['Q1'].values, df['P2'].values,
                                        df['Q2'].values):
            offset_x = .015
            offset_y = 0.1
            arrow_t = 'arc3'
            rad = .3
            if prot == 'Q07866':
                offset_x = -.045
                offset_y = -0.3
                rad = -.1
            if prot == 'Q99613':
                offset_x = -.05
                offset_y = -0.3
                # arrow_t = 'rad2'
            if prot == 'Q02790':
                # offset_x = -.05
                offset_y = 0
                # arrow_t = 'rad2'
            if prot == 'P00568':
                offset_x = -6.5
                offset_y = -.15
                # arrow_t = 'rad2'
            if prot == 'Q02218':
                offset_x = -.2
                offset_y = .05
                # arrow_t = 'rad2'
            if prot == 'P02652':
                offset_x = -.5
                offset_y = .05
                # arrow_t = 'rad2'



            x, y = x0 + offset_x * xlim, y0 + offset_y * ylim
            ax.annotate(f'{prot}', xy=(x, y), fontsize=9)

            ax.annotate('', xy=(x1, y1), xytext=(x0, y0),
                        arrowprops={'arrowstyle': '->', 'connectionstyle': f'{arrow_t}, rad={rad}'}
                        , horizontalalignment='center')

        RolePlot.mark_x(ax, 0)
        RolePlot.mark_y(ax, 0)
        RolePlot.plot_roles(ax)
        RolePlot.postprocess(ax, title='Role change')

    def plot_ctr_vs_sample(ax, data, gene_names):
        '''
            Plots Q/P for ctr and mg conditions in a 1*2 subplot
        '''
        comic_font()


        RolePlot.scatter(ax, P=data['P'], Q=data['Q'], roles=data['Role'])
        # Q_t, P_t = RolePlot.thresholds(df['Q'], df['P'])
        RolePlot.mark_x(ax, 0)
        RolePlot.mark_y(ax, 0)
        RolePlot.plot_prot_names(ax=ax, Q=data['Q'], P=data['P'], gene_names=gene_names, roles=data['Role'])
        RolePlot.postprocess(ax)
        RolePlot.plot_roles(ax)

        # handles = []
        # for i, color in enumerate(RolePlot.roles_colors):
        #     handles.append(plt.scatter([], [], marker='o', label=RolePlot.roles[i],
        #                               edgecolor='black', color=color, linewidth=.2))
        # ll = axes[0].legend(handles=handles,
        #                 bbox_to_anchor=(1.4, 1), prop={'size': 8}
        #                 # title='Enriched Term'
        #                 )
        # axes[0].add_artist(ll)


    def postprocess(ax, title='', show_axis_names=True):
        # serif_font()
        if show_axis_names:
            # ax.set_xlabel('P (the line marks top 25 %)')
            # ax.set_ylabel('Q (the line marks top 25 %)')
            ax.set_xlabel('Product', fontweight='bold')
            ax.set_ylabel('Quotient', fontweight='bold')

        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        xlen, ylen = xlim[1] - xlim[0], ylim[1] - ylim[0]
        ax.set_title(title, fontweight='bold')
        ax.margins(.2)
        ax.set_xlim([-1.1, 1.1])
        ax.set_ylim([-1.1, 1.1])
        ax.set_xticks([])
        ax.set_yticks([])

    def mark_y(ax, y_t):
        """
        Mark threshold line for the sig
        """
        line_color = 'grey'
        dx = ax.get_xlim()
        ax.axhline(y_t, color='grey', linestyle='-.')

    @staticmethod
    def mark_x(ax, x_t):
        """
        Mark threshold line for fold change
        """
        line_color = 'grey'
        ax.axvline(x_t, color=line_color, linestyle='--', linewidth=RolePlot.linewidth)

    @staticmethod
    def scatter(ax, P, Q, roles):

        sc = ax.scatter(P, Q,
                        color=[RolePlot.roles_colors[i] for i in roles],
                        #                    alpha=[1 if flag else .5 for flag in flags]
                        #                    alpha = .7,
                        linewidths=.1,
                        edgecolors='black',
                        s=100
                        )
        # markers = [RolePlot.markers[i] for i in roles]
        # paths = []
        # for marker in markers:
        #     marker_obj = mmarkers.MarkerStyle(marker)
        #     path = marker_obj.get_path().transformed(
        #                 marker_obj.get_transform())
        #     paths.append(path)
        # sc.set_paths(paths)

    @staticmethod
    def plot_roles(ax):
        '''
            Plots the roles on the corners of the scatter plot
        '''
        # fontsize = 9

        ax.text(.01, .99, 'Active', ha='left', va='top', transform=ax.transAxes, fontweight='bold')
        ax.text(.01, .01, 'Buffering', ha='left', va='bottom', transform=ax.transAxes, fontweight='bold')
        ax.text(.99, .01, 'Passive', ha='right', va='bottom', transform=ax.transAxes, fontweight='bold')
        ax.text(.99, .99, 'Critical', ha='right', va='top', transform=ax.transAxes, fontweight='bold')

    @staticmethod
    def plot_prot_names(ax, Q, P, gene_names, roles):
        '''
            Prints names of critical genes on the scatter plot
        '''
        xlim = ax.get_xlim()[1] - ax.get_xlim()[0]
        ylim = ax.get_ylim()[1] - ax.get_ylim()[0]

        arrow_specs = {'arrow_type': 'arc3', 'rad': .2}
        for i, gene_name in enumerate(gene_names):
            offset_x = 0.15
            offset_y = 0.15
            rad = .2
            role = roles[i]
            if role != 3:  # just critical
                continue
            x = P[i]
            y = Q[i]
            if gene_name in ['Q02790']:
                # offset_y = 0
                offset_x = 0.05
                rad = 0
            if gene_name in ['P13667']:
                # offset_y = 0
                offset_x = -0.05
                rad = -.2
            if gene_name in ['Q99613']:
                offset_y = 0.25
                offset_x = +0.07
                rad = -.2
            if gene_name in ['Q07866']:
                offset_y = 0.15
                offset_x = +0.1
                # rad = 0
            if gene_name in ['P14174']:
                offset_y = 0.05
                offset_x = +0.2
            if gene_name in ['Q02218']:
                offset_y = 0.2
                offset_x = +0.1
            if gene_name in ['P00568']:
                offset_y = 0.3
                offset_x = -.15
                rad = -.3
            if gene_name in ['P02652']:
                offset_y = 0.3
                offset_x = +0.1
                rad = -.2



            ax.annotate(gene_name, xy=(x, y), xytext=(x + offset_x * xlim, y + offset_y * ylim), fontsize=9,
                        arrowprops={'arrowstyle': '->', 'lw': 0.9, 'connectionstyle': f'arc3, rad={rad}',
                                    'color': 'black', 'alpha': .7}
                        , horizontalalignment='center')
            # ax.annotate(gene_name, (x+offset_x,y+offset_y), fontsize=fontsize)



class NoiseAnalysis:
    def __init__(self, data_ctr, data_sample, target_genes, run_grn_func, kwargs_grn_func, kwargs_role_analysis):
        self.data_ctr = data_ctr
        self.data_sample = data_sample
        self.target_genes = target_genes
        self.run_grn_func = run_grn_func
        self.kwargs_grn_func = kwargs_grn_func
        self.kwargs_role_analysis = kwargs_role_analysis
        self.studies = ['ctr','sample']
        self.results_mp = None
        self.results_ad = None
        self.std_mpnoise = None
        self.std_adnoise = None
    def analyse_mp_noise(self,  n_rep=100, sigma_noise=1, std_noise = .05):
        # - multiplicative noise
        self.std_mpnoise = std_noise
        self.results_mp = self._analyse_noise(n_rep=n_rep, sigma_noise=sigma_noise, std_noise = std_noise, noise_type='mp')
        return self.results_mp
    def analyse_ad_noise(self,  n_rep=100, sigma_noise=1, std_noise = .05):
        # - additive noise
        self.std_adnoise = std_noise
        self.results_ad = self._analyse_noise(n_rep=n_rep, sigma_noise=sigma_noise, std_noise = std_noise, noise_type='ad')
        return self.results_ad
    def _analyse_noise(self,  n_rep=100, sigma_noise=1, std_noise = .05, noise_type='mp'):
        assert noise_type in ['mp','ad']
        vsa_results_stack_studies = []
        for i_data, data in enumerate([self.data_ctr, self.data_sample]):
            if noise_type == 'mp':
                noisy_data_stack = self._add_mpnoise(data, n_rep=n_rep, sigma=sigma_noise, std=std_noise)
            else:
                noisy_data_stack = self._add_adnoise(data, n_rep=n_rep, sigma=sigma_noise, std=std_noise)
            noisy_links_stack = []
            with tqdm(total=n_rep, desc='Run GRN for noisy data') as pbar:
                for data in noisy_data_stack:
                    noisy_links_stack.append(self.run_grn_func(data=data, i_data=i_data, **self.kwargs_grn_func))
                    pbar.update(1)
            vsa_results_stack = [role_analysis(links, **self.kwargs_role_analysis) for links in noisy_links_stack]
            vsa_results_stack_filtered = [oo.loc[oo['Entry'].isin(self.target_genes), :] for oo in vsa_results_stack]
            vsa_results_stack_studies.append(vsa_results_stack_filtered)
        results = self._organized(vsa_results_stack_studies[0], vsa_results_stack_studies[1])
        return results
    def _plot(self, axes, data, genenames, preferred_titles=None):
        '''
            Plots Q/P for ctr and mg conditions, for n different runs of sensitivity analysis, one window for each prot

            oo_prots: Q/P data for ctr and sample, for each gene
        '''
        assert len(data) == len(genenames)
        comic_font()
        arrow_specs = {'arrow_type': 'arc3', 'rad': .2}

        for idx, (gene, data) in enumerate(zip(genenames, data)): #- plot for each prot on a seperate window
            ax = axes[idx]
            for idxx, (study, study_data) in enumerate(data.items()):
                Ps = [item[1] for item in study_data]
                Qs = [item[0] for item in study_data]
                ax.scatter(Ps, Qs,
                                color=RolePlot.ctr_sample_colors[idxx],
                                alpha=.6,
                                linewidths=.2,
                                edgecolors='black',
                                s=35
                                )
            title = preferred_titles[idx] if (preferred_titles is not None) else gene
            RolePlot.postprocess(ax, title=title, show_axis_names=False)
            RolePlot.mark_x(ax, 0)
            RolePlot.mark_y(ax, 0)
            RolePlot.plot_roles(ax)

            #- arrow
            data_ctr = data['ctr']
            data_sample = data['sample']
            Ps_ctr = [item[1] for item in data_ctr]
            Qs_ctr = [item[0] for item in data_ctr]

            Ps_sample = [item[1] for item in data_sample]
            Qs_sample = [item[0] for item in data_sample]
            x0, y0, x1, y1 = np.mean(Ps_ctr), np.mean(Qs_ctr), np.mean(Ps_sample), np.mean(
                Qs_sample)
            arrow_t = arrow_specs['arrow_type']
            rad = arrow_specs['rad']
            ax.annotate('', xy=(x1, y1), xytext=(x0, y0),
                        arrowprops={'arrowstyle': '->', 'lw': 1.5, 'connectionstyle': f'{arrow_t}, rad={rad}',
                                    'color': 'black', 'alpha': .7}
                        , horizontalalignment='center')
            ax.set_xmargin(.4)
            ax.set_ymargin(.4)

        tags = ['ctr','mg']
        handles = []
        for i, color in enumerate(RolePlot.ctr_sample_colors):
            handles.append(ax.scatter([],[], marker='o', label=tags[i], s=50, edgecolor='black', color=color, linewidth=.2))
        ax.legend(handles=handles,
                  loc = 'upper center',
                  bbox_to_anchor=(1.3,1),
            )
    def _test_significance(self, results):
        assert len(results) == len(self.target_genes)
        def do_it(ctr, sample):
            ctr_x = [point[0] for point in ctr]
            ctr_y = [point[1] for point in ctr]
            sample_x = [point[0] for point in sample]
            sample_y = [point[1] for point in sample]

            # Conduct two-sample t-tests for x and y coordinates
            _, p1 = ttest_ind(ctr_x, sample_x, equal_var=False)
            _, p2 = ttest_ind(ctr_y, sample_y, equal_var=False)
            if (p1<0.05) | (p2<0.05):
                return True
            else:
                False
        sig_flags = []
        for gene_i,_ in enumerate(self.target_genes):
            ctr = results[gene_i]['ctr']
            sample = results[gene_i]['sample']

            sig_flags.append(do_it(ctr,sample))
            # t_statistic, p_value = ttest_ind(r2_ctr, r2_sample)
            # p_values.append(p_value)
        return sig_flags

    def plot_results(self):
        ncols, nrows = len(self.target_genes), 2
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(2 * ncols, 2 * nrows))

        if self.results_mp is not None:
            p_values = self._test_significance(self.results_mp)
            signs = ['*' if (p ) else '' for p in p_values]
            # signs = ['*' if (p<0.05) else '' for p in p_values]
            preferred_titles = [f'{gene} {sign}' for gene, sign in zip(self.target_genes, signs)]
            self._plot(axes[0], self.results_mp, self.target_genes, preferred_titles=preferred_titles)
        axes[0][0].set_ylabel(f'Multiplicative noise ({int(100*self.std_mpnoise)} %)', fontweight='bold')
        if self.results_ad is not None:
            p_values = self._test_significance(self.results_ad)
            signs = ['*' if (p) else '' for p in p_values]
            # signs = ['*' if (p < 0.05) else '' for p in p_values]
            preferred_titles = [f'{gene} {sign}' for gene, sign in zip(self.target_genes, signs)]
            self._plot(axes[1], self.results_ad, self.target_genes, preferred_titles=preferred_titles)
        # axes[1][0].set_ylabel(f'Additive noise ({int(100*self.std_adnoise)} %)', fontweight='bold')
        return fig
    def __add_mpnoise(self, data, n_rep=100, sigma=1, std=.3) -> Tuple[np.array]: # deprecated
        """ Multiplicitive noise
             Creates n_relica noised data
        """
        data_std = np.std(data)
        applied_std = std*data_std
        noisy_data_stack = []
        for i in range(n_rep):
            rand_values = np.random.normal(loc=sigma, scale=applied_std, size=data.shape)
            noisy_data = rand_values*data
            noisy_data_stack.append(noisy_data)

        return noisy_data_stack

    def _add_mpnoise(self, data, n_rep=100, sigma=1, std=.3) -> Tuple[np.array]:
        """ Multiplicative noise
             Creates n_replica noised data
        """
        data_std = np.std(data)
        applied_std = std * data_std
        noisy_data_stack = []
        for i in range(n_rep):
            # frac_noise = .05
            frac_noise =  np.random.rand()
            noise_mask = np.random.choice([0, 1], size=data.shape, p=[1 - frac_noise, frac_noise])
            rand_values = np.random.normal(loc=sigma, scale=applied_std, size=data.shape)
            noisy_data = data * (1 + (rand_values-1) * noise_mask)
            noisy_data_stack.append(noisy_data)

        return noisy_data_stack

    def _add_adnoise(self, data, n_rep=100, sigma=1, std=.05) -> Tuple[pd.DataFrame]:
        """ Additive noise
                 Creates n_relica noised links
        """
        data_std = np.std(data)
        applied_std = data_std * std
        noisy_data_stack = []
        for i in range(n_rep):
            rand_values = np.random.normal(loc=sigma, scale=applied_std, size=data.shape)
            noisy_data = rand_values + data
            noisy_data_stack.append(noisy_data)

        return noisy_data_stack


    def _organized(self, batch_results_ctr, batch_results_sample):
        """
            Sorts batch VSA results based on target proteins
        """
        batch_results_studies = [batch_results_ctr, batch_results_sample]
        def extract_F(batch_results, prot, tag='Q'):
            #- extract results for prot
            batch_results_filtered = [results.query(f"Entry == '{prot}'") for results in batch_results]
            batch_results_ = np.asarray([results[tag].to_numpy()[0] for results in batch_results_filtered])
            return batch_results_

        results_stack = []
        for prot in self.target_genes:
            prot_results = {}
            for i_study, study in enumerate(self.studies):
                Qs = extract_F(batch_results_studies[i_study], prot=prot, tag='Q')
                Ps = extract_F(batch_results_studies[i_study], prot=prot, tag='P')
                prot_results[study] =  [(Q,P) for Q,P in zip(Qs,Ps)]
            results_stack.append(prot_results)
        return results_stack
        
        
        
        

