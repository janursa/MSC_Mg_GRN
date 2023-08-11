"""
    Sets of functions useful vester's sensitivity analysis
"""
import sys
import os
import numpy as np
import pandas as pd
from typing import List, Tuple


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from common_tools import serif_font
from common_tools.links import choose_top_quantile


def standardize_array(arr: List[float]) -> List[float]: #Deprecated
    """To set 0.75 quntile as the center and standardize to maximum value"""
    q75 = np.quantile(arr, 0.75)
    arr = arr - q75
    max_abs_val = np.max(np.abs(arr))
    arr = arr / max_abs_val
    return arr

def outlier_threshold(vector):
    return np.quantile(vector, 0.9)
def role_analysis(links, gene_names) -> pd.DataFrame:
    '''
        Vester's sensitivity analysis. 
            - active_sum: active sum. Sum along rows of the influence matrix and it indicates how much does a variable influence all the others.
            - passive_sum: passive sum. Its is the sum along columns of the influence matrix and it indicates how sensitive a variable is, how does it react to the influence of others
            - Active: +AS, -PS
            - Passive: -AS, -PS
            - Critical: +AS, +PS
            - neutral: -AS, -PS
    ------------------------------------------------------
    inputs: links (DataFrame)
    output: VSA (DataFrame) -> protname: active_sum, passive_sum, Role
    '''
    # active sum and passive sum
    links = choose_top_quantile(links, 0.75)
    active_sum = np.asarray([sum(links.query(f"Regulator == '{gene}'")['Weight']) for gene in gene_names])
    passive_sum = np.asarray([sum(links.query(f"Target == '{gene}'")['Weight']) for gene in gene_names])

    # define the roles ['Buffering', 'Passive', 'Active', 'Critical'] -> [0, 1, 2, 3]
    active_sum_threhsold = outlier_threshold(active_sum)
    passive_sum_threhsold = outlier_threshold(passive_sum)
    # active_sum, passive_sum = standardize_array(active_sum), standardize_array(passive_sum)
    roles = [2*as_flag+ps_flag for as_flag, ps_flag in zip(active_sum>active_sum_threhsold, passive_sum>passive_sum_threhsold)]
    vsa_results = pd.DataFrame(data={'Entry': gene_names, 'AS': active_sum, 'PS': passive_sum, 'Role':roles })
    return vsa_results
def determine_critical_role_change(df_ctr, df_sample, target_role:int=3):
    '''
        Finds those genes with a change in the role from ctr to sample.

    Inputs: target_role -> critical:3
    '''
    # - find those genes in both ctr and sample that have target role, i.e. critical
    indices_critical = list(map(lambda df: df.index[df['Role']==target_role].values, [df_ctr, df_sample]))
    indices_critical = [i for j in indices_critical for i in j] #flatten
    # - remove those that are present in both groups
    unique_elements, counts = np.unique(indices_critical, return_counts=True)
    indices_target = unique_elements[counts == 1].tolist()

    df_ctr_c = df_ctr.loc[indices_target,:].reset_index(drop=True)
    df_sample_c = df_sample.loc[indices_target,:].reset_index(drop=True)

    df_role_change = pd.DataFrame()
    df_role_change['Entry'] = df_ctr_c['Entry']
    df_role_change[['AS1','PS1']] = df_ctr_c[['AS','PS']]
    df_role_change[['AS2','PS2']] = df_sample_c[['AS','PS']]
    return df_role_change
def determine_top_role_change(control: pd.DataFrame, sample: pd.DataFrame, top_quantile: float = 0.95) -> pd.DataFrame:

    '''
    Find the distance in the role from control to sample and return the top_quantile of the distances.

    Args:
        control (pd.DataFrame): The control DataFrame with columns 'Entry', 'AS', and 'PS'.
        sample (pd.DataFrame): The sample DataFrame
        top_quantile (float): The quantile to use for selecting the top distances.

    Returns:
        pd.DataFrame: A DataFrame with columns 'Entry' and 'Distance', containing the entries with the top distances.

    '''
    # - calculate distance
    as_distance = (sample['AS'].to_numpy(float) - control['AS'].to_numpy(float))**2
    ps_distance = (sample['PS'].to_numpy(float) - control['PS'].to_numpy(float))**2
    distance = as_distance + ps_distance
    df_distance = pd.DataFrame({'Entry': control['Entry'], 'Distance':distance})
    df_distance = df_distance.assign(PS1=control['PS'], AS1=control['AS'])
    df_distance = df_distance.assign(PS2=sample['PS'], AS2=sample['AS'])

    # - define 0.75 top percentile as the cut-off value to determine critical role
    cut_off = np.quantile(df_distance['Distance'].values.tolist(), q=top_quantile)
    short_list = df_distance.loc[df_distance['Distance'] >= cut_off, :].reset_index(drop=True)
    return short_list
class RolePlot:
    """
    Plots VSA results
    """
    linewidth = 1
    markers = ['o', 'o', 'o', '*']
    roles_colors = ['lightgreen', 'lightblue', 'yellow', 'red']
    ctr_sample_colors = ['cyan', 'purple']
    sizes = [40 * i for i in [1, 1.5, 2, 5]]
    roles = ['Neutral', 'Passive', 'Active', 'Critical']
    def __repr__(self):
        return f'Plot tools for VSA roles: {self.roles}'
    @classmethod
    def create_role_legends(cls, ax):
        """
        Generates handles for role legends
        """
        handles = []
        for i, color in enumerate(cls.roles_colors):
            handles.append(ax.scatter([], [], marker='o', label=cls.roles[i], s=100,
                                                  edgecolor='black', color=color, linewidth=.5))
        return handles
    @staticmethod
    def annotate_gene_names(ax, name, x0, y0, offset_x, offset_y):
        """
            Annotate the gene name with an oragne arrow offsetting the location of the scatter point
        """
        xlim = ax.get_xlim()[1] - ax.get_xlim()[0]
        ylim = ax.get_ylim()[1] - ax.get_ylim()[0]
        offset_x = offset_x or 0.015
        offset_y = offset_y or 0.1
        x, y = x0 + offset_x * xlim, y0 + offset_y * ylim
        ax.annotate(f'{name}', xy=(x0, y0), xytext=(x, y),
                        arrowprops={'arrowstyle': '->', 'lw': .6, 'color': 'orange', 'linestyle': '--', 'alpha': .7}, fontsize=8)

    @classmethod
    def plot_role_change(cls, df, ax, custom_annotation=None, active_sum_range=None, passive_sum_range=None):
        '''
            Plots Q/P for role change betwen ctr and sample.
            It draws an arrow from ctr to sample.
        '''

        # - first ctr
        for i in [1,2]:
            ax.scatter(df[f'PS{i}'], df[f'AS{i}'],
                            color=cls.ctr_sample_colors[i-1],
                            alpha=.5,
                            linewidths=.2,
                            edgecolors='black',
                            s=100
                            )

        for gene, x0, y0, x1, y1 in zip(df['Entry'].values, df['PS1'].values, df['AS1'].values, df['PS2'].values,
                                        df['AS2'].values):
            if custom_annotation is not None:
                offset_x, offset_y, rad_arrow = custom_annotation(gene)
            rad_arrow = rad_arrow or 0.3
            # - draw the arrow to indicate gene name (orange)
            cls.annotate_gene_names(ax, gene, x0, y0, offset_x, offset_y)
            # - draw the role change arrow
            ax.annotate(f'', xy=(x1, y1), xytext=(x0, y0),
                        arrowprops={'arrowstyle': '->', 'connectionstyle': f'arc3, rad={rad_arrow}', 'alpha': .7}
                        , horizontalalignment='center', fontsize=9)

        # cls.mark_x(ax, 0)
        # cls.mark_y(ax, 0)
        cls.plot_roles(ax)
        cls.postprocess(ax, title='Role change', xlim=passive_sum_range, ylim=active_sum_range)

    @classmethod
    def plot_ctr_vs_sample(cls, ax, data, gene_names, custom_annotation, active_sum_range=None, passive_sum_range=None):
        '''
            Plots Q/P for ctr and mg conditions in a 1*2 subplot
        '''
        serif_font()
        active_sum = data['AS']
        passive_sum = data['PS']
        cls.scatter(ax, active_sum=active_sum, passive_sum=passive_sum, roles=data['Role'])
        cls.mark_x(ax, outlier_threshold(passive_sum))
        cls.mark_y(ax, outlier_threshold(active_sum))
        cls.postprocess(ax, xlim=passive_sum_range, ylim=active_sum_range)
        cls.plot_roles(ax)

        # annotate target proteins names on the graph
        roles = data['Role']
        for i, gene_name in enumerate(gene_names):
            if custom_annotation is not None:
                offset_x, offset_y = custom_annotation(gene_name)

            offset_x = offset_x or 0.015
            offset_y = offset_y or 0.1

            if roles[i] != 3:  # just critical
                continue

            x = passive_sum[i]
            y = active_sum[i]

            cls.annotate_gene_names(ax, gene_name, x, y, offset_x, offset_y)

    @staticmethod
    def postprocess(ax, title='', show_axis_names=True, xlim=None, ylim=None):
        if show_axis_names:
            ax.set_xlabel('Passive Sum', fontweight='bold', fontsize=9)
            ax.set_ylabel('Active Sum', fontweight='bold',  fontsize=9)
        ax.set_title(title, fontweight='bold')
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        ax.set_xticks([])
        ax.set_yticks([])
    @classmethod
    def mark_y(cls, ax, y_t):
        """
        Mark threshold line for the sig
        """
        ax.axhline(y_t, color='grey', linestyle='-.', linewidth=cls.linewidth)

    @classmethod
    def mark_x(cls, ax, x_t):
        """
        Mark threshold line for fold change
        """
        line_color = 'grey'
        ax.axvline(x_t, color=line_color, linestyle='--', linewidth=cls.linewidth)

    @classmethod
    def scatter(cls, ax, active_sum, passive_sum, roles):

        ax.scatter(passive_sum, active_sum,
                        color=[cls.roles_colors[i] for i in roles],
                        alpha = .7,
                        linewidths=.1,
                        edgecolors='black',
                        s=100
                        )

    @staticmethod
    def plot_roles(ax):
        '''
            Plots the roles on the corners of the scatter plot
        '''
        ax.text(.01, .99, 'Active', ha='left', va='top', transform=ax.transAxes, fontweight='bold', fontsize=8)
        ax.text(.01, .01, 'Neutral', ha='left', va='bottom', transform=ax.transAxes, fontweight='bold', fontsize=8)
        ax.text(.99, .01, 'Passive', ha='right', va='bottom', transform=ax.transAxes, fontweight='bold', fontsize=8)
        ax.text(.99, .99, 'Critical', ha='right', va='top', transform=ax.transAxes, fontweight='bold', fontsize=8)

