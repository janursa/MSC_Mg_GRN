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


def standardize_array(arr: List[float]) -> List[float]:
    """To set 0.75 quntile as the center and standardize to maximum value"""
    q75 = np.quantile(arr, 0.75)
    arr = arr - q75
    max_abs_val = np.max(np.abs(arr))
    arr = arr / max_abs_val
    return arr
def role_analysis(links, gene_names) -> pd.DataFrame:
    '''
        Vester's sensitivity analysis. 
            - active_sum: active sum. Sum along rows of the influence matrix and it indicates how much does a variable influence all the others.
            - passive_sum: passive sum. Its is the sum along columns of the influence matrix and it indicates how sensitive a variable is, how does it react to the influence of others
            - Q: active_sum/passive_sum -> how dominant
            - P: active_sum.passive_sum -> how participative a variable is
            - Active: +Q
            - Passive: -Q, -P
            - Critical: +Q, +P
            - Buffering: -Q, -P
    ------------------------------------------------------
    inputs: links (DataFrame)
    output: VSA (DataFrame) -> protname: active_sum, passive_sum, Role
    '''
    # active sum and passive sum
    # links = choose_top_quantile(links, 0.5)
    active_sum = np.asarray([sum(links.query(f"Regulator == '{gene}'")['Weight']) for gene in gene_names])
    passive_sum = np.asarray([sum(links.query(f"Target == '{gene}'")['Weight']) for gene in gene_names])

    # Quotient and Product
    quotient, product = active_sum/passive_sum, passive_sum*passive_sum

    # define the roles ['Buffering', 'Passive', 'Active', 'Critical'] -> [0, 1, 2, 3]
    quotient, product = standardize_array(quotient), standardize_array(product)
    quotient_flags = quotient>0
    product_flags = product>0
    roles = [2*quotient_flag+product_flag for quotient_flag, product_flag in zip(quotient_flags,product_flags)]
    vsa_results = pd.DataFrame(data={'Entry': gene_names, 'Q': quotient, 'P': product, 'Role':roles })
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
    q_distance = (sample['Q'].to_numpy(float) - control['Q'].to_numpy(float))**2
    p_distance = (sample['P'].to_numpy(float) - control['P'].to_numpy(float))**2
    distance = q_distance + p_distance
    df_distance = pd.DataFrame({'Entry': control['Entry'], 'Distance':distance})
    df_distance = df_distance.assign(P1=control['P'], Q1=control['Q'])
    df_distance = df_distance.assign(P2=sample['P'], Q2=sample['Q'])

    # - define 0.75 top percentile as the cut-off value to determine critical role
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

    @classmethod
    def plot_role_change(cls, df, ax, custom_annotation=None):
        '''
            Plots Q/P for role change betwen ctr and sample.
            It draws an arrow from ctr to sample.
        '''
        xlim = ax.get_xlim()[1] - ax.get_xlim()[0]
        ylim = ax.get_ylim()[1] - ax.get_ylim()[0]

        # - first ctr
        for i in [1,2]:
            ax.scatter(df[f'P{i}'], df[f'Q{i}'],
                            color=cls.ctr_sample_colors[i-1],
                            alpha=1,
                            linewidths=.2,
                            edgecolors='black',
                            s=100
                            )

        for prot, x0, y0, x1, y1 in zip(df['Entry'].values, df['P1'].values, df['Q1'].values, df['P2'].values,
                                        df['Q2'].values):

            if custom_annotation is not None:
                offset_x, offset_y, rad_arrow, offset_text_location = custom_annotation(prot)

            offset_x = offset_x or 0.015
            offset_y = offset_y or 0.1
            rad_arrow = rad_arrow or 0.3

            x, y = x0 + offset_x * xlim, y0 + offset_y * ylim
            if offset_text_location:
                ax.annotate(f'{prot}', xy=(x0, y0), xytext=(x, y), arrowprops={'arrowstyle': '->', 'lw':.6, 'color':'orange', 'linestyle':'--'}, fontsize=9)
            else:
                ax.annotate(f'{prot}', xy=(x, y), fontsize=9)

            ax.annotate(f'', xy=(x1, y1), xytext=(x0, y0),
                        arrowprops={'arrowstyle': '->', 'connectionstyle': f'arc3, rad={rad_arrow}'}
                        , horizontalalignment='center', fontsize=9)

        cls.mark_x(ax, 0)
        cls.mark_y(ax, 0)
        cls.plot_roles(ax)
        cls.postprocess(ax, title='Role change')

    @classmethod
    def plot_ctr_vs_sample(cls, ax, data, gene_names, custom_annotation):
        '''
            Plots Q/P for ctr and mg conditions in a 1*2 subplot
        '''
        serif_font()
        cls.scatter(ax, product=data['P'], quotient=data['Q'], roles=data['Role'])
        cls.mark_x(ax, 0)
        cls.mark_y(ax, 0)
        cls.postprocess(ax)
        cls.plot_roles(ax)

        xlim = ax.get_xlim()[1] - ax.get_xlim()[0]
        ylim = ax.get_ylim()[1] - ax.get_ylim()[0]

        # annotate target proteins names on the graph
        quotient = data['Q']
        product = data['P']
        roles = data['Role']
        for i, gene_name in enumerate(gene_names):
            if custom_annotation is not None:
                offset_x, offset_y, rad = custom_annotation(gene_name)

            offset_x = offset_x or 0.015
            offset_y = offset_y or 0.1
            rad = rad or 0.3

            if roles[i] != 3:  # just critical
                continue

            x = product[i]
            y = quotient[i]
            ax.annotate(gene_name, xy=(x, y), xytext=(x + offset_x * xlim, y + offset_y * ylim), fontsize=9,
                        arrowprops={'arrowstyle': '->', 'lw': 0.9, 'connectionstyle': f'arc3, rad={rad}',
                                    'color': 'black', 'alpha': .7}
                        , horizontalalignment='center')

    @staticmethod
    def postprocess(ax, title='', show_axis_names=True, xlim=[-1.4, 1.4], ylim=[[-1.4, 1.4]]):
        if show_axis_names:
            ax.set_xlabel(r'Product (AS $\times$ PS)', fontweight='bold', fontsize=9)
            ax.set_ylabel('Quotient (AS / PS)', fontweight='bold',  fontsize=9)
        ax.set_title(title, fontweight='bold')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks([])
        ax.set_yticks([])

    @staticmethod
    def mark_y(ax, y_t):
        """
        Mark threshold line for the sig
        """
        ax.axhline(y_t, color='grey', linestyle='-.')

    @classmethod
    def mark_x(cls, ax, x_t):
        """
        Mark threshold line for fold change
        """
        line_color = 'grey'
        ax.axvline(x_t, color=line_color, linestyle='--', linewidth=cls.linewidth)

    @classmethod
    def scatter(cls, ax, product, quotient, roles):

        ax.scatter(product, quotient,
                        color=[cls.roles_colors[i] for i in roles],
                        #                    alpha=[1 if flag else .5 for flag in flags]
                        #                    alpha = .7,
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
        ax.text(.01, .01, 'Buffering', ha='left', va='bottom', transform=ax.transAxes, fontweight='bold', fontsize=8)
        ax.text(.99, .01, 'Passive', ha='right', va='bottom', transform=ax.transAxes, fontweight='bold', fontsize=8)
        ax.text(.99, .99, 'Critical', ha='right', va='top', transform=ax.transAxes, fontweight='bold', fontsize=8)

