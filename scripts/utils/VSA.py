"""
    Sets of functions useful vester's sensitivity analysis
"""
import sys
import os
import typing
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import math

import scipy

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from scripts.utils import serif_font, comic_font

def Q_P_thresholds(AS, PS) -> typing.Tuple[float,float]:
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
def VestersSA(links, protnames) -> pd.DataFrame:
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
    AS = np.asarray([sum(links.query(f"Regulator == '{gene}'")['Weight']) for gene in protnames])
    PS = np.asarray([sum(links.query(f"Target == '{gene}'")['Weight']) for gene in protnames])
    #- Quotient and Product
    Q, P = AS/PS, PS*PS
    # Q, P = AS , PS

    #- normalize links
    # AS, PS = AS/np.std(AS), PS/np.std(PS)

    #- define the roles: ['Buffering', 'Passive', 'Active', 'Critical'] -> [0, 1, 2, 3] 
    Q_t, P_t = Q_P_thresholds(Q, P)
    Q = Q - Q_t
    P = P - P_t
    Q_flags = Q>0
    P_flags = P>0
    roles = [2*Q_flag+P_flag for Q_flag, P_flag in zip(Q_flags,P_flags)]
    vsa_results = pd.DataFrame(data={'Entry': protnames, 'Q': Q, 'P': P, 'Role':roles })
    return vsa_results
def role_change(df_ctr, df_sample, target_role:int=3):
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

    df_role_changed = pd.DataFrame()
    df_role_changed['Entry'] = df_ctr_c['Entry']
    df_role_changed[['Q_1','P_1','Role_1']]=df_ctr_c[['Q','P','Role']]
    df_role_changed[['Q_2','P_2','Role_2']]=df_sample_c[['Q','P','Role']]
    return df_role_changed
class VSA_plot:
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
        for i, color in enumerate(VSA_plot.roles_colors):
            handles.append(ax.scatter([], [], marker='o', label=VSA_plot.roles[i], s=100,
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
            ax.scatter(df[f'P_{i}'], df[f'Q_{i}'],
                            color=VSA_plot.ctr_sample_colors[i-1],
                            alpha=1,
                            linewidths=.2,
                            edgecolors='black',
                            s=100
                            )

        for prot, x0, y0, x1, y1 in zip(df['Entry'].values, df['P_1'].values, df['Q_1'].values, df['P_2'].values,
                                        df['Q_2'].values):
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

        VSA_plot.mark_x(ax, 0)
        VSA_plot.mark_y(ax, 0)
        VSA_plot.plot_roles(ax)
        VSA_plot.postprocess(ax, title='Role change')


    def plot_noise_analysis(axes, oo_prots, study_1='ctr', study_2='sample', show_title=False):
        '''
            Plots Q/P for ctr and mg conditions, for n different runs of sensitivity analysis, one window for each prot

            oo_prots: Q/P data for ctr and sample, for each protein
        '''
        comic_font()
        # matplotlib.rcParams.update({'font.size': 10})
        arrow_specs = {'arrow_type': 'arc3', 'rad': .2}

        for idx, (prot, data) in enumerate(oo_prots.items()): #- plot for each prot on a seperate window
            ax = axes[idx]
            xlim, ylim = ax.get_xlim(), ax.get_ylim()
            xlen, ylen = xlim[1] - xlim[0], ylim[1] - ylim[0]

            sig_flags = [False, False] # both need to be True to annotate * on the change
            roles = [] #to check if the role changes; its -1, 0, 1
            xy = []
            for idxx, study in enumerate(data.keys()):
                sc = ax.scatter(data[study]['P'], data[study]['Q'],
                                color=VSA_plot.ctr_sample_colors[idxx],
                                alpha=.6,
                                linewidths=.2,
                                edgecolors='black',
                                s=35
                                )

                P_vector = data[study]['P']
                Q_vector = data[study]['Q']
                xy.append([np.mean(P_vector),np.mean(Q_vector)])
                if (np.mean(P_vector)>0) & (np.mean(Q_vector)>0):
                    roles.append(1)
                else:
                    roles.append(0)
                _, p1 = scipy.stats.ttest_1samp(P_vector, 0)
                _, p2 = scipy.stats.ttest_1samp(Q_vector, 0)

                if (p1<0.05) & (p2<0.05):
                    sig_flags[idxx] = True

            if (roles[0] != roles[1]) & (sig_flags[0] & sig_flags[1]):
                ind = roles.index(1)
                x = xy[ind][0]
                y = xy[ind][1] + .1*ylen

                ax.annotate(r'$*$', xy=(x, y), fontsize=15, color='red', fontweight='bold')
            title = prot if show_title else ''
            VSA_plot.postprocess(ax, title=title, show_axis_names=False)
            VSA_plot.mark_x(ax, 0)
            VSA_plot.mark_y(ax, 0)
            VSA_plot.plot_roles(ax)

            # - arrow
            x0, y0, x1, y1 = np.mean(data[study_1]['P']), np.mean(data[study_1]['Q']), np.mean(data[study_2]['P']), np.mean(
                data[study_2]['Q'])
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
        for i, color in enumerate(VSA_plot.ctr_sample_colors):
            handles.append(ax.scatter([],[], marker='o', label=tags[i], s=50, edgecolor='black', color=color, linewidth=.2))
        ax.legend(handles=handles,
                  loc = 'upper center',
                  bbox_to_anchor=(1.3,1),
            )


    def plot_ctr_vs_sample(axes, data_stack, preferred_names, gene_names):
        '''
            Plots Q/P for ctr and mg conditions in a 1*2 subplot
        '''
        comic_font()

        for j, df in enumerate(data_stack):
            ax = axes[j]
            VSA_plot.scatter(ax, P=df['P'], Q=df['Q'], roles=df['Role'])
            # Q_t, P_t = VSA_plot.thresholds(df['Q'], df['P'])
            VSA_plot.mark_x(ax, 0)
            VSA_plot.mark_y(ax, 0)
            VSA_plot.plot_prot_names(ax=ax, Q=df['Q'], P=df['P'], gene_names=gene_names, roles=df['Role'])
            VSA_plot.postprocess(ax, title=preferred_names[j])
            VSA_plot.plot_roles(ax)

        # handles = []
        # for i, color in enumerate(VSA_plot.roles_colors):
        #     handles.append(plt.scatter([], [], marker='o', label=VSA_plot.roles[i],
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
        ax.axvline(x_t, color=line_color, linestyle='--', linewidth=VSA_plot.linewidth)

    @staticmethod
    def scatter(ax, P, Q, roles):

        sc = ax.scatter(P, Q,
                        color=[VSA_plot.roles_colors[i] for i in roles],
                        #                    alpha=[1 if flag else .5 for flag in flags]
                        #                    alpha = .7,
                        linewidths=.1,
                        edgecolors='black',
                        s=100
                        )
        # markers = [VSA_plot.markers[i] for i in roles]
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



