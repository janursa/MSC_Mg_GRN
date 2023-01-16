"""
    Sets of functions useful vester's sensitivity analysis
"""
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
import utils
from ._imports import *
linewidth = 2
markers = ['o','o','o','*']
colors = ['lightgreen','lightblue','yellow','r']
sizes = [40*i for i in [1,1.5,2,5]]
roles = ['Buffering','Passive','Active','Critical']
@staticmethod
def thresholds(AS, PS) -> typing.Tuple[float,float]:
    '''
        Calculate the threshold to mark fold change and pvalue 
    '''
    # def pool(df_1, df_2, key):
    #     return list(df_1[key])+list(df_2[key])+
    # PS = pool(df_ctr, df_sample, 'PS')
    # AS = pool(df_ctr, df_sample, 'AS')
    # PS_t, AS_t = np.mean(PS), np.mean(AS)
    PS_t, AS_t = np.quantile(PS,.75), np.quantile(AS,.75)
    # AS_t = 0
    # PS_t = 0
    return AS_t, PS_t

@staticmethod
def analyse(links, protnames) -> pd.DataFrame:
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
    

    AS = [sum(links.loc[links['Regulator']==gene,:]['Weight']) for gene in protnames]
    PS = [sum(links.loc[links['Target']==gene,:]['Weight']) for gene in protnames] 
    #- normalize links
    AS, PS = AS/np.std(AS), PS/np.std(PS)
    vsa_results = pd.DataFrame(data={'Entry':protnames,'AS':AS, 'PS':PS})
    #- define the roles: ['Buffering', 'Passive', 'Active', 'Critical'] -> [0, 1, 2, 3] 
    AS_t, PS_t = VSA.thresholds(vsa_results['AS'].values, vsa_results['PS'].values)

    vsa_results.loc[:,'AS'] = vsa_results['AS'] - AS_t 
    vsa_results.loc[:,'PS'] = vsa_results['PS'] - PS_t 
    AS_flags = vsa_results['AS'].values>0
    PS_flags = vsa_results['PS'].values>0
    roles = [2*AS_flag+PS_flag for AS_flag, PS_flag in zip(AS_flags,PS_flags)]
    vsa_results.loc[:,'Role'] = roles
    return vsa_results
@staticmethod
def mark_y(ax, y_t):
    """
    Mark threshold line for the sig
    """
    line_color = 'grey'
    dx = ax.get_xlim()
    ax.axhline(y_t,color='grey', linestyle='-.')
@staticmethod
def mark_x(ax, x_t):
    """
    Mark threshold line for fold change
    """
    line_color = 'grey'
    ax.axvline(x_t,color=line_color, linestyle='--',linewidth=VSA.linewidth)
@staticmethod
def scatter(ax,PS, AS, roles):
    
    sc = ax.scatter(PS,AS,
               color=[VSA.colors[i] for i in roles],
#                    alpha=[1 if flag else .5 for flag in flags]
#                    alpha = .7,
               linewidths=.1,
               edgecolors='black',
               s=50
              # s = [VSA.sizes[i] for i in roles],
              )
    # markers = [VSA.markers[i] for i in roles]
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
def plot_prot_names(ax, AS, PS, gene_names, roles):
    '''
        Prints names of critical genes on the scatter plot
    '''
    xlim = ax.get_xlim()[1]-ax.get_xlim()[0]
    ylim = ax.get_ylim()[1]-ax.get_ylim()[0]
    
    arrow_specs = {'arrow_type':'arc3', 'rad':.2}
    for i, gene_name in enumerate(gene_names):
        offset_x = 0.2
        offset_y = 0.2
        rad = .2
        role = roles[i]
        if role != 3: #just critical 
            continue
        x = PS[i]
        y = AS[i]
        if gene_name == 'Q02218' or gene_name=='P26447':
            offset_y = -.05
            offset_x = .3
            rad=.2
        if gene_name == 'P62854':
            offset_y = -.25
            offset_x = .15
            rad=.2
        if gene_name == 'Q07065':
            offset_y = .2
            offset_x = -.25
            rad=-.2
        if gene_name == 'Q9UJZ1':
            offset_y = .15
            offset_x = 0.1
            rad=-.2
        if gene_name == 'P21281':
            offset_y = -.05
            offset_x = 0.22
            rad=.2
        # if gene_name == 'P21281', 'P26447', 'P62854','Q02218','Q07065','Q9UJZ1'


        ax.annotate(gene_name,xy=(x,y), xytext=(x+offset_x*xlim,y+offset_y*ylim), fontsize=7,
         arrowprops={'arrowstyle':'->','lw':0.9, 'connectionstyle':f'arc3, rad={rad}','color':'black', 'alpha':.7}
         ,horizontalalignment='center')
        # ax.annotate(gene_name, (x+offset_x,y+offset_y), fontsize=fontsize)
    

    
    # fontsize = 7
    
    

@staticmethod
def postprocess(ax, title, show_axis_names=True):
    if show_axis_names:
        ax.set_xlabel('PS')
        ax.set_ylabel('AS')

    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    xlen, ylen = xlim[1]-xlim[0], ylim[1]-ylim[0]
    ax.set_title(title)
    ax.margins(.3)
    ax.set_xticks([])
    ax.set_yticks([])
@staticmethod
def role_change(df_ctr, df_sample, target_roles=[3]):    
    '''
        Finds those genes with a change in the role from ctr to sample.

    Inputs: target_role -> one of 0,1,2,3
    '''
    flags = []
    for target_role in target_roles:
        flags.append(list(df_ctr['Role'] == target_role))
        flags.append(list(df_sample['Role'] == target_role))
    flags = np.array(flags).T  
    flags = [np.any(row) for row in flags]
    
    df_ctr_c = df_ctr.loc[flags,:].reset_index(drop=True)
    df_sample_c = df_sample.loc[flags,:].reset_index(drop=True)

    df_role_changed = pd.DataFrame()
    df_role_changed['Entry'] = df_ctr_c['Entry']
    df_role_changed[['AS_1','PS_1','Role_1']]=df_ctr_c[['AS','PS','Role']]
    df_role_changed[['AS_2','PS_2','Role_2']]=df_sample_c[['AS','PS','Role']]
    return df_role_changed
def plot_ctr_vs_sample(df_ctr, df_sample, preferred_names, OUTPUT_DIR=''):
    '''
        Plots AS/PS for ctr and mg conditions in a 1*2 subplot
    '''
    utils.comic_font()
    DIR = os.path.join(OUTPUT_DIR, 'VSA')
    rows = 1
    cols = 2
    fig, axes = plt.subplots(rows, cols, tight_layout=True, figsize=(cols*3.5, rows*3))
    dfs = [df_ctr, df_sample]
    titles = ['Ctr','Mg']
    for j in range(cols):
        df = dfs[j]
        if df is None:
            continue
        ax = axes[j]
        VSA.scatter(ax, PS=df['PS'], AS=df['AS'], roles=df['Role'])
        # AS_t, PS_t = VSA.thresholds(df['AS'], df['PS'])
        VSA.mark_x(ax, 0)
        VSA.mark_y(ax, 0)
        VSA.plot_prot_names(ax=ax, AS=df['AS'], PS=df['PS'], gene_names=preferred_names, roles=df['Role'])
        VSA.postprocess(ax, title=titles[j])
        VSA.plot_roles(ax)

    handles = []
    for i, color in enumerate(VSA.colors):
        handles.append(ax.scatter([],[],marker='o', label=VSA.roles[i],
         edgecolor='black', color=color, linewidth=.2))
    ll = plt.legend(handles=handles, 
        bbox_to_anchor=(1,1), prop={'size': 8}
        # title='Enriched Term'
        )
    # ax.add_artist(ll)

    fig.savefig(os.path.join(DIR, 'roles_ctr_mg.pdf'),bbox_extra_artists=(ll,), bbox_inches='tight')
    fig.savefig(os.path.join(DIR, 'roles_ctr_mg.png'), dpi=300, transparent=True, bbox_extra_artists=(ll,), bbox_inches='tight')
def plot_SA(VSAs, name, OUTPUT_DIR):
    '''
        Plots AS/PS for ctr and mg conditions, for n different runs of sensitivity analysis, one window for each prot
    '''
    utils.comic_font() 
    matplotlib.rcParams.update({'font.size': 8})
    arrow_specs = {'arrow_type':'arc3', 'rad':.2}
    ncols = 3
    nrows = 2
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, layout='tight', figsize=(1.6 * ncols, 1.6 * nrows))
    
    colors = ['lightgreen','purple']
    for idx,(prot, data) in enumerate(VSAs.items()):
        i = int(idx/(ncols))
        j = idx%ncols
        # print(i,j)
        ax = axes[i][j]

        # ax = axes[idx]
        #- plot role change 
        
        for idxx, study in enumerate(data.keys()):
            sc = ax.scatter(data[study]['PS'], data[study]['AS'],
               color=colors[idxx],
               alpha = .7,
               linewidths=.2,
               edgecolors='black',
               s = 25
              # s = [VSA.sizes[i] for i in roles],
              )
            if np.mean(data[study]['PS'])>0 and np.mean(data[study]['AS'])>0:
                x = np.mean(data[study]['PS'])
                y = np.mean(data[study]['AS'])+.8

                ax.annotate(r'$*$', xy=(x,y), fontsize = 12, color='red')

        VSA.postprocess(ax, title=prot, show_axis_names=False)
        VSA.mark_x(ax, 0)
        VSA.mark_y(ax, 0)
        VSA.plot_roles(ax)
        ax.set_xlim([-3,3])
        ax.set_ylim([-3,3])
        #- arrow
        fontsize = 8
        x0, y0, x1, y1  = np.mean(data['ctr']['PS']), np.mean(data['ctr']['AS']), np.mean(data['mg']['PS']), np.mean(data['mg']['AS'])
        arrow_t = arrow_specs['arrow_type']
        rad = arrow_specs['rad']
        ax.annotate('',xy=(x1,y1), xytext=(x0,y0),
         arrowprops={'arrowstyle':'->','lw':1.5, 'connectionstyle':f'{arrow_t}, rad={rad}','color':'black', 'alpha':.7}
         ,horizontalalignment='center')

    # tags = ['Ctr','Mg']
    # handles = []
    # for i, color in enumerate(colors):
    #     handles.append(ax.scatter([],[],marker='o', label=tags[i], edgecolor='black', color=color, linewidth=.2))
    # ll = plt.legend(handles=handles, 
    #     bbox_to_anchor=(1.3,3.4), 
    #     )
    # fig.savefig(os.path.join(OUTPUT_DIR, f'VSA/{name}.png'), dpi=300, transparent=True, bbox_extra_artists=(ll,), bbox_inches='tight')
    # fig.savefig(os.path.join(OUTPUT_DIR, f'VSA/{name}.pdf'),bbox_extra_artists=(ll,), bbox_inches='tight')
    fig.savefig(os.path.join(OUTPUT_DIR, f'VSA/{name}.pdf'), bbox_inches='tight')
    fig.savefig(os.path.join(OUTPUT_DIR, f'VSA/{name}.png'), bbox_inches='tight',dpi=300)

def plot_role_change(df, OUTPUT_DIR=''):
    '''
        Plots AS/PS for role change betwen ctr and sample
    '''
    DIR = os.path.join(OUTPUT_DIR, 'VSA')
    #- define default specs for arrows
    arrow_specs = {}
    for prot in df['Entry'].values:
        arrow_specs[prot] = {'arrow_type':'arc3', 'rad':.3}
    #- manual changes to arrow specs
    # arrow_specs['P46379']['rad'] = -.1
    # gene names shown on ctr and not mg
    # an arrow from ctr to mg only for any change from or to critical role
    # names of ctr or mg on the symbols
    rows = 1
    cols = 1
    fig, axes = plt.subplots(rows, cols, tight_layout=True, figsize=(cols*4, rows*4))
    ax = axes 
    #- first ctr
    VSA.scatter(ax, PS=df['PS_1'], AS=df['AS_1'], roles=df['Role_1'])
    #- second sample
    VSA.scatter(ax, PS=df['PS_2'], AS=df['AS_2'], roles=df['Role_2'])
    #- arrow from ctr to sample
    fontsize = 8
    for prot, x0, y0, x1, y1 in zip(df['Entry'].values, df['PS_1'].values, df['AS_1'].values, df['PS_2'].values, df['AS_2'].values): 
        print(prot)
        x, y = x0-.15, y0-.012
        # if prot == 'Q07954' or prot=='Q14444':
        #     x, y = x0-.14, y0+.01
        # if prot == 'P12956':
        #     x, y = x1-.1, y1-.012
        ax.annotate(f'{prot}', xy=(x,y), fontsize = 7)
        arrow_t = arrow_specs[prot]['arrow_type']
        rad = arrow_specs[prot]['rad']
        ax.annotate('',xy=(x1,y1), xytext=(x0,y0),
         arrowprops={'arrowstyle':'->', 'connectionstyle':f'{arrow_t}, rad={rad}'}
         ,horizontalalignment='center')
    
    VSA.mark_x(ax, 0)
    VSA.mark_y(ax, 0)
    VSA.plot_roles(ax)
    VSA.postprocess(ax, title='')

    fig.savefig(os.path.join(DIR, 'roles_ctr_to_mg.png'), dpi=300, transparent=True)

