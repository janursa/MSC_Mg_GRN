# extract time series sig proteins: at least 1 overexpression of e = 0.5, and passed sig test
import importlib
importlib.reload(utils)
df_ov = utils.test_overexpression(df, **specs)
print('Overexpressed proteins: {}'.format(len(df_ov)))
# impute missing values: available techniques to try: interpolate, univariate and multivariate: https://scikit-learn.org/stable/modules/impute.html
df_interp = df_ov.interpolate() 
df_interp = utils.listwise_deletion(df_interp)
print('Inter proteins: {}'.format(len(df_interp)))
df_sig  = utils.test_sig(df_interp, **specs)
print('Sig proteins: {}'.format(len(df_sig)))
print(df_sig.head())

utils.plot_time_series(df, prots = ['P06865'], **specs)
utils.plot_time_series(df_interp, prots = ['P06865'], **specs)

# extract pairwise sig proteins
df_sig = utils.sig_df(df, **specs)
print('Sig proteins: {}'.format(len(df_sig)))
# impute missing values: available techniques to try: interpolate, univariate and multivariate: https://scikit-learn.org/stable/modules/impute.html
df_interp = df_sig.interpolate() 
df_interp = utils.listwise_deletion(df_interp)
print('Sig interpolated proteins: {}'.format(len(df_interp)))

# stats of sig proteins: different time points, unknown proteins

# number of sig proteins for different time points
sig_prot_names = {} 
for day in specs['time']:
    sig_prot_names[day] = utils.sig_test(df,day,**specs)[specs['p_ID']].tolist()
# plot number of sig proteins for measurement times 
fig, axes = plt.subplots(figsize = (15,4), nrows = 1, ncols = 3, tight_layout = True)
def plot_stats_time_1(ax):
    ax.bar(range(len(sig_prot_names)),[len(x) for x in sig_prot_names.values()])
    ax.set_ylabel('Number of proteins')
    ax.set_xlabel('Measurement day')
    ax.set_xticks(range(len(sig_prot_names)))
    _ = ax.set_xticklabels(sig_prot_names.keys())
    ax.set_title('Differentially expressed proteins')
plot_stats_time(axes[0])
#------
# sig incidence number: how many times a protein is detected sig across different time points
limit_n = 3 # if the incidences are below this, dont plot it
sig_prots_e = [x for xx in sig_prot_names.values() for x in xx] # duplicates included
sig_prots = list(set(sig_prots_e)) # duplicates excluded
print('Number of prots upregulated in at least one time point: ',len(sig_prots))
sig_prots_i = {i:sig_prots_e.count(i) for i in sig_prots} # incidences 
sig_prots_i_f = {}
for key,value in sig_prots_i.items():
    if value >= limit_n:
        sig_prots_i_f[key] = value
print('Number of prots upregulated in at least 3 time points: ',len(sig_prots_i_f))
def plot_stats_time_2(ax):
    ax.bar(range(len(sig_prots_i_f)),sig_prots_i_f.values())
    ax.set_ylabel('Number of incidences')
    ax.set_xlabel('Measurement day')
    ax.set_xticks(range(len(sig_prots_i_f)))
    ax.set_xticklabels(sig_prots_i_f.keys())
    ax.set_title('Multiple differentially expressed proteins')
plot_stats_time_2(axes[1])
#----
# sig incidence number for unknown params
unknown_sigs_prots = [x for x in sig_prots if 'p_' in x ]
unknown_sigs_prots_i = {i:sig_prots_e.count(i) for i in unknown_sigs_prots} # incidences 
print('Total number of unknown paramters: ',len([x for x in df[specs['p_ID']].tolist() if 'p_' in x ]))
print('Sig number of unknown paramters: ',len(unknown_sigs_prots))
def plot_stats_time_3(ax):
    ax.bar(range(len(unknown_sigs_prots_i)),unknown_sigs_prots_i.values())
    ax.set_ylabel('Number of incidences')
    ax.set_xticks(range(len(unknown_sigs_prots_i)))
    _ = ax.set_xticklabels(unknown_sigs_prots_i.keys())
    ax.set_title('Multiple differentially expressed proteins')
plot_stats_time_3(axes[2])
#----
# output the name of sig proteins 
def remove_unknown(sig_prot_names):
    sig_prots_time_f = {}
    for key,values in sig_prot_names.items():
        values_f = [x for x in values if x != np.nan]
        _str = ''
        for value in values_f:
            _str += value + '\n'
        sig_prots_time_f[key] = _str
    return sig_prots_time_f
sig_prots_time_f = remove_unknown(sig_prot_names)
with open('results/sig_prots.txt','w') as f:
    for key,value in sig_prots_time_f.items():
        f.write('\n \nday '+ str(key)+': \n'+ value)

# plot table of FC and pvalues: TODO:this needs work
from matplotlib.ticker import FormatStrFormatter
fig, ax = plt.subplots(1, 1, tight_layout=True, figsize=(4, 3))
fig.patch.set_visible(False)
ax.axis('off')
ax.axis('tight')
df_chosen = df_chosen.iloc[0:10,:]
df_chosen_r = df_chosen.round(4)
table= plt.table(cellText=df_chosen_r.values, colLabels=df_chosen_r.columns, loc='center')
from matplotlib.font_manager import FontProperties
for (row, col), cell in table.get_celld().items():
    if (row == 0):
        cell.set_text_props(fontproperties=FontProperties(weight='bold'))
fig.savefig(os.path.join(OUTPUT_DIR,'table.pdf'))

#------ plot for calibration pool -----# currently we only run once
def process(study):
    #- read best oo for all runs
    scores_pool, scores = utils.pool_scores(study, OUTPUT_DIR)
    bestparams_pool, bestparams, bestparams_poolofpool = utils.pool_bestparams(study, OUTPUT_DIR)
    #- convert oo to a df
    def create_oo_df(protnames, scores_pool, scores, bestparams_pool):
        scores_df = pd.DataFrame(data={'Protein':protnames, 'ScoreMean':scores, 'ScorePool': list(scores_pool), 'bestparams_pool':bestparams_pool})
        return scores_df
    oo_df = create_oo_df(protnames, scores_pool, scores, bestparams_pool)
    #- choose top prots based on oob scores: for visualization
    def choose_top(df, cut_off):
        df_sorted = df.sort_values('ScoreMean',ascending=False).reset_index(drop=True)
        return df_sorted.loc[df_sorted['ScoreMean']>cut_off,:]
#     oo_short = choose_top(oo_df, cut_off=.5)
    oo_short = oo_df.sort_values('ScoreMean',ascending=False).reset_index(drop=True)
    return oo_short, bestparams_poolofpool, scores
oo_short_ctr, bestparams_poolofpool_ctr, scores_ctr = process('ctr')
oo_short_sample, bestparams_poolofpool_sample, scores_sample = process('mg')

np.savetxt(os.path.join(OUTPUT_DIR,'calibration', 'scores_ctr.csv'), scores_ctr, delimiter=",")
np.savetxt(os.path.join(OUTPUT_DIR,'calibration', 'scores_mg.csv'), scores_sample, delimiter=",")#- plot scores
importlib.reload(utils)
#- box plot of pool scores
fig = utils.plot_bestscores_indivitualprots(data_ctr=oo_short_ctr['ScorePool'], 
                      data_sample=oo_short_sample['ScorePool'],
                      preferred_names_ctr= oo_short_ctr['Protein'],
                      preferred_names_sample= oo_short_sample['Protein'],
                      xlabel='Proteins')
fig.savefig(os.path.join(OUTPUT_DIR,'postprocess','bestscores.png'), dpi=300, transparent=False, facecolor='white')

#- box plot for seperate params
fig = utils.plot_bestparams_individualprots(oo_short_ctr['bestparams_pool'].values, 
                            oo_short_sample['bestparams_pool'].values,
                            priors=param_grid, 
                            preferred_names_ctr=oo_short_ctr['Protein'],
                            preferred_names_sample=oo_short_sample['Protein'],
                           
                           )
#- box plot for poolofpool
importlib.reload(utils)
fig = utils.plot_bestparams_poolofpool(bestparams_poolofpool_ctr,bestparams_poolofpool_sample, priors=param_grid)
fig.savefig(os.path.join(OUTPUT_DIR,'postprocess','bestparams_pool.png'), dpi=300, transparent=False, facecolor='white')


# volcanic plot: -log10(p values) vs log2(foldchange)
class Plot:
    linewidth = 2
    max_diff = np.max(df_diff.iloc[0:,1:].values)
    fold_change_line = np.log2(fold_change_t)
    @staticmethod
    def mark_p_value(ax,sig_cut_off):
        """
        Markline for the sig
        """
        line_color = 'grey'
        dx = ax.get_xlim()
        ax.axhline(sig_cut_off,color='grey', linestyle='-.')
    @staticmethod
    def mark_FC(ax,FC_cut_off):
        """
        Markline for fold change
        """
        line_color = 'grey'
        FC_cut_offs = [-FC_cut_off,FC_cut_off]
        for FC_cut_off in FC_cut_offs:
            ax.axvline(FC_cut_off,color=line_color, linestyle='--',linewidth=Plot.linewidth)
    @staticmethod
    def scatter(ax,log2FC,negLog10pv,flags):
        ax.scatter(log2FC,negLog10pv,
                   color=['red' if flag else 'black' for flag in flags],
#                    decay_coeff=[1 if flag else .5 for flag in flags]
                   decay_coeff = .7,
                   linewidths=.5,
                   edgecolors='cyan'
                  )
    def print_names(ax, log2FC, flag):
        fontsize = 7
        log2FC_s = log2FC[flag].values
        gene_names_s = df_flags.loc[flag, 'Protein']
        sign_flag = log2FC_s>=0
        gene_names_right = gene_names_s[sign_flag]
        gene_names_negative = gene_names_s[~sign_flag]
        for i, name in enumerate(gene_names_right):
            ax.annotate(name,(1.37,3.4-.3*i),fontsize=fontsize)
        for i, name in enumerate(gene_names_negative):
            ax.annotate(name,(-1.75,3.4-.3*i),fontsize=fontsize)
#         ax.annotate(text_right, (1.1,1),fontsize=fontsize)
    @staticmethod
    def postprocess(ax,time):
        ax.set_xlabel('Log₂ (fold-change)')
        ax.set_ylabel('- Log₁₀ (p-value)')
        ax.set_xlim([-1.3*Plot.max_diff, 1.3*Plot.max_diff])
        ax.set_title(f'Day {time}')
        ax.set_xticks([-Plot.max_diff,-Plot.fold_change_line,int(0),Plot.fold_change_line, Plot.max_diff])
        ax.set_xticklabels([-round(Plot.max_diff,2),-round(Plot.fold_change_line,2),int(0),round(Plot.fold_change_line,2),round(Plot.max_diff,2)])
def plot_multiple(axes, days, cols):
    for i, time in enumerate(days):
        ax = axes[int(i/cols)][i%cols]
        log2FC = df_diff[time]
        negLog10pv = -np.log10(df_sig['P.Value'].values)
        DE_flags_time = (df_flags['Sig'] * df_flags[f'FC_{time}']).values 
        Plot.scatter(ax,log2FC,negLog10pv, DE_flags_time)
        Plot.mark_p_value(ax,-np.log10(sig_t))
        Plot.mark_FC(ax,np.log2(fold_change_t))
        Plot.postprocess(ax,time)
        Plot.print_names(ax,log2FC,DE_flags_time)
def plot_selected():
    days = [1,3,10,14]
    rows = 2
    cols = 2
    fig, axes = plt.subplots(rows, cols, tight_layout=True, figsize=(cols*4, rows*3))
    plot_multiple(axes, days, cols)
    fig.savefig(os.path.join(OUTPUT_DIR,'statAnalysis','volcanic_selected.pdf'))
    fig.savefig(os.path.join(OUTPUT_DIR,'statAnalysis','volcanic_selected.png'),facecolor='white',dpi=300)

def plot_all():
    rows = 6
    cols = 2
    fig, axes = plt.subplots(rows, cols, tight_layout=True, figsize=(cols*4, rows*3))
    plot_multiple(axes, specs['time'], cols)
    fig.savefig(os.path.join(OUTPUT_DIR,'statAnalysis','volcanic_all.pdf'))
    fig.savefig(os.path.join(OUTPUT_DIR,'statAnalysis','volcanic_all.png'),facecolor='white',dpi=300)
# plot_all()
plot_selected()

# DE Based on both sig and fold change
# calculate fold change: log(mg)-log(ctr)
df_diff = pd.DataFrame() #columns: ['Protein','1','2', etc (time points)]
df_diff['Protein'] = df_imput['Protein']
for t in specs['time']:
    ctr = f'ctr_{t}'
    mg = f'mg_{t}'
    df_diff[t] = (df_imput.loc[:,mg]-df_imput.loc[:,ctr]).values
# select 
fold_change_t = 1.5
for t in specs['time']:
    tag = f'FC_{t}'
    df_flags[tag] = abs(df_diff[t]) > np.log2(fold_change_t)
df_flags['DE'] = df_flags['Sig'] * np.any(df_flags.iloc[0:,2:].values,axis=1)
print(list(df_flags['DE'].values).count(True))

#---- to format data for r plot 
def reireive_string_data_F(enrich_type):
    '''
        Reads string EA data such as enrichment.Function
    '''
    FILE =  os.path.join(MAIN_DIR, 'results/enrichment_analysis/enrichment.all.tsv')
    data = pd.read_csv(FILE,sep='\t',index_col=False)
    # keep the specified enrichment type
    filter_ = data.loc[:,'#category'] == enrich_type
    data = data.loc[filter_,:]
    data.drop('#category',axis=1, inplace=True)
    # remove some columns
    data.drop('matching proteins in your network (IDs)',axis=1,inplace=True) 
    data.drop('background gene count',axis=1,inplace=True) 
    data.columns = ['ID','Description', 'Count', 'Strength','FDR','genenames']
    # set gene ratio is in the format of x/total
    GeneRatio = []
    for item in data['Count']:
        GeneRatio.append(str(item)+'/'+str(len(DE_protnames)))
    data.loc[:,'GeneRatio'] = GeneRatio
    # not useful but required as input to the graph
    data['pvalue'] = data['FDR'] 
    data['qvalue'] = data['pvalue'] 
    data['BgRatio'] = data['GeneRatio'] #TODO: fix this
    data.rename(columns={'genenames':'geneID'}, inplace=True)
    cut_t = 40
    short_descript = []
    for item in data.loc[:,'Description']:
        if len(item)>cut_t:
            aa =  item[0:cut_t] + '...'
        else:
            aa = item
        short_descript.append(aa)
    
        
    data.loc[:,'Description'] = short_descript
    data.sort_values('FDR', inplace=True)

#     if enrich_type == 'GO Process': # select top n ones
#         n = len(data)
#     elif enrich_type == 'GO Component':
#         n = 10
#     elif enrich_type == 'UniProt Keywords':
#         n = 10
#     elif enrich_type == 'UniProt Keywords':
#         n = 10
#     elif enrich_type == 'Reactome':
#         n = 10
#     else:
    n = len(data)
    data.sort_values('pvalue', ascending=True, inplace=True)
    data = data.iloc[0:n,:]
        
    return data

def plot_bestparams_poolofpool(data_ctr, data_sample, priors):
    '''
        Plots box plot for all params in one window. 
        For each param, the pool of best params obtained for all genes are pooled together.
        This plot shows how variabled are the inferred params across all genes. 
    '''
    fig, axes = plt.subplots(1, 2, tight_layout=True, figsize=(7,4))
    def normalize(xx, priors):
        xx = {key: np.array(list(set(values))) for key, values in xx.items()}
        return {key: (values - min(priors[key])) / (max(priors[key]) - min(priors[key])) for key, values in
                    xx.items()}
    datas = [data_ctr, data_sample]
    titles = ['ctr','mg']
    for i, data in enumerate(datas):
        data = normalize(data, priors)
        axes[i].boxplot(data.values())
        axes[i].set_ylabel('Normalized quantity')
        axes[i].set_yticks([])
        axes[i].set_yticklabels([])
        axes[i].set_xticklabels(labels=data.keys(), rotation=45)
        axes[i].set_title(titles[i])
    return fig
def plot_bestparams_individualprots(data_ctr,data_sample, priors, preferred_names_ctr, preferred_names_sample):
    """
        Plots boxplot for indivual param in seperate window. In each window, the variation in 
        best params are given for each gene seperately, obtained during different runs of tunning.
    """
    nrows = len(priors)
    ncols = 2
    fig, axes = plt.subplots(nrows, ncols, gridspec_kw={'width_ratios': [3, 2]}, tight_layout=True, figsize=(4*ncols,nrows*3))
    datas = [data_ctr,data_sample]
    preferred_names_s = [preferred_names_ctr, preferred_names_sample]
    for i, key in enumerate(priors.keys()):
        for j,data in enumerate(datas):
            pool = []
            for item in data:
                pool.append(item[key])
            axes[i][j].boxplot(pool)
            axes[i][j].set_title(key)
            axes[i][j].set_ylabel('Quantity')
            axes[i][j].set_xticklabels(preferred_names_s[j], rotation=45)
    return fig

def plot_bestscores_indivitualprots(data_ctr, data_sample, preferred_names_ctr,preferred_names_sample, xlabel=''):
    """plots scores as a box plot for a set"""
    fig, axes = plt.subplots(1, 2, tight_layout=True, figsize=(12,4), gridspec_kw={'width_ratios': [2, 2]})
    data_s = [data_ctr, data_sample]
    titles = ['ctr','mg']
    preferred_names_s = [preferred_names_ctr, preferred_names_sample]
    for i, data in enumerate(data_s):
        axes[i].boxplot(data)
        axes[i].set_ylabel('Testing score')
        axes[i].set_xlabel(xlabel)
        axes[i].set_title(titles[i])
        # axes[i].set_xticks(range(1, len(tags) + 1))
        axes[i].set_xticklabels(preferred_names_s[i],rotation=45)
    return fig
#- show all installed fonts
import matplotlib.font_manager
from IPython.core.display import HTML

def make_html(fontname):
    return "<p>{font}: <span style='font-family:{font}; font-size: 24px;'>{font}</p>".format(font=fontname)

code = "\n".join([make_html(font) for font in sorted(set([f.name for f in matplotlib.font_manager.fontManager.ttflist]))])

# HTML("<div style='column-count: 2;'>{}</div>".format(code))