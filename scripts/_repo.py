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