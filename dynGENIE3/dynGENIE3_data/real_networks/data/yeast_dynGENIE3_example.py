import _pickle

with open('yeast_data.pkl','rb') as f:
    (TS_data, time_points, genes, TFs, alphas) = _pickle.load(f)

# TS_data: time series data
# time_points: corresponding time points
# genes: gene names
# TFs: transcription factors
# alphas: gene degradation rates

(VIM, alphas, prediction_score, stability_score, treeEstimators) = dynGENIE3(TS_data,time_points,alpha=alphas,gene_names=genes,regulators=TFs)
