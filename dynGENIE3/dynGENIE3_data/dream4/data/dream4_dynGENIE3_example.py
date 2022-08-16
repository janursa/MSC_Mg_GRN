import _pickle

with open('size10_1_data.pkl', 'rb') as f:
    (TS_data, time_points, SS_data) = _pickle.load(f)

# TS_data: time series data
# time_points: corresponding time points
# SS_data: steady-state data

(VIM, alphas, prediction_score, stability_score, treeEstimators) = dynGENIE3(TS_data,time_points,SS_data=SS_data)
get_link_list(VIM,file_name='DREAM4_dynGENIE3_InSilico_Size10_1.txt')
