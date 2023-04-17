import os
from subprocess import run
from pathlib import Path
import numpy as np
from scripts.imports import MAIN_DIR, DATA_DIR, STATISTICAL_ANALYSIS_DIR, ENRICH_DIR, CALIBRATION_DIR, GRN_DIR, \
    MODELSELECTION_DIR, VSA_DIR, VSA_NOISE_DIR, GRN_VISUALIZE_DIR, F_selected_models, F_model_name_2_method_and_DE_type, \
    F_DE_data,RANDOM_REGULATORS_DIR, MODELSELECTION_BASELINE_DIR

periods = ['early', 'late']
periods_days = ['day1_11', 'day1_21']
imputs = ['MinProb', 'KNN']
studies = ['ctr', 'mg']
GRN_methods = ['RF', 'ridge', 'portia']

def selected_models():
    if os.path.exists(Path(MODELSELECTION_DIR)/ 'selected_models.txt'):
        return F_selected_models()
    else:
        return []
rule edit_original_proteomics_data:
    input:
        Path(MAIN_DIR)/'data'/'original_omics.xlsx'
    output:
        Path(DATA_DIR)/'edited_data.csv'
    shell:
        "python scripts/edit_data/edit.py --original_df_dir {input} --df_dir {output}"
rule extract_differntial_expression_data:
    input:
        imput_data = expand(Path(STATISTICAL_ANALYSIS_DIR) / '{period}' / 'ProteinAbundance_tables/ProtAbundance__Norm_n_Imp_{imput}.csv', period=periods_days, imput=imputs),
        sig_table = expand(Path(STATISTICAL_ANALYSIS_DIR) / '{period}' / 'DiffExp_tables/TopTable_TimeCourse_{imput}.csv', period=periods_days, imput=imputs),
    output:
        data = expand(Path(DATA_DIR)/"data_{period}_{imput}_{study}.csv", period=periods, imput=imputs, study=studies),
        DE_protnames = Path(DATA_DIR)/'DE_protnames.txt',
        DE_data = Path(DATA_DIR) / 'DE_data.csv'
    message:
        "Extract differential expressed proteins and associated expression data"
    # params:
    #     prefix=lambda wildcards, output: output[0][:-4]
    priority: 0
    threads: 1
    shell:
        "python scripts/post_statistical_analysis/process_PSA.py --periods {periods} --periods_days {periods_days} --imputs {imputs} --studies {studies}"
rule enrichment_analysis:
    input:
        rules.extract_differntial_expression_data.output.DE_protnames
    output:
        annotations = expand(Path(ENRICH_DIR)/ 'enrichment_all_{period}_{imput}.csv', period=periods, imput=imputs),
        network = expand(Path(ENRICH_DIR) / 'network_{period}_{imput}.csv', period=periods, imput=imputs)
    message:
        "Conduct enrichment analysis using STRING to obtain enriched netwrok and biological functions"
    script:
        "scripts/enrichment_analysis/string_enquiry.py"
rule calibrate_RF_models:
    input:
        ancient(rules.extract_differntial_expression_data.output.DE_data)
    output:
        best_scores = expand(Path(CALIBRATION_DIR)/ "{method}"/ "{period}_{imput}"/ 'best_scores_{study}.npy', method='RF', period=periods, imput=imputs, study=studies),
        best_params = expand(Path(CALIBRATION_DIR)/ "{method}"/ "{period}_{imput}"/ 'best_params_{study}.npy', method='RF', period=periods, imput=imputs, study=studies)
    params:
        n_replica = 10
    message:
        "Conduct calibration of random forest based models"
    shell:
        "python scripts/calibration/RF_tune.py --studies {studies} --n_replica {params.n_replica}"
rule calibrate_ridge_models:
    input:
        rules.extract_differntial_expression_data.output.DE_data
    output:
        best_scores = expand(Path(CALIBRATION_DIR)/ "{method}"/ "{period}_{imput}"/ 'best_scores_{study}.npy', method='ridge', period=periods, imput=imputs, study=studies),
        best_params = expand(Path(CALIBRATION_DIR)/ "{method}"/ "{period}_{imput}"/ 'best_params_{study}.npy', method='ridge', period=periods, imput=imputs, study=studies),
    message:
        "Conduct calibration of ridge based models"
    shell:
        "python scripts/calibration/ridge_tune.py --studies {studies}"
rule plot_calibration_scores:
    input:
        ancient(rules.calibrate_RF_models.output.best_scores),
        ancient(rules.calibrate_ridge_models.output.best_scores),
    output:
        expand(Path(CALIBRATION_DIR)/ "{method}"/ "scores.pdf", method=['RF','ridge'], period=periods, imput=imputs, study=studies),
    message:
        "Plot testing scores for ridge and rf models"
    script:
        "scripts/calibration/visualize_scores.py"
rule infer_GRN_RF:
    input:
        ancient(rules.extract_differntial_expression_data.output.DE_data),
        ancient(rules.calibrate_RF_models.output.best_params)
    output:
        links = expand(Path(GRN_DIR)/ "{method}"/ "Links_{period}_{imput}_{study}.csv", method='RF', period=periods, imput=imputs, study=studies),
        links_pool = expand(Path(GRN_DIR) / "{method}" / "Links_pool_{period}_{imput}_{study}.csv", method='RF',period=periods,imput=imputs,study=studies),
    params:
        n_replica = 100
    message:
        "Conduct GRN inference using RF"
    shell:
        "python scripts/GRN/RF_GRN.py --studies {studies} --n_replica {params.n_replica}"
rule infer_GRN_ridge:
    input:
        rules.extract_differntial_expression_data.output.DE_data,
        rules.calibrate_ridge_models.output.best_params
    output:
        links = expand(Path(GRN_DIR)/ "{method}"/ "Links_{period}_{imput}_{study}.csv", method='ridge', period=periods, imput=imputs, study=studies),
    message:
        "Conduct GRN inference using ridge"
    shell:
        "python scripts/GRN/ridge_GRN.py --studies {studies}"
rule infer_GRN_portia:
    input:
        rules.extract_differntial_expression_data.output.DE_data,
    output:
        links = expand(Path(GRN_DIR)/ "{method}"/ "Links_{period}_{imput}_{study}.csv", method='portia', period=periods, imput=imputs, study=studies)
    message:
        "Conduct GRN inference using portia"
    shell:
        "python scripts/GRN/portia_GRN.py --studies {studies}"
rule create_baseline_models:
    input:
        ancient(rules.infer_GRN_RF.output.links),
        rules.infer_GRN_ridge.output.links,
        rules.infer_GRN_portia.output.links,
    output:
        expand(Path(MODELSELECTION_BASELINE_DIR) / 'ep_scores_{DE_type}.csv', DE_type=F_DE_data().keys()),
        expand(Path(MODELSELECTION_BASELINE_DIR) / 'ep_scores_series_{DE_type}.csv', DE_type=F_DE_data().keys()),
    params:
        n_random_links = 1000
    message:
        "Creates random baseline models using inferred links and calculate ep for random models"
    shell:
        "python scripts/model_selection/create_baseline_models.py --GRN_methods {GRN_methods} --n_random_links {params.n_random_links}"
rule calculate_model_selection_scores:
    input:
        rules.create_baseline_models.output
    output:
        Path(MODELSELECTION_DIR)/'scores.json'
    message:
        "Calculates model selections scores of early precision rate, sig flags, etc."
    shell:
        "python scripts/model_selection/calculate_scores.py --GRN_methods {GRN_methods} "
rule model_selection_violin_plot:
    input:
        rules.calculate_model_selection_scores.output
    output:
        'results/model_selection/violinplot.pdf'
    params:
        n_random_links = 1000
    message:
        "Violin plot for all models"
    shell:
        "python scripts/model_selection/violin_plot.py  --n_random_links {params.n_random_links}"
rule select_best_models:
    input:
        rules.calculate_model_selection_scores.output
    output:
        shortlisted_models = 'results/model_selection/shortlisted_models.txt',
        selected_models = 'results/model_selection/selected_models.txt'
    message:
        "Select best models for early and late phases by filtering based on R2 score, sig_flag and using epr score"
    shell:
        "python scripts/model_selection/model_selection.py"
rule line_plot_shortlisted_models:
    input:
        rules.calculate_model_selection_scores.output,
        rules.select_best_models.output.selected_models
    output:
        'results/model_selection/lineplot.pdf',
    message:
        "Plots epr vs top quantiles for the shortlisted models"
    shell:
        "python scripts/model_selection/line_plot.py"
rule model_selection:
    input:
        'results/model_selection/violinplot.pdf',
        'results/model_selection/shortlisted_models.txt',
        'results/model_selection/selected_models.txt',
        'results/model_selection/lineplot.pdf',
rule map_protnames_to_genenames:
    input:
        rules.extract_differntial_expression_data.output.DE_protnames
    output:
        Path(DATA_DIR)/ 'protnames_to_genenames.json'
    message:
        "Maps protnames to genenames using https://rest.uniprot.org"
    script:
        "scripts/protname_to_genename.py"
rule VSA_role_analysis:
    input:
        rules.select_best_models.output.selected_models,
        rules.infer_GRN_ridge.output.links,
        rules.infer_GRN_RF.output.links,
        rules.infer_GRN_portia.output.links,
        rules.map_protnames_to_genenames.output
    output:
        vsa = expand(Path(VSA_DIR) / 'vsa_{selected_model}_{study}.csv', selected_model=selected_models(), study=studies),
        top_role_change = expand(Path(VSA_DIR) / 'top_role_change_{selected_model}.csv', selected_model = selected_models()),
        critical_role_change = expand(Path(VSA_DIR) / 'critical_role_change_{selected_model}.csv', selected_model = selected_models()),
        plot_role_changes = expand(Path(VSA_DIR) / 'role_changes_{selected_model}.pdf', selected_model = selected_models())
    params:
        top_quantile_role_change = 0.9 # top 10 percent with the largest role change are chosen for top role change
    message:
        "Vester's sensitivity analysis to study protein roles in the network"
    shell:
        "python scripts/VSA/analyse_roles.py --studies {studies} --top_quantile_role_change {params.top_quantile_role_change}"
rule uncertainity_analysis_VSA_noise:
    input:
        rules.select_best_models.output.selected_models,
        rules.infer_GRN_ridge.output.links,
        rules.infer_GRN_RF.output.links,
        rules.infer_GRN_portia.output.links,
        rules.VSA_role_analysis.output.top_role_change,
        rules.VSA_role_analysis.output.critical_role_change,
        rules.map_protnames_to_genenames.output
    output:
        target_genes = Path(VSA_NOISE_DIR)/ "target_genes.json"
    params:
        n_noisy_datasets=100,# number of noisy datasets to be created for each noise type and std
        warm_start=True,#To use previous results
    message:
        "Robustness analysis for the results of VSA"
    shell:
        "python scripts/uncertainity_analysis/VSA_noise/analyse_noise.py --studies {studies} --warm_start {params.warm_start} --n_noisy_datasets {params.n_noisy_datasets}"
rule visualize_protein_network:
    input:
        rules.uncertainity_analysis_VSA_noise.output.target_genes
    output:
        expand(Path(GRN_VISUALIZE_DIR) / "GRN_{model_name}_{study}.pdf", model_name=selected_models(), study=studies),
    params:
        top_n_links = 2
    message:
        "Plots protein-protein connections for the selected proteins"
    shell:
        "python scripts/GRN/visualize_network.py --studies {studies} --top_n_links {params.top_n_links}"
rule plot_enriched_annotations:
    input:
        rules.enrichment_analysis.output.annotations
    output:
        expand(Path(ENRICH_DIR) / 'plots' / '{DE_type}.png', DE_type=[F_model_name_2_method_and_DE_type(model_name)[1] for model_name in selected_models()])
    params:
        top_n = 10,
        length_limit = 40
    message:
        "Plots enriched annotations of GO for the selected models"
    shell:
        "python scripts/enrichment_analysis/plot_EA.py --top_n {params.top_n} --length_limit {params.length_limit}"
rule uncertainity_analysis_random_regulators:
    input:
        rules.extract_differntial_expression_data.input.imput_data,
        rules.select_best_models.output.selected_models
    output:
        expand(Path(RANDOM_REGULATORS_DIR) / 'ep_scores_random_{selected_model}.csv', selected_model=selected_models())
    params:
        n_features = 50,
        n_repeat = 100
    message:
        "Runs uncertainity analysis on the DE proteins as regulators by randomly sampling proteins as regulating and calculating early precision scores"
    shell:
        "python scripts/uncertainity_analysis/regulators_random/calculate_random_scores.py --n_features {params.n_features} --n_repeat {params.n_repeat}"
rule uncertainity_analysis_plot_violin:
    input:
        rules.uncertainity_analysis_random_regulators.output,
        rules.calculate_model_selection_scores.output
    output:
        Path(RANDOM_REGULATORS_DIR) / f'violinplot.png'
    params:
        n_repeat = 100
    message:
        "Plots the results of uncertainity analysis on the DE proteins as regulators versus DE proteins"
    shell:
        "python scripts/uncertainity_analysis/regulators_random/violin_plot.py --n_repeat {params.n_repeat}"

rule all:
    input:
        Path(VSA_NOISE_DIR) / "target_genes.json",
        expand(Path(GRN_VISUALIZE_DIR) / "GRN_{model_name}_{study}.pdf",model_name=selected_models(),study=studies),
