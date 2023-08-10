
configfile: "config.yaml"
include: "data_preprocessing.smk"



subworkflow workflow1:
    snakefile: "data_preprocessing.smk"
    configfile: "config.yaml"



# rule calibrate_RF_models:
#     input:
#         ancient(rules.extract_differntial_expression_data.output.DE_data)
#     output:
#         best_scores = expand(Path(CALIBRATION_DIR)/ "{method}"/ "{period}_{imput}"/ 'best_scores_{study}.npy', method='RF', period=periods, imput=imputs, study=studies),
#         best_params = expand(Path(CALIBRATION_DIR)/ "{method}"/ "{period}_{imput}"/ 'best_params_{study}.npy', method='RF', period=periods, imput=imputs, study=studies)
#     params:
#         n_replica = 10
#     message:
#         "Conduct calibration of random forest based models"
#     shell:
#         "python scripts/calibration/RF_tune.py --studies {studies} --n_replica {params.n_replica}"
# rule calibrate_ridge_models:
#     input:
#         rules.extract_differntial_expression_data.output.DE_data
#     output:
#         best_scores = expand(Path(CALIBRATION_DIR)/ "{method}"/ "{period}_{imput}"/ 'best_scores_{study}.npy', method='ridge', period=periods, imput=imputs, study=studies),
#         best_params = expand(Path(CALIBRATION_DIR)/ "{method}"/ "{period}_{imput}"/ 'best_params_{study}.npy', method='ridge', period=periods, imput=imputs, study=studies),
#     message:
#         "Conduct calibration of ridge based models"
#     shell:
#         "python scripts/calibration/ridge_tune.py --studies {studies}"
# rule plot_calibration_scores:
#     input:
#         ancient(rules.calibrate_RF_models.output.best_scores),
#         ancient(rules.calibrate_ridge_models.output.best_scores),
#     output:
#         expand(Path(CALIBRATION_DIR)/ "{method}"/ "scores.pdf", method=['RF','ridge'], period=periods, imput=imputs, study=studies),
#     message:
#         "Plot testing scores for ridge and rf models"
#     script:
#         "scripts/calibration/visualize_scores.py"
# rule infer_GRN_RF:
#     input:
#         ancient(rules.extract_differntial_expression_data.output.DE_data),
#         ancient(rules.calibrate_RF_models.output.best_params)
#     output:
#         links = expand(Path(GRN_DIR)/ "{method}"/ "Links_{period}_{imput}_{study}.csv", method='RF', period=periods, imput=imputs, study=studies),
#         links_pool = expand(Path(GRN_DIR) / "{method}" / "Links_pool_{period}_{imput}_{study}.csv", method='RF',period=periods,imput=imputs,study=studies),
#     params:
#         n_replica = 100
#     message:
#         "Conduct GRN inference using RF"
#     shell:
#         "python scripts/GRN/RF_GRN.py --studies {studies} --n_replica {params.n_replica}"
# rule infer_GRN_ridge:
#     input:
#         rules.extract_differntial_expression_data.output.DE_data,
#         rules.calibrate_ridge_models.output.best_params
#     output:
#         links = expand(Path(GRN_DIR)/ "{method}"/ "Links_{period}_{imput}_{study}.csv", method='ridge', period=periods, imput=imputs, study=studies),
#     message:
#         "Conduct GRN inference using ridge"
#     shell:
#         "python scripts/GRN/ridge_GRN.py --studies {studies}"
# rule infer_GRN_portia:
#     input:
#         rules.extract_differntial_expression_data.output.DE_data,
#     output:
#         links = expand(Path(GRN_DIR)/ "{method}"/ "Links_{period}_{imput}_{study}.csv", method='portia', period=periods, imput=imputs, study=studies)
#     message:
#         "Conduct GRN inference using portia"
#     shell:
#         "python scripts/GRN/portia_GRN.py --studies {studies}"
# rule create_baseline_models:
#     input:
#         ancient(rules.infer_GRN_RF.output.links),
#         rules.infer_GRN_ridge.output.links,
#         rules.infer_GRN_portia.output.links,
#     output:
#         expand(Path(MODELSELECTION_BASELINE_DIR) / 'ep_scores_{DE_type}.csv', DE_type=F_DE_data().keys()),
#         expand(Path(MODELSELECTION_BASELINE_DIR) / 'ep_scores_series_{DE_type}.csv', DE_type=F_DE_data().keys()),
#     params:
#         n_random_links = 1000
#     message:
#         "Creates random baseline models using inferred links and calculate ep for random models"
#     shell:
#         "python scripts/model_selection/create_baseline_models.py --GRN_methods {GRN_methods} --n_random_links {params.n_random_links}"
# rule calculate_model_selection_scores:
#     input:
#         rules.create_baseline_models.output
#     output:
#         Path(MODELSELECTION_DIR)/'scores.json'
#     message:
#         "Calculates model selections scores of early precision rate, sig flags, etc."
#     shell:
#         "python scripts/model_selection/calculate_scores.py --GRN_methods {GRN_methods} "
# rule model_selection_violin_plot:
#     input:
#         rules.calculate_model_selection_scores.output
#     output:
#         'results/model_selection/violinplot.pdf'
#     params:
#         n_random_links = 1000
#     message:
#         "Violin plot for all models"
#     shell:
#         "python scripts/model_selection/violin_plot.py  --n_random_links {params.n_random_links}"
# rule select_best_models:
#     input:
#         rules.calculate_model_selection_scores.output
#     output:
#         shortlisted_models = 'results/model_selection/shortlisted_models.txt',
#         selected_models = 'results/model_selection/selected_models.txt'
#     message:
#         "Select best models for early and late phases by filtering based on R2 score, sig_flag and using epr score"
#     shell:
#         "python scripts/model_selection/model_selection.py"
# rule line_plot_shortlisted_models:
#     input:
#         rules.calculate_model_selection_scores.output,
#         rules.select_best_models.output.selected_models
#     output:
#         'results/model_selection/lineplot.pdf',
#     message:
#         "Plots epr vs top quantiles for the shortlisted models"
#     shell:
#         "python scripts/model_selection/line_plot.py"
# rule model_selection:
#     input:
#         'results/model_selection/violinplot.pdf',
#         'results/model_selection/shortlisted_models.txt',
#         'results/model_selection/selected_models.txt',
#         'results/model_selection/lineplot.pdf',
# rule uncertainity_analysis_random_regulators:
#     input:
#         rules.extract_differntial_expression_data.input.imput_data,
#         rules.select_best_models.output.selected_models
#     output:
#         expand(Path('results/uncertainity_analysis/random_regulators') / 'ep_scores_random_{selected_model}.csv', selected_model=selected_models())
#     params:
#         n_features = 50,
#         n_repeat = 100
#     message:
#         "Runs uncertainity analysis to evaluate if DE proteins are significantly better than reandom regulators"
#     shell:
#         "python scripts/uncertainity_analysis/regulators_random/calculate_random_scores.py --n_features {params.n_features} --n_repeat {params.n_repeat}"
# rule uncertainity_analysis_random_regulators_plot_violin:
#     input:
#         rules.uncertainity_analysis_random_regulators.output,
#         rules.calculate_model_selection_scores.output
#     output:
#         Path('results/uncertainity_analysis/random_regulators') / f'violinplot.png'
#     params:
#         n_repeat = 100
#     message:
#         "Plots the results of uncertainity analysis on the DE proteins as regulators versus DE proteins"
#     shell:
#         "python scripts/uncertainity_analysis/regulators_random/violin_plot.py --n_repeat {params.n_repeat}"
# rule VSA_role_analysis:
#     input:
#         rules.select_best_models.output.selected_models,
#         rules.map_protnames_to_genenames.output
#     output:
#         vsa = expand(Path(VSA_DIR) / 'vsa_{selected_model}_{study}.csv', selected_model=selected_models(), study=studies),
#         top_role_change = expand(Path(VSA_DIR) / 'top_role_change_{selected_model}.csv', selected_model = selected_models()),
#         critical_role_change = expand(Path(VSA_DIR) / 'critical_role_change_{selected_model}.csv', selected_model = selected_models()),
#     message:
#         "Vester's sensitivity analysis to study protein roles in the network"
#     script:
#         "scripts/VSA/analyse_roles.py"
# rule uncertainity_analysis_VSA_noise:
#     input:
#         rules.select_best_models.output.selected_models,
#         rules.infer_GRN_ridge.output.links,
#         rules.infer_GRN_RF.output.links,
#         rules.infer_GRN_portia.output.links,
#         rules.VSA_role_analysis.output.top_role_change,
#         rules.VSA_role_analysis.output.critical_role_change,
#         rules.map_protnames_to_genenames.output
#     output:
#         target_genes = expand(Path("results/uncertainity_analysis/VSA_noise")/ "target_genes_{model_name}.txt", model_name=selected_models())
#     message:
#         "Robustness analysis for the results of VSA"
#     script:
#         "scripts/uncertainity_analysis/VSA_noise/analyse_noise.py"
# rule plot_VSA_analysis:
#     input:
#         rules.VSA_role_analysis.output,
#         rules.uncertainity_analysis_VSA_noise.output
#     output:
#         plot_role_changes=expand(Path("results/VSA/plots") / 'role_changes_{selected_model}.pdf',selected_model=selected_models())
#     message:
#         "Plot Vester's sensitivity analysis results with sig role change"
#     script:
#         "scripts/VSA/plot_roles.py"
# rule identify_top_links:
#     input:
#         rules.uncertainity_analysis_VSA_noise.output.target_genes,
#         rules.infer_GRN_portia.output.links,
#         rules.map_protnames_to_genenames.output,
#     output:
#         expand(Path('results/network_analysis/top_links') / 'top_links_{model_name}.csv', model_name=selected_models()),
#         expand(Path('results/network_analysis/top_links') / 'top_links_around_target_genes_{model_name}.csv', model_name=selected_models())
#     message:
#         "Identify the top links in the network and the top links around the target proteins"
#     script:
#         "scripts/analyze_network/identify_top_links.py"
# rule uncertainity_analysis_top_links:
#     input:
#         rules.identify_top_links.output
#     output:
#         expand(Path('results/uncertainity_analysis/network_noise/') / '{model_name}/links_divergence_scores.csv', model_name=selected_models()),
#     message:
#         "Runs uncertainity analysis on the links in the network and determine divergence score"
#     script:
#         "scripts/uncertainity_analysis/network_noise/calculate_divergence_scores.py"
# rule visualize_gene_network:
#     """visualize the final shortlisted network. It requires cytoscape to be opened"""
#     input:
#         rules.uncertainity_analysis_top_links.output,
#     message:
#         "Plots GRN using cytoscape and igraph, It only sends the graphs to cytoscape but has to be manually taken care of after"
#     script:
#         "scripts/network_analysis/plot_network.py"
# rule enrichment_analysis:
#     input:
#         rules.extract_differntial_expression_data.output.DE_protnames
#     output:
#         annotations = expand(Path(ENRICH_DIR)/ 'enrichment_all_{period}_{imput}.csv', period=periods, imput=imputs),
#         network = expand(Path(ENRICH_DIR) / 'network_{period}_{imput}.csv', period=periods, imput=imputs)
#     message:
#         "Conduct enrichment analysis using STRING to obtain enriched netwrok and biological functions"
#     script:
#         "scripts/enrichment_analysis/string_enquiry.py"
# rule plot_enriched_annotations:
#     input:
#         rules.enrichment_analysis.output.annotations
#     output:
#         expand(Path(ENRICH_DIR) / 'plots' / '{DE_type}.png', DE_type=[F_model_name_2_method_and_DE_type(model_name)[1] for model_name in selected_models()])
#     params:
#         top_n = 10,
#         length_limit = 40
#     message:
#         "Plots enriched annotations of GO for the selected models"
#     shell:
#         "python scripts/enrichment_analysis/plot_EA.py --top_n {params.top_n} --length_limit {params.length_limit}"
