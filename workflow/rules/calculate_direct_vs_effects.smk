## Calculate probability of direct vs. indirect CRISPR hits

# fit powerlaw models to predict direct hit rates as function of distance to TSS
rule model_direct_effects_rate:
  input:
    cis_rates = "results/{project}_analyses/cis_positive_hit_rates.tsv",
    trans_rates = "results/{project}_analyses/trans_positive_hit_rates.tsv"
  output:
    models = "results/{project}_analyses/direct_effect_models/direct_effects_model.rds",
    direct_rates_datasets = "results/{project}_analyses/direct_effect_models/direct_rates_per_dataset.tsv",
    direct_rates_average  = "results/{project}_analyses/direct_effect_models/direct_rates_average_across_datasets.tsv"
  params:
    seed = 20250508,
    max_distance = 1e6
  conda: "sceptre"
  script:
    "../scripts/model_direct_effects_rate.R"
    
# impute trans positive hit rates for CRISPR datasets without empirically computed hit rates
rule impute_trans_effects_rate:
  input:
    trans_rates = "results/{project}_analyses/trans_positive_hit_rates.tsv",
    dataset_ids = "config/crispr_dataset_ids.tsv"
  output: "results/{project}_analyses/trans_positive_hit_rates_imputed.tsv"
  conda: "sceptre"
  script:
    "../scripts/impute_trans_effects_rate.R"

# plot rates and modelled direct rates across all datasets
rule analyze_modeled_direct_rates:
  input:
    cis_rates = "results/{project}_analyses/cis_positive_hit_rates.tsv",
    trans_rates = "results/{project}_analyses/trans_positive_hit_rates.tsv",
    models = "results/{project}_analyses/direct_effect_models/direct_effects_model.rds",
    direct_rates_datasets = "results/{project}_analyses/direct_effect_models/direct_rates_per_dataset.tsv",
    direct_rates_average  = "results/{project}_analyses/direct_effect_models/direct_rates_average_across_datasets.tsv",
    dataset_ids = "config/crispr_dataset_ids.tsv"
  output: "results/{project}_analyses/model_indirect_effects.html"
  params:
    max_distance = 1e6
  conda: "sceptre"
  resources:
    mem = "16G"
  script:
    "../scripts/model_indirect_effects.Rmd"
    
# predict direct vs indirect effects probability rate for training and held-out CRISPR data
rule predict_direct_vs_indirect_effect:
  input:
    models = "results/{project}_analyses/direct_effect_models/direct_effects_model.rds",
    indirect_rates = "results/{project}_analyses/trans_positive_hit_rates_imputed.tsv",
    dataset_ids = "config/crispr_dataset_ids.tsv",
    crispr  = lambda wildcards: config["crispr_benchmark_files"][wildcards.dataset]["file"],
    gene_tss = config["tss_annot"]
  output: "results/annotated_crispr_data_{project}/EPCrisprBenchmark_{dataset}_element_classes_direct_effects.tsv.gz"
  params:
    dist_col = lambda wildcards: config["crispr_benchmark_files"][wildcards.dataset]["dist_col"],
    include_imputed = True
  conda: "sceptre"
  script:
    "../scripts/predict_direct_vs_indirect_effect.R"
    
# predict direct vs indirect effects probability rate for DC-TAP-seq data used in paper
rule predict_direct_vs_indirect_effect_dctap:
  input:
    crispr = "resources/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_fdr20.tsv",
    models = "results/dctap_analyses/direct_effect_models/direct_effects_model.rds",
    indirect_rates = "results/dctap_analyses/trans_positive_hit_rates_imputed.tsv",
    dataset_ids = "config/crispr_dataset_ids.tsv",
    gene_tss = config["tss_annot"]
  output: "results/annotated_crispr_data_dctap/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_fdr20_direct_effects.tsv"
  params:
    dist_col = "distance_to_gencode_gene_TSS"
  conda: "sceptre"
  script:
    "../scripts/predict_direct_vs_indirect_effect_dctap.R"
