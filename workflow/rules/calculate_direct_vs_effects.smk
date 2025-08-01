## Calculate probability of direct vs. indirect CRISPR hits

# fit powerlaw models to predict direct hit rates as function of distance to TSS
rule model_direct_effects_rate:
  input:
    cis_rates = "results/cis_positive_hit_rates.tsv",
    trans_rates = "results/trans_positive_hit_rates.tsv"
  output:
    models = "results/direct_effect_models/direct_effects_model.rds",
    direct_rates_datasets = "results/direct_effect_models/direct_rates_per_datasets.tsv",
    direct_rates_average  = "results/direct_effect_models/direct_rates_average_across_datasets.tsv"
  params:
    seed = 20250508,
    max_distance = 1e6
  conda: "sceptre"
  script:
    "../scripts/model_direct_effects_rate.R"
    
# impute trans positive hit rates for CRISPR datasets without empirically computed hit rates
rule impute_trans_effects_rate:
  input:
    trans_rates = "results/trans_positive_hit_rates.tsv",
    dataset_ids = "config/crispr_dataset_ids.tsv"
  output: "results/trans_positive_hit_rates_imputed.tsv"
  conda: "sceptre"
  script:
    "../scripts/impute_trans_effects_rate.R"
    
# predict direct vs indirect effects probability rate for training and held-out CRISPR data
rule predict_direct_vs_indirect_effect:
  input:
    models = "results/direct_effect_models/direct_effects_model.rds",
    indirect_rates = "results/trans_positive_hit_rates_imputed.tsv",
    dataset_ids = "config/crispr_dataset_ids.tsv",
    crispr  = lambda wildcards: config["crispr_benchmark_files"][wildcards.dataset]["file"],
    gene_tss = config["tss_annot"]
  output: "results/annotated_crispr_data/EPCrisprBenchmark_{dataset}_element_classes_direct_effects.tsv.gz"
  params:
    dist_col = lambda wildcards: config["crispr_benchmark_files"][wildcards.dataset]["dist_col"],
    include_imputed = True
  conda: "sceptre"
  script:
    "../scripts/predict_direct_vs_indirect_effect.R"

# plot rates and modelled direct rates across all datasets
rule analyze_modeled_direct_rates:
  input:
    cis_rates = "results/cis_positive_hit_rates.tsv",
    trans_rates = "results/trans_positive_hit_rates.tsv",
    models = "results/direct_effect_models/direct_effects_model.rds",
    direct_rates_datasets = "results/direct_effect_models/direct_rates_per_datasets.tsv",
    direct_rates_average  = "results/direct_effect_models/direct_rates_average_across_datasets.tsv",
    dataset_ids = "config/crispr_dataset_ids.tsv"
  output: "results/model_indirect_effects.html"
  params:
    max_distance = 1e6
  conda: "sceptre"
  resources:
    mem = "16G"
  script:
    "../scripts/model_indirect_effects.Rmd"
