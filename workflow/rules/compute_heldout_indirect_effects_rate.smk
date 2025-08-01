## Perform differential expression tests for trans-acting interactions to compute expected rate
## of indirect effects for held-out CRISPR datasets

# annotate cis results for DC-TAP-seq datasets
rule annotate_cis_results_dctap:
  input: 
    cis_results = lambda wildcards: config["datasets"][wildcards.sample]["cis_results"],
    encode_file = lambda wildcards: config["datasets"][wildcards.sample]["encode_file"]
  output: "results/{sample}/annotated_cis_results.tsv.gz"
  params:
    valid_col = "DistalElement_Gene"
  wildcard_constraints:
    sample = "K562_DC_TAPseq|WTC11_DC_TAPseq"
  conda: "sceptre"
  script:
    "../scripts/annotate_cis_results_dctap.R"
    
# annotate cis results for any other ENCODE dataset
rule annotate_cis_results_encode:
  input: 
    cis_results = lambda wildcards: config["datasets"][wildcards.sample]["cis_results"],
    encode_file = lambda wildcards: config["datasets"][wildcards.sample]["encode_file"]
  output: "results/{sample}/annotated_cis_results.tsv.gz"
  conda: "sceptre"
  script:
    "../scripts/annotate_cis_results_encode.R"

# perform de tests for trans-effects
rule perform_trans_acting_de_tests:
  input:
    sceptre_object = lambda wildcards: config["datasets"][wildcards.sample]["sceptre_object"],
    genome = lambda wildcards: config["datasets"][wildcards.sample]["genome"],
    cis_results = "results/{sample}/annotated_cis_results.tsv.gz"
  output: 
    results = "results/{sample}/results_run_discovery_analysis.rds",
    gene_expr = "results/{sample}/gene_expr.csv",
    genes_per_pert = "results/{sample}/genes_per_perturbation.csv"
  log: "results/{sample}/perform_trans_acting_de_tests.log"
  params:
    seed = 20240909,
    include_genes = "cis",  # one of 'cis' or 'all'
    genes_per_pert = 100,
    guide_assignment_method = lambda wildcards: config["datasets"][wildcards.sample]["guide_assignment_method"],
    guide_threshold = lambda wildcards: config["datasets"][wildcards.sample]["guide_threshold"],
    side = "both",
    grna_integration_strategy = "union"
  conda: "sceptre"
  resources:
    mem_mb = "128000",
    runtime = "12h"
  script:
    "../scripts/perform_trans_acting_de_tests.R"

# calculate positive rates for cis and trans effects
rule calculate_positive_hit_rates:
  input:
    trans_results = "results/{sample}/results_run_discovery_analysis.rds",
    cis_results = "results/{sample}/annotated_cis_results.tsv.gz",    
    gene_univ = config["tss_annot"]
  output:
    cis_rate = "results/{sample}/cis_positive_rate.tsv",
    trans_rate = "results/{sample}/trans_positive_rate.tsv"
  params:
    bin_size = 5e4
  conda: "sceptre"
  script:
    "../scripts/calculate_positive_hit_rates.R"

# analyze indirect effects
rule analyze_indirect_effects:
  input:
    trans_results = "results/{sample}/results_run_discovery_analysis.rds",
    cis_results = "results/{sample}/annotated_cis_results.tsv.gz",
    gene_univ = config["tss_annot"],
    gene_expr = "results/{sample}/gene_expr.csv"
  output: "results/{sample}/analyze_indirect_effects.html"
  params:
    seed = 20240909,
    dist_threshold = lambda wildcards: config["datasets"][wildcards.sample]["dist_threshold"]
  conda: "sceptre"
  resources:
    mem_mb = "16000",
    runtime = "6h"
  script:
    "../scripts/analyze_indirect_effects.Rmd"
