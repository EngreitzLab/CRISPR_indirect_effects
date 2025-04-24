## Run Sceptre for trans-acting effects

# save.image(paste0("trans_effects_", snakemake@wildcards$sample, ".rda"))
# stop()

set.seed(snakemake@params$seed)

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(rtracklayer)
  library(sceptre)
})

message("\nLoading input data")

# load table linking guides to perturbation targets
gRNA_groups_table <- read.table(snakemake@input$guides, header = TRUE, sep = '\t')
colnames(gRNA_groups_table)[1:2] <- c("grna_id", "grna_target")

# load gene UMI counts matrix
gene_matrix <- fread(snakemake@input$dge)
gene_matrix <- as(as.matrix(gene_matrix, rownames = 1), "sparseMatrix")
invisible(gc())

# load gRNA UMI counts matrix
gRNA_matrix <- fread(snakemake@input$pert)
gRNA_matrix <- as(as.matrix(gRNA_matrix, rownames = 1), "sparseMatrix")

# load and format cell metadata for sceptre
if (!is.null(snakemake@input$cell_meta)) {
  cell_metadata <- read.table(snakemake@input$cell_meta, header = TRUE, sep = '\t')
  cell_metadata$cell_batches <- as.factor(cell_metadata$cell_batches)
}

# load and format cell metadata for sceptre
# cell_metadata <- read.table(snakemake@input$cell_meta, header = TRUE, sep = '\t')
# extra_covariates <- data.frame(batch = as.factor(cell_metadata$cell_batches))
# rownames(extra_covariates) <- cell_metadata$cell_barcode

# load genome annotations and get chromosome of each gene
annot <- import(snakemake@input$genome)
genes <- annot %>% 
  as.data.frame() %>% 
  select(response_chr = seqnames, response_id = gene_id) %>% 
  mutate(response_id = sub("\\..+$", "", response_id)) %>% 
  distinct()

# remove Y chromosome genes, since some of them are annotated to both X and Y
genes <- filter(genes, response_chr != "chrY")

# load list of genes in cis-analysis
encode_results <- fread(snakemake@input$encode_file)
cis_results <- fread(snakemake@input$cis_results)

# get gene symbol for every gene id from processed ENCODE data
genes_symbols <- encode_results %>%
  select(gene = measuredEnsemblID, gene_symbol = measuredGeneSymbol) %>% 
  distinct()

# add gene symbols to cis results 
cis_results <- cis_results %>% 
  left_join(genes_symbols, by = c("response_id" = "gene")) %>% 
  filter(!is.na(gene_symbol))

# add ValidConnection column to cis results
cis_results <- cis_results %>% 
  mutate(name = paste0(gene_symbol, "|", grna_target, ":.")) %>% 
  left_join(distinct(select(encode_results, name, ValidConnection)), by = "name")

# only retain genes that were used in valid pairs in cis-analysis
cis_results_filt <- filter(cis_results, ValidConnection == "TRUE")
cis_genes <- unique(cis_results_filt[[1]])

message("Creating sceptre object")

# create sceptre object
if (!is.null(snakemake@input$cell_meta)) {
  sceptre_object <- import_data(
    response_matrix = gene_matrix,
    grna_matrix = gRNA_matrix,
    grna_target_data_frame = gRNA_groups_table,
    moi = snakemake@params$moi,
    extra_covariates = cell_metadata
  )
} else {
  sceptre_object <- import_data(
    response_matrix = gene_matrix,
    grna_matrix = gRNA_matrix,
    grna_target_data_frame = gRNA_groups_table,
    moi = snakemake@params$moi
  )
}

# create all perturbation-gene pairs
pairs <- construct_trans_pairs(sceptre_object)

# add chromosomes for each gene in all pairs table and information whether a gene was cis analysis
pairs <- pairs %>% 
  left_join(genes, by = "response_id") %>% 
  mutate(cis_gene = response_id %in% cis_genes)

# only retain pairs involving enhancers and genes on other chromosomes than the enhancer
trans_pairs <- pairs %>% 
  filter(grepl("^chr.+", grna_target)) %>% 
  mutate(grna_target_chr = sub("^(chr.+):.+", "\\1", grna_target)) %>% 
  filter(grna_target_chr != response_chr)

# filter for genes in cis analysis only if specified
if (snakemake@params$include_genes == "cis") {
  message("Building trans-pairs using genes from cis analysis")
  trans_pairs <- filter(trans_pairs, cis_gene == TRUE)
} else if (snakemake@params$include_genes == "all") {
  message("Building trans-pairs using all available genes")
} else {
  stop("Invalid 'include_genes' parameter", call. = FALSE)
}

# count how many genes are available per perturbation
available_genes_per_pert <- count(trans_pairs, grna_target, name = "available_genes")

# pick up to 100 trans genes per enhancer perturbation
message("Sampling ", snakemake@params$genes_per_pert, " trans genes per perturbation")
trans_pairs_sampled <- trans_pairs %>% 
  group_by(grna_target) %>% 
  sample_n(size = snakemake@params$genes_per_pert) %>% 
  ungroup()

# count how many genes are tested per perturbation and add to available gene counts
tested_genes_per_pert <- count(trans_pairs_sampled, grna_target, name = "tested_genes")
genes_per_pert <- left_join(available_genes_per_pert, tested_genes_per_pert, by = "grna_target")

# raise warning if for some perturbations less than 100 genes are available
if (any(genes_per_pert$tested_genes < 100)) {
  warning("Fewer than 100 genes available for testing for some perturbations", call. = FALSE)
}

# set analysis parameters
sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = select(trans_pairs_sampled, grna_target, response_id),
  side = snakemake@params$side,
  grna_integration_strategy = snakemake@params$grna_integration_strategy,
)

message("\nRunning sceptre sceptre:\n")

# assign guides to cells using specified approach
sceptre_object <- assign_grnas(
  sceptre_object,
  method = snakemake@params$guide_assignment_method,
  threshold = snakemake@params$guide_threshold
)

# run QC and filtering
sceptre_object <- run_qc(sceptre_object)

# run calibration check using negative control perturbations
guide_types <- unique(sceptre_object@grna_target_data_frame$gRNA_type)
if ("negative_control" %in% guide_types | "negative_targeting" %in% guide_types) {
  sceptre_object <- run_calibration_check(sceptre_object, parallel = FALSE)
}

# run discovery analysis to test for trans-acting effects
sceptre_object <- run_discovery_analysis(sceptre_object, parallel = FALSE)

# write all results to output directory
write_outputs_to_directory(sceptre_object, dirname(snakemake@output$results))
fwrite(genes_per_pert, file = snakemake@output$genes_per_pert, quote = FALSE, na = "NA")
