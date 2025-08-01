## Annotate cis-results with element and pair annotations for DC-TAP-seq experiments

# save.image("pairs.rda")
# stop()

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# load annotated results
cis_results <- readRDS(snakemake@input$cis_results)
encode_results <- fread(snakemake@input$encode_file)

# only retain data on the given sample
sample <- sub("_DC_TAPseq", "", snakemake@wildcards$sample)
encode_results <- filter(encode_results, cell_type == sample)

# add valid connection columns and distance to TSS to cis-results
cis_results <- encode_results %>% 
  mutate(pair_uid = paste0(gene_id, "|", design_file_target_name)) %>% 
  select(grna_target = design_file_target_name, response_id = gene_id, pair_uid, gene_symbol,
         dist_to_tss = distance_to_gencode_gene_TSS, element_location,
         valid_connection = all_of(snakemake@params$valid_col)) %>% 
  distinct() %>% 
  left_join(cis_results, ., by = c("grna_target", "response_id"))

# get all valid elements that do not overlap any annotated promoters
valid_elements <- cis_results %>% 
  filter(element_location == "distal") %>%
  pull(grna_target) %>% 
  unique()

# create output with column labeling valid elements
output <- mutate(cis_results, valid_element = grna_target %in% valid_elements)

# save output to file
fwrite(output, file = snakemake@output[[1]], sep = ",", quote = FALSE, na = "NA")
