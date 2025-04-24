
# save.image("filter.rda")
# stop()

library(tidyverse)

# load dataset in EPbenchmarking format
crispr <- read_tsv(snakemake@input[[1]], show_col_types = FALSE)

# select pairs from specified datasets...
if (!is.null(snakemake@params$datasets)) {
  crispr <- filter(crispr, Dataset %in% snakemake@params$datasets)
}

# ...and cell types
if (!is.null(snakemake@params$cell_types)) {
  crispr <- filter(crispr, CellType %in% snakemake@params$cell_types)
}

# extract original column names for later formatting
crispr_cols <- colnames(crispr)

# calculate distance to TSS
crispr <- crispr %>%
  mutate(distToTSS = if_else(
    chrom == chrTSS,
    true = abs(((chromStart + chromEnd) / 2) - ((startTSS + endTSS) / 2)),
    false = NA_real_))

# get all positive and negatives
positives <- filter(crispr, Regulated == TRUE)
negatives <- filter(crispr, Regulated == FALSE)

# apply filter for positives based on distance threshold
positives <- positives %>% 
  mutate(filt = distToTSS <= snakemake@params$dist_filter)

# get proportion of filtered pairs
filt_prop <- positives %>% 
  group_by(Dataset) %>% 
  summarize(prop_pos = sum(filt == TRUE) / n())

# only retain positives passing filter
positives_filt <- select(filter(positives, filt == TRUE), -filt)

# add positive ratios per dataset to negatives and randomly sample an equal proportion of negatives
negatives_sample <- negatives %>% 
  left_join(filt_prop, by = "Dataset") %>% 
  group_by(Dataset) %>% 
  sample_n(., size = round(n() * unique(prop_pos))) %>% 
  select(-prop_pos)

# combine filtered positives and sampled negatives to form output
crispr_prop_filt <- bind_rows(positives_filt, negatives_sample) %>% 
  arrange(chrom, chromStart, chromEnd, measuredGeneSymbol) %>% 
  select(all_of(crispr_cols))

# count the number of positives and negatives per dataset
full_counts <- crispr %>% 
  filter(ValidConnection == "TRUE") %>% 
  group_by(Dataset) %>% 
  summarize(positives = sum(Regulated == TRUE), negatives = sum(Regulated == FALSE), total = n())

filtered_counts <- crispr_prop_filt %>% 
  filter(ValidConnection == "TRUE") %>% 
  group_by(Dataset) %>% 
  summarize(positives = sum(Regulated == TRUE), negatives = sum(Regulated == FALSE), total = n())

# combine into one table for output
counts <- bind_rows(full = full_counts, filterd = filtered_counts, .id = "type")

# save tables to output files
write_tsv(crispr, file = snakemake@output$full)
write_tsv(crispr_prop_filt, file = snakemake@output$filtered)
write_tsv(counts, file = snakemake@output$pair_counts)
