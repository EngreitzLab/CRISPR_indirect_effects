## Replace old element categories in CRISPR data with updated categories

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(here)
})

# training and held-out CRISPR data files
training_file <- "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_main/2025_0227_CTCF_and_H3K27ac/results/2025_0516_fix_peak_sizes/categorized_pairs/training.h3k27me3_quantile_50.ratio_quantile_35.h3k27ac_quantiles_60_90.WTC11_as_WTC11.tsv.gz"
heldout_file <- "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_main/2025_0227_CTCF_and_H3K27ac/results/2025_0516_fix_peak_sizes/categorized_pairs/validation_all.h3k27me3_quantile_50.ratio_quantile_35.h3k27ac_quantiles_60_90.WTC11_as_WTC11.tsv.gz"

# file containing final element categories
categories_file <- "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_main/2025_0227_CTCF_and_H3K27ac/results/2025_0620_revisit_categories/recategorized_v2/crispr_for_distal_regulation.recategorized.tsv.gz"

# load all input files
training <- read_tsv(training_file, show_col_types = FALSE)
heldout <- read_tsv(heldout_file, show_col_types = FALSE)
categories <- read_tsv(categories_file, show_col_types = FALSE)

# element category columns to add the CRISPR data
category_cols <- c("CTCF_peak_overlap", "H3K27me3_peak_overlap", "H3K4me1_peak_overlap",
                   "H3K27ac_peak_overlap", "CTCF.RPM", "DHS.RPM", "H3K27ac.RPM", "H3K27me3.RPM",               
                   "H3K4me1.RPM", "DHS.RPM.expandedRegion", "H3K27ac.RPM.expandedRegion",
                   "H3K27me3.RPM.expandedRegion", "H3K4me1.RPM.expandedRegion", "elementName",                
                   "distanceToTSS", "CTCF.H3K27ac.ratio", "element_category")

# extract columns to merge with crispr data
categories_to_add <- categories %>% 
  select(chrom = chr, chromStart = start, chromEnd = end, measuredGeneSymbol,
         CellType = cell_type, data_category, all_of(category_cols))

# old category columns to replace in CRISPR data files
cols_remove <- c("CTCF.H3K27ac.expandedRegion.ratio", "element_category_with_dnase",
                 "element_category_simple", category_cols)

# columns to use to merge new categories into CRISPR data tables
merge_cols <- c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol", "CellType", "data_category")

# replace old category columns with new ones for each CRISPR dataset
training <- training %>% 
  select(-any_of(cols_remove)) %>% 
  left_join(categories_to_add, by = merge_cols) %>% 
  select(-data_category)

heldout <- heldout %>% 
  select(-any_of(cols_remove)) %>% 
  left_join(categories_to_add, by = merge_cols) %>% 
  select(-data_category)

# write reformatted CRISPR datasets to output files
write_tsv(training, file = here("resources/training.h3k27me3_quantile_50.ratio_quantile_35.h3k27ac_quantiles_60_90.WTC11_as_WTC11.reformatted.tsv.gz"))
write_tsv(heldout, file = here("resources/validation_all.h3k27me3_quantile_50.ratio_quantile_35.h3k27ac_quantiles_60_90.WTC11_as_WTC11.reformatted.tsv.gz"))
