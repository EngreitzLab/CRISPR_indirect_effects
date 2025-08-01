## Reformat DC-TAP-seq CRISPR data file before merging elements to correct format to annotate with
## direct vs. indirect effect calculations

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(here)
})

# load DC-TAP-seq CRISPR data
crispr_file <- "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_main/2025_0227_CTCF_and_H3K27ac/results/2025_0516_fix_peak_sizes/categorized_pairs/no_power_filter/all_dc_tap.h3k27me3_quantile_50.ratio_quantile_35.h3k27ac_quantiles_60_90.WTC11_as_WTC11.tsv.gz"
crispr <- read_tsv(crispr_file, show_col_types = FALSE)

# reformat table
crispr <- crispr %>% 
  mutate(Dataset = paste0(cell_type, "_DC_TAP")) %>% 
  rename(Significant = significant, EffectSize = pct_change_effect_size,
         chrom = resized_merged_targeting_chr_hg38, chromStart = resized_merged_targeting_start_hg38,
         chromEnd = resized_merged_targeting_end_hg38, measuredGeneSymbol = gene_symbol)

# save to output file
output_file <- here("resources/all_dc_tap.h3k27me3_quantile_50.ratio_quantile_35.h3k27ac_quantiles_60_90.WTC11_as_WTC11.reformatted.tsv.gz")
write_tsv(crispr, file = output_file)
