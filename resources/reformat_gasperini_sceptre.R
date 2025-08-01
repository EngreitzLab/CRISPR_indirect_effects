## Reformat file for Gasperini et al., 2019 data analyzed with sceptre to correct format to annotate
## with direct vs. indirect effect calculations

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(here)
})

# load Gasperini CRISPR data
crispr_all_file <- "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_main/2025_0227_CTCF_and_H3K27ac/results/2025_0516_fix_peak_sizes/categorized_pairs/no_power_filter/Gasperini_all_power.h3k27me3_quantile_50.ratio_quantile_35.h3k27ac_quantiles_60_90.tsv.gz"
crispr_filt_file <- "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_main/2025_0227_CTCF_and_H3K27ac/results/2025_0516_fix_peak_sizes/categorized_pairs/no_power_filter/Gasperini_DistalElementGene_SCEPTRE_well_powered.es15_over_80pct.h3k27me3_quantile_50.ratio_quantile_35.h3k27ac_quantiles_60_90.tsv.gz"
crispr_all <- read_tsv(crispr_all_file, show_col_types = FALSE)
crispr_filt <- read_tsv(crispr_filt_file, show_col_types = FALSE)

# reformat tables
crispr_all <- crispr_all %>% 
  mutate(Dataset = "Gasperini2019") %>% 
  rename(Significant = significant, EffectSize = pct_change_effect_size,
         chrom = targeting_chr_hg38, chromStart = targeting_start_hg38,
         chromEnd = targeting_end_hg38, measuredGeneSymbol = gene_symbol)

crispr_filt <- crispr_filt %>% 
  mutate(Dataset = "Gasperini2019") %>% 
  rename(Significant = significant, EffectSize = pct_change_effect_size,
         chrom = targeting_chr_hg38, chromStart = targeting_start_hg38,
         chromEnd = targeting_end_hg38, measuredGeneSymbol = gene_symbol)

# save to output files
outfile_all <- here("resources/Gasperini_all_power.h3k27me3_quantile_50.ratio_quantile_35.h3k27ac_quantiles_60_90.reformatted.tsv.gz")
outfile_filt <- here("resources/Gasperini_DistalElementGene_SCEPTRE_well_powered.es15_over_80pct.h3k27me3_quantile_50.ratio_quantile_35.h3k27ac_quantiles_60_90.reformatted.tsv.gz")
write_tsv(crispr_all, file = outfile_all)
write_tsv(crispr_filt, file = outfile_filt)

