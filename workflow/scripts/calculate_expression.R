## Compute average expression per gene from a gene expression matrix

# save.image("calc_expr.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

# load gene expression data
dge <- fread(snakemake@input[[1]])
genes <- dge[[1]]
dge <- as(data.matrix(dge[, -1]), "sparseMatrix")
rownames(dge) <- genes

# calculate average UMIs/cell
avg_umi <- rowMeans(dge)
avg_umi <- data.table(gene = names(avg_umi), avg_expr = avg_umi)

# save to output file
fwrite(avg_umi, snakemake@output[[1]], quote = FALSE, na = "NA")
