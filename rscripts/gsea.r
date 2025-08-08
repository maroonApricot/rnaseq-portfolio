# %%
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(metafor)
library(tidyverse)
library(msigdbr)
# -------------------------------- Annotations ------------------------------- #
# %%
species <- "Homo sapiens"

# Get gene sets from different collections
gs_h <- msigdbr(species = species, collection = "H")
gs_c2 <- msigdbr(species = species, collection = "C2", subcollection = "REACTOME") # add more as necessary

# Select gene sets of interest
selected_sets <- list(
    H = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    C2 = "REACTOME_MITOCHONDRIAL_BIOGENESIS",
)

# Subset each
subset_h <- subset(gs_h, gs_name == selected_sets$H)
subset_c2 <- subset(gs_c2, gs_name == selected_sets$C2)

# Combine into one TERM2GENE data frame
combined_sets <- rbind(subset_h, subset_c2)
term2gene <- combined_sets[, c("gs_name", "gene_symbol")] # gene_symbol, ensembl_gene

# ------------------------ reading in data from deseq2 ----------------------- #
# %%
df <- read.csv("/path/to/deseq/results.csv", header = TRUE)

original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$X
gene_list <- na.omit(original_gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)

# %%
gsea_single <- GSEA(
    geneList = gene_list,
    TERM2GENE = term2gene,
    pvalueCutoff = 1,
    verbose = TRUE,
    nPermSimple = 1000, # increase as necessary to capture highly skewed enrichment
    eps = 0
)

pdf("HS_OP.pdf", width = 8, height = 6)
gseaplot2(gsea_single, geneSetID = "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
dev.off()

res_df <- gsea_single@result
print(res_df)
