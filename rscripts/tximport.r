# %%
# Load libraries
library(tximport)
library(biomaRt)
library(readr)
library(dplyr)
# ---------------------------------------------------------------------------- #
#                      Step 1: Pull all salmon quant files                     #
# ---------------------------------------------------------------------------- #
# %%
quant_dir <- "/path/to/folder/containing/salmon/outputs" # <-- update this path

# Get full paths to quant.sf files inside each sample folder
sample_dirs <- list.dirs(quant_dir, recursive = FALSE, full.names = TRUE)
quant_files <- file.path(sample_dirs, "quant.sf")

# Extract sample names from directory names
sample_names <- basename(sample_dirs)
names(quant_files) <- sample_names

# ---------------------------------------------------------------------------- #
#    Step 2: Generate transcript ID to gene ID mapping for gene-level counts   #
# ---------------------------------------------------------------------------- #
# %%
ensembl <- useMart("ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl",
    host = "https://www.ensembl.org",
    verbose = TRUE
)

# Get tx2gene mapping
tx2gene <- getBM(
    attributes = c("ensembl_transcript_id_version", "ensembl_gene_id"),
    mart = ensembl,
    verbose = TRUE
)
# ---------------------------------------------------------------------------- #
#        Step 3: Run tximport, which can return both raw counts and TPM        #
# ---------------------------------------------------------------------------- #
# %%
txi <- tximport(
    files = quant_files,
    type = "salmon",
    tx2gene = tx2gene
    # countsFromAbundance = "scaledTPM"
)

write.table(txi$abundance, file = "export_tpm.tsv", sep = "\t") # TPM
write.table(txi$counts, file = "export_raw_counts.tsv", sep = "\t") # raw counts
# ---------------------------------------------------------------------------- #
#               Optional: load metadata to run DESeq2 afterwards               #
# ---------------------------------------------------------------------------- #
# %%
metadata <- read.csv("metadata.csv")
rownames(metadata) <- metadata$Run
keep <- colnames(txi$counts)
metadata_filtered <- metadata[rownames(metadata) %in% keep, ]

# %%
library(DESeq2)

dds <- DESeqDataSetFromTximport(txi, colData = metadata_filtered, design = ~tissue)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

res <- results(dds, contrast = c("tissue", "tumor", "normal")) # adjust contrast according to experiment
res <- na.omit(res)
res <- res[order(res$log2FoldChange, decreasing = TRUE), ]
