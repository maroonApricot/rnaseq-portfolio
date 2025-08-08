# %%
library(arrayQualityMetrics)
library(GEOquery)
library(limma)
library(oligo)
library(annotate)
library(readr)
library(dplyr)
library(tidyr)
library(biomaRt)
library(stringr)
library(ggplot2)
library(devtools)
devtools::install_github("bmbolstad/preprocessCore",
    dependencies = TRUE,
    upgrade = "always",
    configure.args = "--disable-threading"
)

# %% prep data
accession <- "GSE13270"
path <- sprintf("/home/hjang620/20250603_Anh_HCC/microarray_micenash/%s", accession)
cel_path <- sprintf("%s/CEL", path)
gse <- getGEO(accession, destdir = path, AnnotGPL = TRUE)
gse <- gse[[1]]
getGEOSuppFiles(accession, baseDir = path, makeDirectory = FALSE)
untar(sprintf("%s/%s_RAW.tar", path, accession), exdir = cel_path)
metadata <- pData(gse)
annotation <- fData(gse)

# %% pick out and rename relevant metadata
metadata <- pData(gse)
metadata <- dplyr::select(metadata, description, supplementary_file)
# Extract strain (the word before "mice")
metadata$mouse_strain <- sub(".* (\\S+) mice,.*", "\\1", metadata$description)
metadata$condition <- sub(".*,\\s*", "", metadata$description)
metadata$description <- paste(metadata$mouse_strain, metadata$condition, sep = "-")
metadata$sample_name <- basename(metadata$supplementary_file)
metadata$sample_name <- stringr::str_remove(metadata$sample_name, "\\.gz$")
metadata$supplementary_file <- NULL
rownames(metadata) <- metadata$sample_name

metadata <- metadata %>%
    mutate(cond = ifelse(
        str_count(condition, " ") >= 2,
        str_split_fixed(condition, " ", 3)[, 3],
        condition
    ))
metadata$condition <- NULL
# metadata <- subset(metadata, mouse_strain == "B6J")

write_csv(metadata, path = sprintf("%s/metadata.csv", path))
# metadata <- metadata %>% separate(condition, into = c("time", "unit", "diet"), sep = ", ")
# metadata <- metadata %>% dplyr::select(-sample, -tissue) # gsub("\\)", "", as.character(metadata$cond))

# %% Load and preprocess data
cel_files <- list.celfiles(cel_path, full.names = TRUE)
# cel_files <- cel_files[sapply(strsplit(basename(cel_files), "_"), function(x) x[2]) == "B6J"]
data <- read.celfiles(filenames = cel_files)
eset <- rma(data, target = "core")

# %% running quality control
pData(data)$description <- metadata$description
arrayQualityMetrics(
    expressionset = data,
    outdir = sprintf("%s/array_quality_report", path),
    force = TRUE,
    do.logtransform = FALSE,
    intgroup = "description",
    reporttitle = sprintf("Quality Report %s", accession)
)

# %% visualize normalization levels
pdf(sprintf("%s/boxplot_lfc2.pdf", path), width = 8, height = 6)
boxplot(exprs(eset), outline = FALSE)
dev.off()

# %% i want 2 give up and become a housewife
ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", verbose = TRUE)
annot <- getBM(
    attributes = c("affy_mogene_1_0_st_v1", "ensembl_gene_id", "external_gene_name"),
    filter = "affy_mogene_1_0_st_v1",
    values = rownames(exprs(eset)),
    mart = ensembl
)

# %% Prepare expression data as a data.frame for merging
expr_data <- data.frame(PROBEID = rownames(exprs(eset)), exprs(eset))
full_output <- merge(annot, expr_data, by.x = "affy_mogene_1_0_st_v1", by.y = "PROBEID")
# full_output <- avereps(full_output, full_output$ensembl_gene_id, quote=FALSE)
# full_output <- cbind(annot, exprs(eset))
# annotated <- full_output %>% select(ID, "Gene symbol", "Gene ID", ends_with("CEL.gz"))
write_csv(full_output, path = sprintf("%s/full_B6J.csv", path))
# write_csv(annotated, path = sprintf("%s/annotated_expr.csv", path))
# %%
# pData(eset)$sample_name <- rownames(pData(eset))
# design <- cbind(sample_number = rownames(design), design)
# design$sample_number <- rownames(design)
# design <- merge(pData(eset), design, by = "sample_number")
# rownames(design) <- design$sample_name
# design$sample_name <- NULL
# %%
metadata_subset <- subset(metadata, mouse_strain == "B6J")
exprs_subset <- expr_data[, rownames(metadata_subset)]

# Re-factor disease group
metadata_subset$cond <- factor(metadata_subset$cond)

# Build design matrix
design_subset <- model.matrix(~ 0 + cond, data = metadata_subset)
colnames(design_subset) <- levels(metadata_subset$cond)

# Fit model
fit_subset <- lmFit(exprs_subset, design_subset)

# Define contrasts (e.g. early vs control)
contrasts <- makeContrasts(HFD - control, levels = design_subset)
fit <- lmFit(exprs_subset, design_subset)

fit_subset2 <- contrasts.fit(fit, contrasts)
fit_subset2 <- eBayes(fit_subset2)

full_results <- topTable(fit_subset2, number = Inf)
full_results <- tibble::rownames_to_column(full_results, "ID")
# p_cutoff <- 0.05

# full_results <- filter(full_results, adj.P.Val < p_cutoff)
full_results <- merge(annot, full_results, by.x = "affy_mogene_1_0_st_v1", by.y = "ID")
write_csv(full_results, path = sprintf("%s/B6J_results.csv", path))

# if there are duplicate, use the one with the lowest p-value
summarized <- full_results %>%
    group_by(ensembl_gene_id) %>%
    slice_min(order_by = adj.P.Val, with_ties = FALSE) %>%
    ungroup()

write_csv(summarized, path = sprintf("%s/B6J_results_filtered.csv", path))
# %%
# %%
# %%
# %%
# %%
# %%
# %%
# %%
# %%
# %%
# %%
# %%
# %%
# %%
# %%
# %%
metadata_subset <- subset(metadata, mouse_strain == c("B6J", "B6N"))
exprs_subset <- expr_data[, rownames(metadata_subset)]

# Re-factor disease group
metadata_subset$cond <- factor(metadata_subset$cond, levels = c("control", "HFD"))

# Build design matrix
design_subset <- model.matrix(~ mouse_strain + cond, data = metadata_subset)
# colnames(design_subset) <- levels(metadata_subset$cond)

# Fit model
fit_subset <- lmFit(exprs_subset, design_subset)

# Define contrasts (e.g. early vs control)
contrasts <- makeContrasts(condHFD, levels = design_subset)
fit <- lmFit(exprs_subset, design_subset)

fit_subset2 <- contrasts.fit(fit, contrasts)
fit_subset2 <- eBayes(fit_subset2)

full_results <- topTable(fit_subset2, coef = "condHFD", number = Inf)
full_results <- tibble::rownames_to_column(full_results, "ID")
# p_cutoff <- 0.05

# full_results <- filter(full_results, adj.P.Val < p_cutoff)
full_results <- merge(annot, full_results, by.x = "affy_mogene_1_0_st_v1", by.y = "ID")
write_csv(full_results, path = sprintf("%s/B6X_covar_results.csv", path))

# if there are duplicate, use the one with the lowest p-value
summarized <- full_results %>%
    group_by(ensembl_gene_id) %>%
    slice_min(order_by = adj.P.Val, with_ties = FALSE) %>%
    ungroup()

write_csv(summarized, path = sprintf("%s/B6X_covar_results_filtered.csv", path))
