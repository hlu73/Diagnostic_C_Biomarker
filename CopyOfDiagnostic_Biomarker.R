
# Load the required libraries
library(TCGAbiolinks)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(stringr)
library(limma)



# Function to combine TCGA sample ID with clinical data
merge_data_frames <- function(df_A, df_B) {
  # Function to extract the matching prefix from column names of TPM data
  extract_prefix <- function(x) {
    str_extract(x, "TCGA-\\d{2}-\\d{4}")
  }
  
  # Create the new data frame
  result <- data.frame(
    ID = colnames(df_A),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      # Create a temporary "matching_id" column that contains the extracted prefix
      matching_id = extract_prefix(ID)
    ) %>%
    left_join(
      df_B %>% select(bcr_patient_barcode, gender, days_to_birth),
      by = c("matching_id" = "bcr_patient_barcode")
    ) %>%
    # Remove the temporary matching column
    select(-matching_id)  
  
  # Rename columns if needed
  names(result) <- c("SampID", "gender", "days_to_birth")
  
  return(result)
}



# Function for processing TCGA metadata
process_dataframe <- function(df) {
  df %>%
    # Add Batch column
    mutate(
      Batch = str_extract(SampID, "(?<=-)([A-Za-z0-9]+)(?=-\\d{2}$)"),
      .after = SampID
    ) %>%
    # Calculate and group age
    mutate(
      Age = case_when(
        abs(days_to_birth) < 40*365.25 ~ "< 40",
        abs(days_to_birth) < 50*365.25 ~ "40-49",
        abs(days_to_birth) < 60*365.25 ~ "50-59",
        abs(days_to_birth) < 70*365.25 ~ "60-69",
        abs(days_to_birth) < 80*365.25 ~ "70-79",
        TRUE ~ "≥ 80"
      )
    ) %>%
    # Transform gender
    mutate(
      gender = case_when(
        gender == "MALE" ~ "1",
        gender == "FEMALE" ~ "2",
        TRUE ~ gender  # Keep original value if not MALE or FEMALE
      )
    ) %>%
  # Remove the days_to_birth column
  select(-days_to_birth)
}



# Function to impute the missing values of categorical variables with the most frequent category.
mode_imputation <- function(df, column_name) {
  # Check if the column exists in the dataframe
  if (!column_name %in% names(df)) {
    stop(paste("Column", column_name, "not found in the dataframe"))
  }
  
  # Find the mode (most frequent category)
  mode_value <- df %>%
    filter(!is.na(!!sym(column_name))) %>%
    count(!!sym(column_name)) %>%
    slice_max(n) %>%
    pull(!!sym(column_name)) %>%
    first()
  
  # Impute missing values with the mode
  df <- df %>%
    mutate(!!sym(column_name) := ifelse(is.na(!!sym(column_name)), mode_value, !!sym(column_name)))
  
  return(df)
}



# Function to align metadata(Sample Info) and log expression data
align_dataframes <- function(dfA, dfB) {
  # Sort the column names of dfA and row names of dfB
  sorted_A_names <- sort(colnames(dfA))
  sorted_B_names <- sort(dfB$SampID)
  
  # Find the intersection of the sorted names
  common_names <- intersect(sorted_A_names, sorted_B_names)
  
  if (length(common_names) == 0) {
    stop("No common names found between dfA columns and dfB rows")
  }
  
  # Subset and reorder dfA based on common_names
  new_dfA <- dfA[, common_names, drop = FALSE]
  
  # Subset dfB to keep only the rows where SampID is in common_names
  new_dfB <- dfB[dfB$SampID %in% common_names, , drop = FALSE]
  
  # Reorder new_dfB to match the order of common_names
  new_dfB <- new_dfB[match(common_names, new_dfB$SampID), , drop = FALSE]
  
  return(list(A = new_dfA, B = new_dfB))
}



# Function to transform the syntax in the covariate data
syntax_transform <- function(df) {
  # Check if the required columns exist
  if (!all(c("Batch", "Age") %in% colnames(df))) {
    stop("The dataframe must contain 'Batch' and 'Age' columns.")
  }
  
  # Remove "-" from Batch column
  df$Batch <- gsub("-", "", df$Batch)
  
  # Transform Age column
  age_mapping <- c(
    "< 40" = 1,
    "40-49" = 2,
    "50-59" = 3,
    "60-69" = 4,
    "70-79" = 5,
    "≥ 80" = 6
  )
  
  df$Age <- age_mapping[df$Age]
  
  # Check if any Age values were not mapped (will be NA)
  if (any(is.na(df$Age))) {
    warning("Some Age values could not be mapped. Check for unexpected values in the Age column.")
  }
  
  return(df)
}



# Function to set thresholds and filter results (adjusted p-value < 0.01, |log2FC| > 1)
filter_de_genes <- function(results, padj_threshold = 0.01, lfc_threshold = 1) {
  return(results[results$adj.P.Val < padj_threshold & abs(results$logFC) > lfc_threshold, ])
}



# ========================================================================
# Obtain NSCLC expression data from TCGA, including both lung adenocarcinoma (LUAD) 
# and lung squamous cell carcinoma (LUSC)).
# ========================================================================

# Query for LUAD RNA-seq data (STAR-Counts)
query_luad <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq",
  sample.type = "Primary Tumor"
)

# Download the LUAD data
GDCdownload(query_luad)

# Prepare the LUAD data (this will include the tpm_unstranded column)
data_luad <- GDCprepare(query_luad)

# Extract the tpm_unstranded data for LUAD
tpm_luad <- as.data.frame(SummarizedExperiment::assay(data_luad, "tpm_unstrand"))

# Extract gene information (such as Ensembl IDs)
gene_ids_luad <- rownames(tpm_luad)

# Assign gene IDs as row names for the LUAD TPM data
rownames(tpm_luad) <- gene_ids_luad

# Query for LUSC RNA-seq data (STAR-Counts)
query_lusc <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq",
  sample.type = "Primary Tumor"
)

# Download the LUSC data
GDCdownload(query_lusc)

# Prepare the LUSC data (this will include the tpm_unstranded column)
data_lusc <- GDCprepare(query_lusc)

# Extract the tpm_unstranded data for LUSC
tpm_lusc <- as.data.frame(SummarizedExperiment::assay(data_lusc, "tpm_unstrand"))

# Extract gene information (such as Ensembl IDs)
gene_ids_lusc <- rownames(tpm_lusc)

# Assign gene IDs as row names for the LUSC TPM data
rownames(tpm_lusc) <- gene_ids_lusc

# View the LUAD TPM dataframe
View(tpm_luad)

# View the LUSC TPM dataframe
View(tpm_lusc)

# ===============================================================================
# Obtain NSCLC clinical data from TCGA, including both lung adenocarcinoma (LUAD) 
# and lung squamous cell carcinoma (LUSC))
# ===============================================================================

# Query for LUAD clinical data
query_clinical_luad <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "bcr xml"
)

# Download LUAD clinical data
GDCdownload(query_clinical_luad)

# Prepare LUAD clinical data
clinical_luad <- GDCprepare_clinic(query_clinical_luad, clinical.info = "patient")

# Query for LUSC clinical data
query_clinical_lusc <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "bcr xml"
)

# Download LUSC clinical data
GDCdownload(query_clinical_lusc)

# Prepare LUSC clinical data
clinical_lusc <- GDCprepare_clinic(query_clinical_lusc, clinical.info = "patient")

# ==========
# GTEx tpm data
# ==========

# Load the GTEx data
gtex_data <- read.csv("GTEx_gene_tpm_lung.csv", header = TRUE)

# View the GTEx data
View(gtex_data)

# ================================================
# Filter biological covariates from clinical data
# ================================================

# Extract only columns of gender and age from TCGA clinical data
luad_cov <- clinical_luad[, c("bcr_patient_barcode", "gender", "days_to_birth")]
lusc_cov <- clinical_lusc[, c("bcr_patient_barcode", "gender", "days_to_birth")]

# Use merge_data_frames function to merge the TCGA sample IDs and the gender and age data
luad_cov_merge <- merge_data_frames(tpm_luad, luad_cov)
lusc_cov_merge <- merge_data_frames(tpm_lusc, lusc_cov)

# Use process_dataframe function to add batch information and transform gender and age data
luad_cov_processed <- process_dataframe(luad_cov_merge)
lusc_cov_processed <- process_dataframe(lusc_cov_merge)

# Use mode_imputation function to fill in missing values
luad_cov_processed <- mode_imputation(luad_cov_processed, "gender")
lusc_cov_processed <- mode_imputation(lusc_cov_processed, "gender")

# Add a "Condition" column
luad_cov_processed <- luad_cov_processed %>%
  mutate(Condition = "LUAD") %>%
  select(1, Condition, everything())

lusc_cov_processed <- lusc_cov_processed %>%
  mutate(Condition = "LUSC") %>%
  select(1, Condition, everything())

# Transform the syntax in "Batch" and "Age"
luad_cov_processed <- syntax_transform(luad_cov_processed)
lusc_cov_processed <- syntax_transform(lusc_cov_processed)

# ====================================================================

# Load GTEx clinical data as dataframes
phenotype_gtex <- read.table("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                             header = TRUE, sep = "\t")
sampAttr_gtex <- read.table("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                            header = TRUE, sep = "\t")

## Extract and combine the information on batch, gender and age
# Extract information on batch, gender and age
clinical_gtex_1 <- phenotype_gtex[, c("SUBJID", "SEX", "AGE")]
clinical_getx_2 <- sampAttr_gtex[, c("SAMPID", "SMNABTCH")]

# Create a new column in clinical_getx_2 with the matching prefix
clinical_getx_2$matching_id <- substr(clinical_getx_2$SAMPID, 1, 10)

# Merge clinical_getx_2 with clinical_getx_1
gtex_cov <- merge(clinical_getx_2, clinical_gtex_1[, c("SUBJID", "SEX", "AGE")], 
                  by.x = "matching_id", 
                  by.y = "SUBJID", 
                  all.x = TRUE)

# Clean up - remove the temporary matching_id column
gtex_cov$matching_id <- NULL

gtex_cov <- gtex_cov %>%
  # Fill blank cells in SMNABTCH column with "BP"
  mutate(SMNABTCH = ifelse(SMNABTCH == "" | is.na(SMNABTCH), "BP", SMNABTCH)) %>%
  # Group AGE
  mutate(AGE = case_when(
    AGE < 40 ~ "< 40",
    AGE >= 40 & AGE < 50 ~ "40-49",
    AGE >= 50 & AGE < 60 ~ "50-59",
    AGE >= 60 & AGE < 70 ~ "60-69",
    AGE >= 70 & AGE < 80 ~ "70-79",
    AGE >= 80 ~ "≥ 80",
    # This catches any NA values or ages outside the expected range
    TRUE ~ "Unknown"  
  )) %>%
  # Rename columns
  rename(
    SampID = SAMPID,
    Batch = SMNABTCH,
    gender = SEX,
    Age = AGE
  )

# Add a "Condition" column with all "healthy"
gtex_cov <- gtex_cov %>%
  mutate(Condition = "healthy") %>%
  select(1, Condition, everything())

# Transform all "-" in sample IDs into "."
gtex_cov <- gtex_cov %>%
  mutate(SampID = str_replace_all(SampID, "-", "."))

# Transform the syntax in "Batch" and "Age"
gtex_cov <- syntax_transform(gtex_cov)

# ==============================
# Check for overlap of gene IDs
# ==============================

# Extract gene IDs from all datasets
gtex_gene_ids <- gtex_data$Name
luad_gene_ids <- rownames(tpm_luad)
lusc_gene_ids <- rownames(tpm_lusc)

# Check for overlap between GTEx and LUAD gene IDs
overlap_gene_ids_getx_tcga <- intersect(gtex_gene_ids, luad_gene_ids)

# Check for overlap between LUSC and LUAD gene IDs
overlap_gene_ids_tcga <- intersect(luad_gene_ids, lusc_gene_ids)

# Display the number of overlapping gene IDs
length(overlap_gene_ids_tcga)
length(overlap_gene_ids_getx_tcga)

# ==========================================
# Filter all datasets for overlapping genes
# ==========================================

# Filter GTEx data to keep only the overlapping gene IDs
filtered_gtex_data <- gtex_data[gtex_data$Name %in% overlap_gene_ids_getx_tcga, ]

# Filter TCGA data to keep only the overlapping gene IDs
filtered_tpm_luad <- tpm_luad[rownames(tpm_luad) %in% overlap_gene_ids_getx_tcga, ]
filtered_tpm_lusc <- tpm_lusc[rownames(tpm_lusc) %in% overlap_gene_ids_getx_tcga, ]

# Verify the number of rows match in both datasets
nrow(filtered_gtex_data)
nrow(filtered_tpm_luad)
nrow(filtered_tpm_lusc)
View(filtered_gtex_data)
View(filtered_tpm_luad)
View(filtered_tpm_lusc)

# ===========================
# Set row names in GTEx data
# ===========================

# Remove row names by setting them to NULL
rownames(filtered_gtex_data) <- NULL
# Convert the gene IDs as row names in GTEx data
filtered_gtex_data <- tibble::column_to_rownames(filtered_gtex_data, var = "Name")

# =========================
# Check for missing values
# =========================

sum(is.na(filtered_gtex_data))
sum(is.na(filtered_tpm_luad))
sum(is.na(filtered_tpm_lusc))

# =============================
# 
# =============================

# 1. Boxplot
tpm_long <- filtered_tpm_lusc %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "TPM")

boxplot <- ggplot(tpm_long, aes(x = sample, y = TPM)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Distribution of log2(TPM+1) values across samples",
       y = "log2(TPM+1)")
print(boxplot)

# 2. PCA plot
pca_data <- prcomp(t(filtered_tpm_lusc[rowMeans(filtered_tpm_lusc) > 0.5, ]), scale. = TRUE)
pca_df <- as.data.frame(pca_data$x)
pca_df$sample <- rownames(pca_df)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = sample)) +
  geom_point() +
  geom_text(hjust = 0, vjust = 0) +
  labs(title = "PCA plot of samples")
print(pca_plot)

# ====================================
# Combine data frames for DE analysis
# ====================================

# Combine GTEx and LUAD data
gtex_luad <- Reduce(function(x, y) merge(x, y, by = "row.names", all = TRUE),
                        list(filtered_gtex_data, filtered_tpm_luad))
# Set gene IDs as row names
rownames(gtex_luad) <- gtex_luad$Row.names
# Remove Row.names column after setting row names
gtex_luad <- gtex_luad[, -1]

# Combine GTEx and LUSC data
gtex_lusc <- Reduce(function(x, y) merge(x, y, by = "row.names", all = TRUE),
                    list(filtered_gtex_data, filtered_tpm_lusc))
# Set gene IDs as row names
rownames(gtex_lusc) <- gtex_lusc$Row.names
# Remove Row.names column after setting row names
gtex_lusc <- gtex_lusc[, -1]

# Combine GTEx and both NSCLC data
gtex_nsclc <- Reduce(function(x, y) {
  merged <- merge(x, y, by = "row.names", all = TRUE)
  rownames(merged) <- merged$Row.names
  merged$Row.names <- NULL
  return(merged)
}, list(filtered_gtex_data, filtered_tpm_luad, filtered_tpm_lusc))

# =================================
# Filter out lowly expressed genes
# =================================

# Filter out lowly expressed genes (e.g., genes with mean TPM < 0.5 across samples)
clean_gtex_luad <- gtex_luad[rowMeans(gtex_luad) > 0.5, ]
clean_gtex_lusc <- gtex_lusc[rowMeans(gtex_lusc) > 0.5, ]
clean_gtex_nsclc <- gtex_nsclc[rowMeans(gtex_nsclc) > 0.5, ]

# ========================
# Log-transform the filtered data to log(TPM+1)
# ========================

log_gtex_luad <- log2(clean_gtex_luad + 1)
log_gtex_lusc <- log2(clean_gtex_lusc + 1)
log_gtex_nsclc <- log2(clean_gtex_nsclc + 1)

# =======================
# Create metadata tables and design matrices
# =======================

## Combine covariate data for each dataset for comparison
# LUAD
gtex_luad_sampInfo <- rbind(gtex_cov, luad_cov_processed)
# LUSC
gtex_lusc_sampInfo <- rbind(gtex_cov, lusc_cov_processed)
# NSCLC
gtex_nsclc_sampInfo <- rbind(gtex_cov, luad_cov_processed, lusc_cov_processed)

## Align covariate data and the log expression data for the same order and sample IDs
# LUAD
gtex_luad_align <- align_dataframes(log_gtex_luad, gtex_luad_sampInfo)
# Extract the aligned dataframes
log_gtex_luad_align <- gtex_luad_align$A
gtex_luad_sampInfo_align <- gtex_luad_align$B
# LUSC
gtex_lusc_align <- align_dataframes(log_gtex_lusc, gtex_lusc_sampInfo)
# Extract the aligned dataframes
log_gtex_lusc_align <- gtex_lusc_align$A
gtex_lusc_sampInfo_align <- gtex_lusc_align$B
# NSCLC
gtex_nsclc_align <- align_dataframes(log_gtex_nsclc, gtex_nsclc_sampInfo)
# Extract the aligned dataframes
log_gtex_nsclc_align <- gtex_nsclc_align$A
gtex_nsclc_sampInfo_align <- gtex_nsclc_align$B

## Constructing the design matrix for limma
# LUAD
design_luad <- model.matrix(~ 0 + Condition + Batch + gender + Age, data = gtex_luad_sampInfo_align)
# LUSC
design_lusc <- model.matrix(~ 0 + Condition + Batch + gender + Age, data = gtex_lusc_sampInfo_align)
# NSCLC
design_nsclc <- model.matrix(~ 0 + Condition + Batch + gender + Age, data = gtex_nsclc_sampInfo_align)

# ========================
# DE analysis using limma
# ========================

# Adjust for unequal variability among samples
weights_luad <- arrayWeights(log_gtex_luad_align, design = design_luad)
weights_lusc <- arrayWeights(log_gtex_lusc_align, design = design_lusc)
weights_nsclc <- arrayWeights(log_gtex_nsclc_align, design = design_nsclc)

## LUAD
# Fit the model
fit_luad <- lmFit(log_gtex_luad_align, design_luad, weights = weights_luad)

# Define contrast
contrast_luad <- makeContrasts(
  LUAD_vs_Healthy = ConditionLUAD - Conditionhealthy,
  levels = design_luad
)

# Apply contrast
fit_luad_2 <- contrasts.fit(fit_luad, contrast_luad)
fit_luad_2 <- eBayes(fit_luad_2, trend=TRUE)

# Get results
results_luad <- topTable(fit_luad_2, coef="LUAD_vs_Healthy", number=Inf)

## LUSC
# Fit the model
fit_lusc <- lmFit(log_gtex_lusc_align, design_lusc, weights = weights_lusc)

# Define contrast
contrast_lusc <- makeContrasts(
  LUSC_vs_Healthy = ConditionLUSC - Conditionhealthy,
  levels = design_lusc
)

# Apply contrast
fit_lusc_2 <- contrasts.fit(fit_lusc, contrast_lusc)
fit_lusc_2 <- eBayes(fit_lusc_2, trend=TRUE)

# Get results
results_lusc <- topTable(fit_lusc_2, coef="LUSC_vs_Healthy", number=Inf)

## NSCLC
# Fit the model
fit_nsclc <- lmFit(log_gtex_nsclc_align, design_nsclc, weights = weights_nsclc)

# Define contrasts
contrast_nsclc <- makeContrasts(
  LUAD_vs_Healthy = ConditionLUAD - Conditionhealthy,
  LUSC_vs_Healthy = ConditionLUSC - Conditionhealthy,
  LUAD_vs_LUSC = ConditionLUAD - ConditionLUSC,
  levels = design_nsclc
)

# Apply contrasts
fit_nsclc_2 <- contrasts.fit(fit_nsclc, contrast_nsclc)
fit_nsclc_2 <- eBayes(fit_nsclc_2)

# Get results
results_LUAD_vs_Healthy <- topTable(fit_nsclc_2, coef="LUAD_vs_Healthy", number=Inf)
results_LUSC_vs_Healthy <- topTable(fit_nsclc_2, coef="LUSC_vs_Healthy", number=Inf)
results_LUAD_vs_LUSC <- topTable(fit_nsclc_2, coef="LUAD_vs_LUSC", number=Inf)

# ==================================================================
# Filter results with thresholds and identify overlapping DE genes
# ==================================================================

## Use filter_de_genes function to gwt filtered results
# LUAD
de_luad <- filter_de_genes(results_luad)
# LUSC
de_lusc <- filter_de_genes(results_lusc)
# NSCLC
de_LUAD_vs_Healthy <- filter_de_genes(results_LUAD_vs_Healthy)
de_LUSC_vs_Healthy <- filter_de_genes(results_LUSC_vs_Healthy)
de_LUAD_vs_LUSC <- filter_de_genes(results_LUAD_vs_LUSC)

# ===

length(intersect(rownames(de_luad), rownames(de_LUAD_vs_Healthy)))
length(intersect(rownames(de_lusc), rownames(de_LUSC_vs_Healthy)))
length(intersect(rownames(de_LUAD_vs_Healthy), rownames(de_LUSC_vs_Healthy)))
length(intersect(intersect(rownames(de_luad), rownames(de_LUAD_vs_Healthy)), 
       intersect(rownames(de_lusc), rownames(de_LUSC_vs_Healthy))))

