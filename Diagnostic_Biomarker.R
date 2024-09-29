
# Load the required libraries
library(TCGAbiolinks)
library(dplyr)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(tidyr)
library(stringr)
library(tibble)
library(limma)
library(ggVennDiagram)
library(AnnotationDbi)
library(biomaRt)
library(clusterProfiler)  # For enrichment analysis
library(org.Hs.eg.db)     # Human genome wide annotation database
library(enrichplot)       # For visualization of enrichment results


### Function to combine TCGA sample ID with clinical data
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
  
  # Rename columns
  names(result) <- c("SampID", "gender", "days_to_birth")
  
  return(result)
}



### Function for processing TCGA metadata
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



### Function to impute the missing values of categorical variables with the most frequent category.
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



### Function to align metadata(Sample Info) and log expression data
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



### Function to transform the syntax in the covariate data
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



### Function to set thresholds and filter results (adjusted p-value < 0.01, |log2FC| > 1)
filter_de_genes <- function(results, padj_threshold = 0.01, lfc_threshold = 1) {
  return(results[results$adj.P.Val < padj_threshold & abs(results$logFC) > lfc_threshold, ])
}



### Function to analyze consistency and rank biomarkers
# Takes a list of biomarker genes and their differential expression results from two comparisons
rank_biomarkers <- function(biomarkers, de_results1, de_results2) {
  # Convert row names to a column named "gene_id" in both DE result datasets
  de_results1 <- de_results1 %>% 
    rownames_to_column(var = "gene_id")
  de_results2 <- de_results2 %>% 
    rownames_to_column(var = "gene_id")
  
  biomarkers %>%
    # Convert biomarkers vector to a dataframe with a gene_id column
    enframe(name = NULL, value = "gene_id") %>%
    # Join with DE results
    left_join(de_results1, by = "gene_id") %>%
    left_join(de_results2, by = "gene_id", suffix = c("_1", "_2")) %>%
    # Filter for consistent direction of change
    filter(sign(logFC_1) == sign(logFC_2)) %>%
    # Calculate ranking metrics to assess the magnitude and consistency of their expression changes
    mutate(
      avg_log2FC = (abs(logFC_1) + abs(logFC_2)) / 2,
      consistency_score = pmin(abs(logFC_1), abs(logFC_2)) / 
        pmax(abs(logFC_1), abs(logFC_2)),
      direction = ifelse(logFC_1 > 0, "Up", "Down")
    ) %>%
    # Rank biomarkers
    arrange(desc(avg_log2FC * consistency_score))
}



### Function to prepare a ranked list of all genes for GSEA
prepare_gsea_input <- function(de_results, ranking_column = "logFC", id_column = NULL) {
  # If id_column is not provided, assume Ensembl IDs are row names
  if (is.null(id_column)) {
    ensembl_ids <- rownames(de_results)
  } else {
    ensembl_ids <- de_results[[id_column]]
  }
  
  # Remove version numbers from Ensembl IDs if present
  ensembl_ids <- sub("\\.[0-9]+$", "", ensembl_ids)
  
  # Set up biomaRt
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Convert Ensembl IDs to Entrez IDs
  id_map <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                  filters = "ensembl_gene_id",
                  values = ensembl_ids,
                  mart = mart)
  
  # Merge with original data
  de_results$ensembl_gene_id <- ensembl_ids
  de_results <- merge(de_results, id_map, by = "ensembl_gene_id", all.x = TRUE)
  
  # Create ranked gene list
  ranked_genes <- de_results[[ranking_column]]
  names(ranked_genes) <- de_results$entrezgene_id
  
  # Remove NA values and duplicates
  ranked_genes <- ranked_genes[!is.na(names(ranked_genes))]
  ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]
  
  # Sort the vector in descending order
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # Return both the updated data frame and the ranked gene list
  return(list(
    updated_de_results = de_results,
    ranked_genes = ranked_genes
  ))
}




### Function to perform GSEA
perform_gsea <- function(ranked_genes, name, pvalue_cutoff = 0.05) {
  # Perform GO enrichment analysis using GSEA
  go_gsea <- tryCatch({
    gseGO(geneList = ranked_genes,
          OrgDb = org.Hs.eg.db,
          ont = "ALL",
          minGSSize = 10,
          maxGSSize = 500,
          pvalueCutoff = pvalue_cutoff,
          verbose = FALSE)
  }, error = function(e) {
    message("Error in GO enrichment: ", e$message)
    return(NULL)
  })
  
  # Perform KEGG pathway analysis using GSEA
  kegg_gsea <- tryCatch({
    gseKEGG(geneList = ranked_genes,
            organism = 'hsa',
            minGSSize = 10,
            maxGSSize = 500,
            pvalueCutoff = pvalue_cutoff,
            verbose = FALSE)
  }, error = function(e) {
    message("Error in KEGG enrichment: ", e$message)
    return(NULL)
  })
  
  # Plot and save top GO term if results are not empty
  if (!is.null(go_gsea) && nrow(go_gsea@result) > 0) {
    png(paste0("GO_GSEA_plot_", name, ".png"), width = 800, height = 600)
    print(gseaplot2(go_gsea, geneSetID = 1, title = go_gsea$Description[1]))
    dev.off()
  } else {
    message("No enriched GO terms found at p-value cutoff ", pvalue_cutoff)
  }
  
  # Plot and save top KEGG pathway if results are not empty
  if (!is.null(kegg_gsea) && nrow(kegg_gsea@result) > 0) {
    png(paste0("KEGG_GSEA_plot_", name, ".png"), width = 800, height = 600)
    print(gseaplot2(kegg_gsea, geneSetID = 1, title = kegg_gsea$Description[1]))
    dev.off()
  } else {
    message("No enriched KEGG pathways found at p-value cutoff ", pvalue_cutoff)
  }
  
  return(list(go_gsea = go_gsea, kegg_gsea = kegg_gsea))
}



# Function to map "tobacco_smoking_history" from clinical data to covariate data used for DE analysis
map_smoking_history <- function(A, B) {
  # Create a new column in A for the extracted prefix of "SampID"
  A$SampID_prefix <- substr(A$SampID, 1, 12)
  
  # Now perform a left join using "SampID_prefix" in A and "bcr_patient_barcode" in B
  A <- merge(A, B[, c("bcr_patient_barcode", "tobacco_smoking_history")], 
             by.x = "SampID_prefix", by.y = "bcr_patient_barcode", 
             all.x = TRUE)
  
  # Drop the "SampID_prefix" and "Batch" column as they're no longer needed
  A$SampID_prefix <- NULL
  A$Batch <- NULL
  
  # Return the updated DataFrame A with the mapped column
  return(A)
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

# Filter out lowly expressed genes (Keep genes with TPM > 0.5 in at least 10% of the samples)
keep_luad <- rowSums(gtex_luad > 0.5) >= ncol(gtex_luad) * 0.1
clean_gtex_luad <- gtex_luad[keep_luad, ]

keep_lusc <- rowSums(gtex_lusc > 0.5) >= ncol(gtex_lusc) * 0.1
clean_gtex_lusc <- gtex_lusc[keep_lusc, ]

keep_nsclc <- rowSums(gtex_nsclc > 0.5) >= ncol(gtex_nsclc) * 0.1
clean_gtex_nsclc <- gtex_nsclc[keep_nsclc, ]

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
# Filter results with thresholds
# ==================================================================

## Use filter_de_genes function for each result set to get filtered results
# LUAD
de_luad <- filter_de_genes(results_luad)
# LUSC
de_lusc <- filter_de_genes(results_lusc)
# NSCLC
de_LUAD_vs_Healthy <- filter_de_genes(results_LUAD_vs_Healthy)
de_LUSC_vs_Healthy <- filter_de_genes(results_LUSC_vs_Healthy)
de_LUAD_vs_LUSC <- filter_de_genes(results_LUAD_vs_LUSC)

# =======================================================
# Find overlapping DE genes for LUAD and LUSC separately
# =======================================================

# Find the intersection of LUAD vs GTEx analysis and the contrast in combined NSCLC analysis
luad_biomarkers <- intersect(rownames(de_luad), rownames(de_LUAD_vs_Healthy))
# Find the intersection of LUSC vs GTEx analysis and the contrast in combined NSCLC analysis
lusc_biomarkers <- intersect(rownames(de_lusc), rownames(de_LUSC_vs_Healthy))

# Visualize the overlap using Venn diagrams
ggVennDiagram(
  x = list(LUAD_vs_GTEx = rownames(de_luad), LUAD_in_Combined = rownames(de_LUAD_vs_Healthy)),
  filename = "luad_biomarkers_venn.png",
  main = "LUAD Biomarkers"
) + scale_fill_gradient(low="grey90",high = "blue")

ggVennDiagram(
  x = list(LUSC_vs_GTEx = rownames(de_lusc), LUSC_in_Combined = rownames(de_LUSC_vs_Healthy)),
  filename = "lusc_biomarkers_venn.png",
  main = "LUSC Biomarkers"
) + scale_fill_gradient(low="grey90",high = "blue")

# ====================================
# Identify potential NSCLC biomarkers
# ====================================

# Find the intersection of the biomarkers of LUAD and LUSC
nsclc_biomarkers <- intersect(luad_biomarkers, lusc_biomarkers)

# Visualize NSCLC biomarkers
ggVennDiagram(
  x = list(LUAD_Biomarkers = luad_biomarkers, LUSC_Biomarkers = lusc_biomarkers),
  filename = "nsclc_biomarkers_venn.png",
  main = "NSCLC Biomarkers"
) + scale_fill_gradient(low="grey90",high = "blue")

# =====================================================
# Filter for genes with consistent expression patterns
# =====================================================

## Use rank_biomarkers function to analyze consistency and rank biomarkers
# Rank LUAD biomarkers
ranked_luad_biomarkers <- rank_biomarkers(luad_biomarkers, results_luad, 
                                          results_LUAD_vs_Healthy)

# Rank LUSC biomarkers
ranked_lusc_biomarkers <- rank_biomarkers(lusc_biomarkers, results_lusc, 
                                          results_LUSC_vs_Healthy)

# Rank NSCLC biomarkers
ranked_nsclc_biomarkers <- rank_biomarkers(nsclc_biomarkers, results_LUAD_vs_Healthy, 
                                           results_LUSC_vs_Healthy)

# Select the first 100 items from "gene_id"
biomarker_subset_nsclc <- ranked_nsclc_biomarkers$gene_id[1:100]

# Save the subset to a file
write.csv(biomarker_subset_nsclc, file = "NSCLC_biomarkers_de_analysis.csv", row.names = FALSE)

# ==========================
# Data preparation for GSEA
# ==========================

# Prepare data for GSEA
gsea_input_luad <- prepare_gsea_input(results_luad)
gsea_input_lusc <- prepare_gsea_input(results_lusc)
gsea_input_LUAD_vs_Healthy <- prepare_gsea_input(results_LUAD_vs_Healthy)
gsea_input_LUSC_vs_Healthy <- prepare_gsea_input(results_LUSC_vs_Healthy)

# Access the updated data frames and ranked gene lists
results_luad_updated <- gsea_input_luad$updated_de_results
ranked_gsea_luad <- gsea_input_luad$ranked_genes

results_lusc_updated <- gsea_input_lusc$updated_de_results
ranked_gsea_lusc <- gsea_input_lusc$ranked_genes

results_LUAD_vs_Healthy_updated <- gsea_input_LUAD_vs_Healthy$updated_de_results
ranked_gsea_LUAD_vs_Healthy <- gsea_input_LUAD_vs_Healthy$ranked_genes

results_LUSC_vs_Healthy_updated <- gsea_input_LUSC_vs_Healthy$updated_de_results
ranked_gsea_LUSC_vs_Healthy <- gsea_input_LUSC_vs_Healthy$ranked_genes

# ===========================================================
# Perform GO enrichment and KEGG pathway analysis using GSEA
# ===========================================================

# Perform GSEA for each dataset and get plots
gsea_results_luad <- perform_gsea(ranked_gsea_luad, "luad", pvalue_cutoff = 0.1)
gsea_results_lusc <- perform_gsea(ranked_gsea_lusc, "lusc", pvalue_cutoff = 0.1)
gsea_results_LUAD_vs_Healthy <- perform_gsea(ranked_gsea_LUAD_vs_Healthy, "LUAD_vs_Healthy", pvalue_cutoff = 0.1)
gsea_results_LUSC_vs_Healthy <- perform_gsea(ranked_gsea_LUSC_vs_Healthy, "LUSC_vs_Healthy", pvalue_cutoff = 0.1)

# ====================================
# Data preparation for model training
# ====================================

## LUAD and LUAS covariate data
# Add "tabacco_smaking history" column to the LUAD and LUSC covariate data
clinical_ml_luad <- map_smoking_history(luad_cov_processed, clinical_luad)
# Use mode_imputation function to fill in missing values
clinical_ml_luad <- mode_imputation(clinical_ml_luad, "tobacco_smoking_history")
clinical_ml_lusc <- map_smoking_history(lusc_cov_processed, clinical_lusc)
clinical_ml_lusc <- mode_imputation(clinical_ml_lusc, "tobacco_smoking_history")

## GTEx covariate data
# Create a copy of GTEx covariate data
gtex_ml_cov <- gtex_cov
# Remove the "Batch" column from the copied data frame
gtex_ml_cov <- gtex_ml_cov[, !colnames(gtex_ml_cov) %in% "Batch"]
# Add a "Condition" column
gtex_ml_cov <- gtex_ml_cov %>%
  mutate(tobacco_smoking_history = "0") %>%
  select(1, Condition, everything())
# Use mode_imputation function to fill in missing values
gtex_ml_cov <- mode_imputation(gtex_ml_cov, "tobacco_smoking_history")

## Combine three clinical data
ml_clinical_data <- rbind(gtex_ml_cov, clinical_ml_luad, clinical_ml_lusc)

# Transpose the combined expression data frame
transposed_log_ml_nsclc <- as.data.frame(t(log_gtex_nsclc))

# Convert row names to a column called "RowNames"
transposed_log_ml_nsclc <- rownames_to_column(transposed_log_ml_nsclc, var = "SampID")

## Align covariate data and the log expression data for the same order and sample IDs
# Sort the row names of expression data and clinical data
sorted_A_names <- sort(transposed_log_ml_nsclc$SampID)
sorted_B_names <- sort(ml_clinical_data$SampID)
    
# Find the intersection of the sorted names
common_names <- intersect(sorted_A_names, sorted_B_names)
    
# Subset dataframes to keep only the rows where SampID is in common_names
log_tpm_transposed <- transposed_log_ml_nsclc[transposed_log_ml_nsclc$SampID %in% common_names, , drop = FALSE]
clinical_cov_nsclc <- ml_clinical_data[ml_clinical_data$SampID %in% common_names, , drop = FALSE]
    
# Reorder dataframes to match the order of common_names
log_tpm_transposed <- log_tpm_transposed[match(common_names, log_tpm_transposed$SampID), , drop = FALSE]
clinical_cov_nsclc <- clinical_cov_nsclc[match(common_names, clinical_cov_nsclc$SampID), , drop = FALSE]

# Save dataframes as CSV files
write.csv(log_tpm_transposed, file = "NSCLC_expression_model_training.csv", row.names = FALSE)
write.csv(clinical_cov_nsclc, file = "NSCLC_labels_model_training.csv", row.names = FALSE)
