#' This script was used to make mammary epithelial cell (MEC) lineage-specific proteomic signatures by:
#' filtering each protein based on significance (p<0.05) based on Tukey's test and fold change (FC > 1) thresholds,
#' rather than cutting dendogram
#' 
#' Note: BC = basal cell, ML = luminal mature, LP = luminal progenitor

# Load required libraries
library(openxlsx) # for writing data

# Load functions and data
load("h.possemato_proteome_matrix.RData")
source("makeSignatures_functions.R")

# Analysis ----
# Prepare expression matrices
exp_matrix <- h.pint.combat.pos[,1:29]

# Average columns according to cell type
h.pint.combat.pos.avg <- average_hmpint_columns(h.pint.combat.pos) #make matrix with average
exp_matrix_avg <- h.pint.combat.pos.avg

# Run ANOVA and Tukey's test on all proteins
significance_df <- compute_ANOVA_tukey_pvals(exp_matrix)
rownames(significance_df) <- rownames(h.pint.combat.pos)

# Calculate log 2 FC of each cell type vs the other
log2FC_df <- compute_logFC_3groups(exp_matrix_avg)
rownames(log2FC_df)<- rownames(h.pint.combat.pos)

# Get data frames for separate analyses
celltypes <- c("BC", "ML", "LP")
subset_dfs <- lapply(celltypes, function(ID){
  subset_df(significance_df, log2FC_df, ID)
})
names(subset_dfs) <- celltypes

#  Assign cell type to row "significant hits"
# Define thresholds
pval_thres <- 0.05
FC <- 1
log2FC_thres_up <- log2(FC) # as data has been log2-transformed

# In a list, for each celltype, make a data frame of significance values and logFC 
sig_dfs <- lapply(1:3, function(i){
  # add column for significance
  add_sig_columns(subset_dfs[[i]], pval_thres, log2FC_thres_up)
})
names(sig_dfs) <- celltypes

# Open file connections
excel_filename <- sprintf("20190703_log2FC_ANOVA_Tukey_tables_thres%s.xlsx", log2FC_thres_up)
wb <- openxlsx::createWorkbook("signatures")

# Initialize list object to hold signature vectors for each cell type
significant_hits <- vector("list", 3)
names(significant_hits) <- celltypes

# For each celltype # TODO do in lapply
for(i in 1:3){
  # Get data frame with sig and FC values
  df <- sig_dfs[[i]]
  celltype <- names(sig_dfs[i])
  
  # Write important data frames and lists to Excel file
  current_sheet <- sprintf("%s Significance Table", celltype)
  addWorksheet(wb, current_sheet)
  print(current_sheet)
 
  # Get information
  description1 <- sprintf("conditions: p_val < %s, log2FC >= %s", pval_thres, log2FC_thres_up)
  description2 <- sprintf("%s proteins meet condition for %s", sum(df$Both_Conditions_Met == "Yes"), celltype)
  
  # Make a data frame of significant hits (ie meet sig and FC thresholds/conditions)
  hits <- data.frame(SignificantHits = df$Gene[df$Both_Conditions_Met == "Yes"])
  
  # Write to same Excel worksheet
  writeData(wb, sheet = current_sheet, x = description1, startRow = 1)
  writeData(wb, sheet = current_sheet, x = description2, startRow = 2)
  writeData(wb, sheet = current_sheet, x = df, startRow = 4)
  writeData(wb, sheet = current_sheet, x = hits, startRow = 4, startCol = ncol(df)+5)
  
  # Add significant hits for celltype to list element
  significant_hits[[i]] <- as.character(hits$SignificantHits)
}
# Sort gene list alphabetically
significant_hits <- lapply(significant_hits, FUN = sort)

# File 1 for our reference
# Write to Excel file
# Hits
write_list_to_Excel("Signatures", my_list = significant_hits, wb)

# p-values and logFCs
addWorksheet(wb, "All Tukey and log2FC Results")
writeData(wb, sheet = "All Tukey and log2FC Results", x = cbind(significance_df, log2FC_df), rowNames = FALSE)

# Combat batch-corrected proteome
addWorksheet(wb, "Metabolic Proteome Combat")
writeData(wb, sheet = "Metabolic Proteome Combat", x = h.pint.combat.pos, rowNames = TRUE)

# "Clean" signatures for Enrichr
signatures <- lapply(significant_hits, function(x){
  x <- gsub(";.*", "", x) #remove any alias after semicolons
  x <- gsub("\\..*","",x) #remove any period and number appended to make gene name unique
  x <- unique(x)
  return(x)
})
write_list_to_Excel("Clean Signatures", signatures, wb)

#Save to excel
saveWorkbook(wb, file = excel_filename, overwrite = TRUE)

# File 2 for Signatures files needed for Luis' ssGSEA analysis (enrichment of signatures across Metabric breast cancer subtypes)
excel_filename2 <- sprintf("3-signatures_FC%s.xlsx", FC)
cell_names <- c("Basal", "Luminal Mature", "Luminal Progenitor")
significant_hits_named <- signatures
names(significant_hits_named) <- cell_names

# Create workbook and write lists to it
wb2 <- createWorkbook()
write_list_to_Excel("Signatures", my_list = significant_hits_named, wb2)
saveWorkbook(wb2, file = excel_filename2, overwrite = TRUE)

# save RData for making heatmap
RData_filename <- sprintf("metabolic.signatures.FC%s.pval%s_mathepan.RData", 2^log2FC_thres_up, pval_thres)
save(significant_hits, log2FC_thres_up, pval_thres, file="signatures.RData")
save(signatures, log2FC_thres_up, pval_thres,  file = "signatures.no.extra.characters.RData")
