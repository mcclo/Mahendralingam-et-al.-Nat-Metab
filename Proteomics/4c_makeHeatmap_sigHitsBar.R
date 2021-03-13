#' This script was used to make a heatmap of the metabolic proteome with an annotation bar ..
#' .. for genes showing the metabolic cell lineage signatures - after ANOVA and Tukey test (p<0.05) and logFC > 0
#' 
#' Note: significant_hits object holds the final signatures
#' BC = basal cell, ML = luminal mature, LP = luminal progenitor
#' 
#' Input: 
#' - gene lists
#' Output:
#' - heat maps - individual for each signature, final heat map with "cluster assignment" or "signature" bar
#' - RData file
#' - heat map

# Call required libraries into environment
library(pheatmap)
library(RColorBrewer)
library(cluster)
library(scales)

# Load files from previous scripts
load("signatures.RData")
load("h.possemato_matrices_and_dendograms.RData")
load("h.pint.patient.pheno.info.RData")

# Read in scripts with plotting functions
source("make_human_heatmap_4c.R") # with "signature" annotation bar
source("make_human_heatmap.R")

# Define function to compute and return matrix with z scores (for heatmap plotting)
compute_z_scores <- function(matrix){
  # Formula for z score: 
  #         [element x - mean of row x is in]/ standard deviation of row x is in
  return (apply(matrix, 1, function(x) (x - mean(x)) / sd(x)))
}

# Make z score matrix
exp_matrix <- h.pint.combat.pos[,1:29]  # metabolic protein expression matrix
exp_matrix_z <- t(compute_z_scores(exp_matrix))

# Make a new vector that indicates which genes belong to which cell type
significant_hits_bar <- ifelse(h.pint.combat.pos$Gene.names %in% significant_hits$BC, "BC",
                                      ifelse(h.pint.combat.pos$Gene.names %in% significant_hits$ML, "ML",
                                             ifelse(h.pint.combat.pos$Gene.names %in% significant_hits$LP, "LP", "Unassigned")))

# Define cell type order
celltype_gene_assignments <- c("BC", "ML", "LP", "Unassigned")

# Plot final heatmap with signature/cluster vector----------------------------------
filename <- sprintf("%s_FINAL_possemato_heat_map_signatures.pdf", format(Sys.Date(), "%Y%m%d"))
final_heatmap_plot <- plot_h_heatmap2(exp_matrix = exp_matrix_z, hc_samples = hclust(possemato.sample_dist), 
                                      hc_proteins = hclust(possemato.protein_dist), filename = filename, title = NA,
                                      pint.pheno=pint.pheno, cluster_vector = significant_hits_bar)

# a) Make heatmap for each 
# Define function to make heat maps for each o the 3 cell signatures-----------------
make_cell_signature_heatmap <- function(cell_signatures, title){
  # cell_signatures is a list of 3 signatures for basal, ML and LP respectively
  plot_titles <- c("BC_population", "ML_population", "LP_population") #names of titles/file names
  celltypes <- c("BC", "ML", "LP")
  sapply(1:3, function(i){
    title <- paste(plot_titles[[i]], title, sep="_")
    signature <- cell_signatures[[i]]
    x <- make_human_heatmap(exp_matrix, signature, title)
    # dev.off()
  })
}
# Run function for each signature list (before, after ANOVA, and after ANOVA+Tukey)
make_cell_signature_heatmap(significant_hits, "significant_hits") #unclustered #cluster_cols = F