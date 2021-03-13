#' This script was used to make a heat map from the protein expression matrix of mammary epithelial markers.
#' 
#' All protein abundance values are LFQ-adjusted IBAQ values.	
#' The 0s in the data were imputed with a random number between 1 and 1.5.	

# Load libraries
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(cluster)

# Load total proteome
load("h.total_proteome_mathepan.RData")

# Load functions to subset proteome and make heat map
source("make_human_heatmap.R")

# Define function to read a list from a text file, delimited by space
get_list_from_txt_file <- function(filename){
  df <- read.table(filename, sep = " ", stringsAsFactors = F, header = T)
  return(colnames(df))
}

# 1. get list of cell-specific mammary epithelial markers
# Genes of interest are the mammary cell markers from A. Casey JCB 2018 paper, Figures 1G and 4C
markers <- get_list_from_txt_file("GeneLists/mammary_epithelial_markers.txt")

# 2. Subset table to genes of interest
# Filter gene list and reformat
gene_list <- markers %>%
  unique %>%
  strsplit("[.]") %>%
  unlist %>%
  .[. != "X"] %>%
  toupper

# 3. Generate heat map 
filename <- "diana_mammary_Casey_markers_combat_human_heatmap.pdf"
title <- "Mammary epithelial cell markers"
mamm_markers_heatmap_plot <- make_human_heatmap(full_proteome = h.pint.combat_matrix, signatures = gene_list,filename = filename, title = title)