#' This script makes "metabolic" proteome by filtering total proteome with metabolic genes.
#' 
#' Metabolic gene list from Supplemental 3 in this paper: 
#' Possemato R, Marks KM, Shaul YD, et al. Functional genomics reveal that the serine synthesis pathway is essential in breast cancer. Nature. 2011;476(7360):346-350. Published 2011 Aug 18. doi:10.1038/nature10350

# Load data
load("h.total_proteome_mathepan.RData") 

# Load package for diana function
library(cluster) # install.packages("cluster")

# Read in helper functions
source("make_human_heatmap.R") # get_present_gene_name(), subset_proteome_table2()

# Read in metabolic/Possemato(= pos) genes
possemato_table <- read.csv("GeneLists/S3 Possemato et al Nature 2011 curated metabolic genes.csv", header = TRUE, stringsAsFactors = FALSE)

# Subset total proteome data to metabolic genes
h.pint.combat.pos <- subset_proteome_table2(
                        table =  h.pint.combat,
                        gene_column = h.pint.combat$Gene.names,
                        signatures = possemato_table$Gene.Symbol
                     )

# Get matrix of the 29 samples
exp_matrix <- as.matrix(h.pint.combat.pos[,1:29])

# Make sample/column distance matrix and dendrogram 
possemato.sample_dist <- as.dist(1-cor(exp_matrix, method = "pearson"))
possemato.sample_dend <- as.dendrogram(diana(possemato.sample_dist)) #DIvisive ANAlysis Clustering

# Make protein/row distance matrix and dendrogram 
possemato.protein_dist <- as.dist(1-cor(t(exp_matrix), method = "pearson"))
possemato.protein_dend <- as.dendrogram(diana(possemato.protein_dist)) #DIvisive ANAlysis Clustering

# Find column indices for each celltype
cell_columns <- list(BC = grep("BC", colnames(h.pint.combat.pos)), 
                     ML = grep("LC", colnames(h.pint.combat.pos)), 
                     LP = grep("LP", colnames(h.pint.combat.pos)))

# Save as RData
save(h.pint.combat.pos, possemato_table, possemato.protein_dend, possemato.protein_dist, possemato.sample_dend, possemato.sample_dist, cell_columns, file ="h.possemato_matrices_and_dendograms.RData")
save(possemato.sample_dend, possemato.sample_dist, file="h.possemato.sample_dendogram.RData")
save(h.pint.combat.pos, cell_columns, file="h.possemato_proteome_matrix.RData")
