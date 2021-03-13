#' This script was used to make PCA plots for the metabolic and total human mammary proteome.
#' 
#' Annotations: hormone statuses are different shapes, cells are different colours
#' Package ggbiplot based on ggplot2 was used to make ellipses around clusters
#' 
#' Note: BC = basal cell, ML = luminal mature, LP = luminal progenitor
#' 
#' A nice explanation of PCA: https://stats.stackexchange.com/questions/2691/making-sense-of-principal-component-analysis-eigenvectors-eigenvalues

# Load libraries
library(svglite)
library(devtools) # install.packages("devtools")
library(ggbiplot) # devtools::install_github("vqv/ggbiplot") 
library(scales)

# Load data
load("h.total_proteome_mathepan.RData") #total and metabolic proteome
load("h.pint.patient.pheno.info.RData") #patient annotations
rm(h.pint.nonimputed.avg)

# In patient info (pint.pheno), replace Luminal Cell (LC) name with "Mature Luminal" (ML) for Mathepan 
pint.pheno$cell_type <- gsub("LC", "ML", pint.pheno$cell_type)
all(colnames(h.pint.combat_matrix) == pint.pheno$sample_name) #should be TRUE

# Prepare PCA dataframe
df_pca <- prcomp(t(h.pint.combat_matrix), scale. = TRUE) #scaling is recommended for untransformed/unscaled data
df_out <- as.data.frame(df_pca$x)

# Add annotations
df_out$CellType <- pint.pheno$cell_type
df_out$HormoneStatus <- pint.pheno$hormone_status
head(df_out)

# Calculate percentages for axes labels
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)

# Format axes titles - note: gsub removes spaces inside string
percentage <- gsub(" ", "", paste(colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="")))

# Make annotations for sample cell type (color) and hormone status of patients (shape)
# ggpubr::show_point_shapes()
cell_colors = c(BC="red", ML="blue", LP="skyblue")#, S="darkgray")
hormone.status_shapes <- c(Follicular=2, Luteal=5, Post.Menopausal=0)

# Rescale data so points are between -25 and 25 on both axes
df_pca$x[,1] <- scales::rescale(df_pca$x[,1], to=c(-25,25))
df_pca$x[,2] <- scales::rescale(df_pca$x[,2], to=c(-25,25))

windows()
# Plot PCA with ellipses around cell type clusters
ggbiplot(df_pca_backup, choice = c(1,2), obs.scale = 1, var.scale = 1, 
         groups = pint.pheno$cell_type, ellipse = TRUE, circle = FALSE, var.axes = FALSE) +
  geom_point(aes(shape = pint.pheno$hormone_status, color = pint.pheno$cell_type), size = 3) +
  scale_shape_manual(values = hormone.status_shapes) +
  scale_color_manual(values=cell_colors) +
  # geom_text(size=1.5) +
  xlab(percentage[1]) + 
  ylab(percentage[2]) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key=element_blank()) 

# Save plot
ggsave("Original_PCA_plot_total_proteome_scaled_withellipses.svg", last_plot()) #_withstroma
