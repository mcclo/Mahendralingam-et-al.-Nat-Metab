
library(CytoTRACE)
library(Nebulosa)
library(cowplot)
library(dittoSeq)
library(ggpubr)
library(scClustViz)


#get raw count matrix from your Seurat object
mat <- as.matrix(GetAssayData(b2, slot = "counts"))

#(optional) mat_norm <- as.matrix(NormalizeData(mat))
results <- CytoTRACE(mat, ncores = 1, enableFast = T, subsamplesize = 1000)

#Phenotype (below) refers to a metadata column
celltype <- as.character(b2$celltype)
col <- colnames(b2)
c <- cbind(col, celltype)
rownames(c) <- c[,1]
c <- c[,-1]

#use "emb flag for my own reduction ie/batch correction
# to use my umapharmony embeddings, I need to extract the 2D coordinates as a matrix and supply to the emb flag
umapharmony <- b2@reductions$umapharmony@cell.embeddings
#note: the below two calls save .pdfs into your working directory
plotCytoTRACE(results, phenotype = c, emb = umapharmony)
plotCytoGenes(results, numOfGenes = 20)

#add CytoTRACE scores to your Seurat object
b2 <- AddMetaData(
  object = b2,
  metadata = results$CytoTRACE,
  col.name = "cytoTRACE"
)

b2$celltype <- factor(b2$celltype, levels = c("Basal1", "Basal2","LP1", "LP2", "LP3", "LP4", "ML1", "ML2", "ML3", "ML4", "ML5"))

#to view cytoTRACE scores on umap
w1 <- dittoDimPlot(b2, var = "cytoTRACE", reduction.use = "umapharmony", legend.size = 10, 
                   main = "cytoTRACE score", do.label = F, labels.repel = F, do.ellipse = F, 
                   labels.highlight = T) + 
  theme_classic() +
  theme(legend.text = element_text(size = 25, color = "black", family = "sans"),
        axis.text = element_text(size = 25, color = "black", family = "sans"),
        axis.title = element_text(size = 25, color = "black", family = "sans"),
        title = element_text(size = 25, color = "black", family = "sans"))

w2 <- dittoDimPlot(b2, var = "celltype", reduction.use = "umapharmony", legend.size = 10, 
             main = "Cell states", do.label = F, labels.repel = F, do.ellipse = F,
             labels.highlight = T) + 
  theme_classic() +
  theme(legend.text = element_text(size = 25, color = "black", family = "sans"),
        axis.text = element_text(size = 25, color = "black", family = "sans"),
        axis.title = element_text(size = 25, color = "black", family = "sans"),
        title = element_text(size = 25, color = "black", family = "sans"))

plot_grid(w2, w1, ncol = 2, align = "hvd")

#To run CytoTRACE on each individual cell lineage or even within a cell cluster of interest, you need to 
#subset your object, pull our a new matrix and embeddings, and re-run CytoTRACE.
#BASAL 
#attempt cytotrace on specific populations
Idents(b2) <- "broadcelltype"
basal <- subset(b2, ident = c("Basal"))
mat_basal <- as.matrix(GetAssayData(basal, slot = "counts"))
#"Phenotype refers to a metadata column
celltype_basal <- as.character(basal$celltype)
col_basal <- colnames(basal)
c_basal <- cbind(col_basal, celltype_basal)
rownames(c_basal) <- c_basal[,1]
c_basal <- c_basal[,-1]

#use "emb flag for my own reduction
results_b <- CytoTRACE(mat_basal, ncores = 1, enableFast = T, subsamplesize = 1000)
# to use my umapharmony embeddings, I need to extract the 2D coordinates as a matrix and supply to the emb flag
umapharmony_basal <- basal@reductions$umapharmony@cell.embeddings
plotCytoTRACE(results_b, phenotype = c_basal, gene = "PROCR", emb = umapharmony_basal)
plotCytoGenes(results_b, numOfGenes = 20)

#to plot continuous metadata
plotCytoTRACE(results_b, phenotype = c_basal, otherName = "Hallmark_oxphos", otherValue = basal$HALLMARK_OXIDATIVE_PHOSPHORYLATION, emb = umapharmony_basal)

#to add cytoTRACE results back to seurat metadata
basal <- AddMetaData(
  object = basal,
  metadata = results_b$CytoTRACE,
  col.name = "cytoTRACE"
)


basal$celltype <- factor(basal$celltype, levels = c("Basal1", "Basal2"))


w1 <- dittoDimPlot(basal, var = "celltype", reduction.use = "umapharmony", legend.size = 10, 
                   main = "Cell state", do.label = F, labels.repel = F, do.ellipse = F, size = 2,
                   labels.highlight = T) + 
  theme_classic() +
  theme(legend.text = element_text(size = 15, color = "black", family = "sans"),
        axis.text = element_text(size = 15, color = "black", family = "sans"),
        axis.title = element_text(size = 15, color = "black", family = "sans"),
        title = element_text(size = 15, color = "black", family = "sans"))

w2 <- dittoDimPlot(basal, var = "cytoTRACE", reduction.use = "umapharmony", legend.size = 10, size = 2,
                   main = "CytoTRACE score", do.label = F, labels.repel = F, do.ellipse = F, 
                   labels.highlight = T) + 
  theme_classic() +
  theme(legend.text = element_text(size = 15, color = "black", family = "sans"),
        axis.text = element_text(size = 15, color = "black", family = "sans"),
        axis.title = element_text(size = 15, color = "black", family = "sans"),
        title = element_text(size = 15, color = "black", family = "sans"))

w3 <- dittoDimPlot(basal, var = "HALLMARK_OXIDATIVE_PHOSPHORYLATION", reduction.use = "umapharmony", legend.size = 10, size = 2,
                   main = "OXPHOS", do.label = F, labels.repel = F, do.ellipse = F,
                   labels.highlight = T) + 
  theme_classic() +
  theme(legend.text = element_text(size = 15, color = "black", family = "sans"),
        axis.text = element_text(size = 15, color = "black", family = "sans"),
        axis.title = element_text(size = 15, color = "black", family = "sans"),
        title = element_text(size = 15, color = "black", family = "sans")) 


par(mar = c(2.5, 5,0, 7))
plot_grid(w1, w2, w3, ncol = 1, align = "hvd")

#To correlate CytoTRACE scores with other metadata such as HALLMARK ssGSEA scores

basal_meta <- getMD(basal)
par(mar = c(2.5, 5,0, 7))
ggscatter(basal_meta, x = "cytoTRACE", y = "HALLMARK_OXIDATIVE_PHOSPHORYLATION", # specify data and aesthetics 
          add = "reg.line", add.params = list(color = "red", size = 2), conf.int = TRUE, # add regression to plot
          cor.coef = F, cor.method = "pearson", 
          size = 2, font.label = c(20, "bold", "red"), corr.coef.size = 20,#label correlation coeffient
          xlab = "cytoTRACE", ylab = "Oxidative Phosphorylation ssGSEA score") +
  stat_cor(show.legend = F, size = 6.5) +
  theme(text = element_text(size=16))


#LUMINAL PROGENITOR
LP <- subset(b2, ident = c("Luminal progenitor"))
mat_LP <- as.matrix(GetAssayData(LP, slot = "counts"))
#"Phenotype refers to a metadata column
celltype_LP <- as.character(LP$celltype)
col_LP <- colnames(LP)
c_LP <- cbind(col_LP, celltype_LP)
rownames(c_LP) <- c_LP[,1]
c_LP <- c_LP[,-1]

#use "emb flag for my own reduction
results_LP <- CytoTRACE(mat_LP, ncores = 1, enableFast = T, subsamplesize = 1000)
# to use my umapharmony embeddings, I need to extract the 2D coordinates as a matrix and supply to the emb flag
umapharmony_LP <- LP@reductions$umapharmony@cell.embeddings
plotCytoTRACE(results_LP, phenotype = c_LP, gene = "GATA3", emb = umapharmony_LP)
plotCytoGenes(results_LP, numOfGenes = 20)

LP <- AddMetaData(
  object = LP,
  metadata = results_LP$CytoTRACE,
  col.name = "cytoTRACE"
)

#FOR LP
#paper figure for LP
LP$celltype <- factor(LP$celltype, levels = c("LP1", "LP2", "LP3", "LP4"))


w1 <- dittoDimPlot(LP, var = "celltype", reduction.use = "umapharmony", legend.size = 10, 
                   main = "Cell state", do.label = F, labels.repel = F, do.ellipse = F, size = 2,
                   color.panel = dittoColors(), colors = c(3:7),
                   labels.highlight = T) + 
  theme_classic() +
  theme(legend.text = element_text(size = 15, color = "black", family = "sans"),
        axis.text = element_text(size = 15, color = "black", family = "sans"),
        axis.title = element_text(size = 15, color = "black", family = "sans"),
        title = element_text(size = 15, color = "black", family = "sans"))

w2 <- dittoDimPlot(LP, var = "cytoTRACE", reduction.use = "umapharmony", legend.size = 10, size = 2,
                   main = "CytoTRACE score", do.label = F, labels.repel = F, do.ellipse = F, 
                   labels.highlight = T) + 
  theme_classic() +
  theme(legend.text = element_text(size = 15, color = "black", family = "sans"),
        axis.text = element_text(size = 15, color = "black", family = "sans"),
        axis.title = element_text(size = 15, color = "black", family = "sans"),
        title = element_text(size = 15, color = "black", family = "sans"))

w3 <- dittoDimPlot(LP, var = "HALLMARK_OXIDATIVE_PHOSPHORYLATION", reduction.use = "umapharmony", legend.size = 10, size = 2,
                   main = "OXPHOS", do.label = F, labels.repel = F, do.ellipse = F,
                   labels.highlight = T) + 
  theme_classic() +
  theme(legend.text = element_text(size = 15, color = "black", family = "sans"),
        axis.text = element_text(size = 15, color = "black", family = "sans"),
        axis.title = element_text(size = 15, color = "black", family = "sans"),
        title = element_text(size = 15, color = "black", family = "sans")) 


par(mar = c(2.5, 5,0, 7))
plot_grid(w1, w2, w3, ncol = 1, align = "hvd")

#for correlation analysis

LP_meta <- getMD(LP)
par(mar = c(2.5, 5,0, 7))
l1 <- ggscatter(basal_meta, x = "cytoTRACE", y = "HALLMARK_OXIDATIVE_PHOSPHORYLATION", # specify data and aesthetics 
          add = "reg.line", add.params = list(color = "red", size = 2), conf.int = TRUE, # add regression to plot
          cor.coef = F, cor.method = "pearson", 
          size = 2, font.label = c(20, "bold", "red"), corr.coef.size = 20,#label correlation coeffient
          xlab = "cytoTRACE", ylab = "OXPHOS ssGSEA score", title = "Basal epithelial cells") +
  stat_cor(show.legend = F, size = 6.5) +
  theme(axis.title = element_text(size=25, color = "black", family = "sans"),
        axis.title.x = element_blank(),
        title = element_text(size=25, color = "#ED1C24", family = "sans"))

l2 <- ggscatter(LP_meta, x = "cytoTRACE", y = "HALLMARK_OXIDATIVE_PHOSPHORYLATION", # specify data and aesthetics 
                add = "reg.line", add.params = list(color = "red", size = 2), conf.int = TRUE, # add regression to plot
                cor.coef = F, cor.method = "pearson", 
                size = 2, font.label = c(20, "bold", "red"), corr.coef.size = 20,#label correlation coeffient
                xlab = "cytoTRACE", ylab = "OXPHOS ssGSEA score", title = "Luminal progenitor cells") +
  stat_cor(show.legend = F, size = 6.5) +
  theme(axis.title = element_text(size=25, color = "black", family = "sans"),
        title = element_text(size=25, color = "#00AEEF", family = "sans"))

plot_grid(l1, l2, ncol = 1, align = "hvd")

