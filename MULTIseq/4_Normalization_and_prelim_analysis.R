#playing around with batch effect correction
library(Seurat)
library(ggplot2)
library(sctransform)
library(affy)
library(dplyr)
library(stringr)
library(pheatmap)
library(umap)
library(harmony)
library(SingleCellExperiment)
library(pheatmap)
library(scater)
library(scran)
library(DESeq2)
library(cellassign)


#load each individual RUN.RData file we each RUN/batch you are adding together. These RUN.RData files are files saved 
# after individual quality-control (demultiplexing, doublet removal) file

#_____Now subset all runs to only have the RM samples
RUN1_SCT_final$Genotype <- plyr::mapvalues(
  x = RUN1_SCT_final$barcode, 
  from = c("sample1", "sample2", "sample3", "sample4"), 
  to = c("RM", "RM", "mutneg", "BRCA1")) 
#check to make sure renaming worked
table(RUN1_SCT_final$Genotype)
Idents(RUN1_SCT_final) <- "Genotype"
RUN1_RM <- subset(RUN1_SCT_final, idents = "RM")
#check
table(RUN1_RM$Genotype)

#Run5
table(RUN5_SCT_final$Genotype)
Idents(RUN5_SCT_final) <- "Genotype"
RUN5_RM <- subset(RUN5_SCT_final, idents = "RM")
#check
table(RUN5_RM$Genotype)

#Run6
table(RUN6_SCT_final$Genotype, RUN6_SCT_final$barcode)
Idents(RUN6_SCT_final) <- "barcode"
RUN6_RM <- subset(RUN6_SCT_final, idents = "25")
Idents(RUN6_RM) <- "Genotype"
#check
table(RUN6_RM$Genotype)
#Now merge all runs without batch correction 
merged_runs <- merge(RUN1_RM, y = c(RUN5_RM, RUN6_RM), add.cell.ids = NULL, project = "Breast")
head(colnames(merged_runs))
win.graph()
VlnPlot(merged_runs, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#make a backup object! could also save this separately. 
merged_back <- merged_runs

FeatureScatter(merged_runs, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(merged_runs, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#I know cluster 5 is from one -patient and has virtually no genes detected...thus increase nFeature_RNA (you have to go through entire code and then come back and rerun with different cutoffs)
merged_runs <- subset(merged_runs, subset = nFeature_RNA > 1000)

merged_runs <- NormalizeData(merged_runs, normalization.method = "LogNormalize", scale.factor = 10000)

merged_runs <- FindVariableFeatures(merged_runs, selection.method = "vst", nfeatures = 2000)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(merged_runs), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(merged_runs)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#score cell cycle after 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

merged_runs <- CellCycleScoring(merged_runs, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

#you can also compute cell cycle earlier and use the CC.Difference below to regress out the difference between S and G2M score.
#merged_runs$CC.Difference <- merged_runs$S.Score - merged_runs$G2M.Score
all.genes <- rownames(merged_runs)

merged_runs <- ScaleData(merged_runs, features = all.genes)

merged_runs <- RunPCA(merged_runs, verbose = F)

plot1 <- DimPlot(merged_runs, reduction = "pca", group.by = "Genotype")
plot2 <- DimPlot(merged_runs, reduction = "pca", group.by = "Phase")
plot1 + plot2

merged_runs <- JackStraw(merged_runs, num.replicate = 100)
merged_runs <- ScoreJackStraw(merged_runs, dims = 1:20)
ElbowPlot(merged_runs)

#Using harmony for batch correction in classic Seurat workflow, NOT SCTransform (you can test variations of batch correction using the CellMixS pacakge)
merged_runs <- RunHarmony(merged_runs, "run", assay.use = "RNA")
merged_runs <- FindNeighbors(merged_runs, reduction = "harmony", dims = 1:30, verbose = FALSE)
merged_runs <- FindClusters(merged_runs, resolution = 0.8, verbose = FALSE)
merged_runs <- RunUMAP(merged_runs, reduction = "harmony", assay = "RNA", dims = 1:30, verbose = FALSE, reduction.name = "umapharmony")

DimPlot(merged_runs, reduction = "umapharmony", label = T, pt.size = 0.5, label.size = 4, group.by = "seurat_clusters", ncol = 3)

win.graph()
FeaturePlot(merged_runs, reduction = "umapharmony", features = c("nFeature_RNA"), cols = c("lightblue1", "blue"), ncol = 3, label = TRUE)


#cell type assignment using CellAssign or SCSA in Linux/python
#now going onto cellassign phase need to first convert my seurat object into singlecellexperiment 

#for CellAssign, DO NOT USE SCTransform 
SCE <- as.SingleCellExperiment(merged_runs, assay = "RNA")
#need to first compute size factors for each cell but can i use this on SCT slot or does it take the raw counts from SCT? 
SCE <-computeSumFactors(SCE)
#store size factors -this step can take some time. 
s <- sizeFactors(SCE)


#the following list is curated from multiple sources within CellMarker database
breast_marker_list <- list(
  'Basal epithelial cell' = c('ACTA2', 'KRT14','KRT5', 'MYLK', 'TP63'),
  'Luminal epithelial cell' = c('KRT8', 'ANKRD30A','SYTL2'),
  'Luminal progenitor cell' = c('KRT8', 'KRT19', 'PROM1','SLPI'),
  'Preadipocyte' = c("CD34", "PDGFRA", "DLK1", "PECAM1"),
  'B cell' = c('BLNK', 'CD19','CD79A', 'CD79B', 'MS4A1'),
  'CD4+ T cell' = c('CD4', 'CTLA4', 'FOXP3', 'IL2RA'),
  'CD8+ T cell' = c('CD8A', 'CD8B','GZMK'),
  'Endothelial cell' = c('CDH5', 'SELE','VWF', 'CD34', 'PECAM1'),
  'Fibroblast' = c('COL1A1', 'COL3A1', 'FAP','THY1'),
  'Macrophage' = c('CD14', 'CD163', 'CD68', 'CSF1R', 'FCGR3A'),
  'Natural killer cell' = c('FCGR3A', 'KLRB1', 'KLRB1','KLRC1', 'KLRD1', 'KLRF1','KLRK1', 'NCAM1'),
  'Pericyte' = c('MCAM', 'CSPG4')
)

#now using one of the methods above, make the marker list into a binary matrix
#convert gene lists to binary matrix for input into cellassign. 
marker_gene_mat <- marker_list_to_mat(breast_marker_list, include_other = FALSE)

stopifnot(all(!is.na(marker_gene_mat)))
#view a heatmap of the binary matrix 

pheatmap(marker_gene_mat, fontsize = 8)

#make sure the marker genes are present/named correctly in my SCE object 

sce_marker <- SCE[intersect(rownames(marker_gene_mat), rownames(SCE)),]

result1 <- setdiff(rownames(marker_gene_mat), rownames(sce_marker))

stopifnot(all(!is.na(sce_marker)))


#not sure about this step
counts(sce_marker) <- as.matrix(counts(sce_marker))

print(dim(sce_marker))

#make sure the dims of the marker matrix are the same in terms of gene number 
print(dim(marker_gene_mat))

#adding batch as a covariate
N <- ncol(SCE)
x1 <- SCE$run

X <-model.matrix(~0 + x1)
#call cellassign 

fit <- cellassign(exprs_obj = sce_marker,
                  marker_gene_info = marker_gene_mat,
                  min_delta = 2,
                  X = X,
                  s = s,
                  learning_rate = 1e-2,
                  shrinkage = TRUE,
                  verbose = FALSE)
#check the fit
fit

#add the cellassign calls to the original SCE object 
SCE$cell_type <- fit$cell_type
plotReducedDim(SCE, dimred = "UMAPHARMONY", colour_by = "cell_type", text_by = "cell_type", text_size = 6)

#bring cell_type calls directly into seurat object
merged_runs$cell_type <- fit$cell_type

Idents(merged_runs) <- 'seurat_clusters'
DimPlot(merged_runs, reduction = "umapharmony", group.by = "seurat_clusters", ncol = 1, label = T)


#for input into SCSA in UBUNTU (running python3)
Idents(merged_runs) <- "seurat_clusters"
mark <- FindAllMarkers(merged_runs)
mark2 <- mark
rownames(mark2) <- NULL
write.table(mark, file = "C:/Users/home/SCSA_celltype/SCSA/mark.csv", col.names = T, row.names = F, sep = ",")

curated <- read.csv("C:/Users/home/SCSA_celltype/SCSA/curated_markers.csv")
write.table(curated, file = "C:/Users/home/SCSA_celltype/SCSA/curated.table", row.names = F, sep = "\t", quote = F)

#Ran SCSA in Ubuntu using the "ALL" tissues and my own input list of cell types/ 
#SCSA takes cluster numbers as input. 
#Can name accordingly:
merged_runs$SCSA <-  plyr::mapvalues(
  x = merged_runs$seurat_clusters, 
  from = as.character(c(seq(0, 18, by = 1))), 
  to = c("Luminal progenitor", "Fibroblast", "Mature luminal", "Vascular endothelial", "Luminal progenitor", "Basal", "Vascular endothelial", "Vascular endothelial"
         ,"Pericyte", "CD8+ T cell/NKT", "Mature luminal", "Fibroblast", "Fibroblast", "B cell", "Macrophage", "Lymphatic endothelial"
         , "Luminal progenitor","Mature luminal", "Vascular endothelial"))

DimPlot(merged_runs, reduction = "umapharmony", label = T, pt.size = 0.5, label.size = 5, group.by = c("SCSA"), ncol = 1) + NoLegend()
DimPlot(s, reduction = "umapharmony", label = T, pt.size = 0.5, label.size = 4, split.by = "Genotype", group.by = "cell_type", ncol = 3)


#subset just basal and luminal cells 
Idents(merged_runs) <- "SCSA"
basal_luminal_list <- WhichCells(merged_runs, idents = c("Luminal progenitor", "Mature luminal", "Basal"))

#-----------------------------------------------------------------------------------------------------------------
#pull put ID and recluster from scratch 

b2 <- subset(merged_back, cells = basal_luminal_list)

b2 <- subset(b2, subset = nFeature_RNA > 1000)

#need to remove the cells that will end up falling into ML5 since they express PTPRC (failed classification methods)
ML5cells <- WhichCells(b2, idents = "ML5")
b2 <- subset(b2, cell = ML5cells, invert = T)

#now on to normalization
b2 <- NormalizeData(b2, normalization.method = "LogNormalize", scale.factor = 10000)

b2 <- FindVariableFeatures(b2, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(b2), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(b2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#score cell cycle after 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

b2 <- CellCycleScoring(b2, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

#merged_runs$CC.Difference <- merged_runs$S.Score - merged_runs$G2M.Score
all.genes <- rownames(b2)
b2 <- ScaleData(b2, features = all.genes)

b2 <- RunPCA(b2, verbose = F)

b2 <- RunHarmony(b2, "run", assay.use = "RNA")
b2 <- FindNeighbors(b2, reduction = "harmony", dims = 1:30, verbose = FALSE)
b2 <- FindClusters(b2, resolution = 0.8, verbose = FALSE)
b2 <- RunUMAP(b2, reduction = "harmony", assay = "RNA", dims = 1:30, verbose = FALSE, reduction.name = "umapharmony")

DimPlot(b2, reduction = "umapharmony", label = T, pt.size = 0.5, label.size = 4, group.by = c("seurat_clusters", "barcode", "run"), ncol = 3)

FeaturePlot(b2, reduction = "umapharmony", features = c("percent.mt"), cols = c("yellow", "blue"), ncol = 3, label = TRUE)


b2$celltype <-  plyr::mapvalues(
  x = b2$seurat_clusters, 
  from = as.character(c(seq(0, 10, by = 1))), 
  to = as.character(c("ML1", "LP1", "LP2","LP3", "LP4", "Basal1", "Basal2", "ML2", "ML3", "ML4", "ML5")))


b2$broadcelltype <-  plyr::mapvalues(
  x = b2$celltype, 
  from = as.character(c("LP1", "ML1", "LP2", "LP3", "Basal1", "Basal2", "LP4", 
                        "ML2", "ML3", "ML4", "ML5")), 
  to = as.character(c("Luminal progenitor", "Mature luminal", "Luminal progenitor", "Luminal progenitor", "Basal", "Basal", "Luminal progenitor", 
                      "Mature luminal", "Mature luminal", "Mature luminal", "Mature luminal")))

DimPlot(b2, reduction = "umapharmony", label = T, pt.size = 1, label.size = 4, 
        group.by = c("celltype"), ncol = 1,
        order = rev(c("Basal1", "Basal2", "LP1", "LP2", "LP3", "LP4", "ML1", "ML2", "ML3", "ML4", "ML5")))


