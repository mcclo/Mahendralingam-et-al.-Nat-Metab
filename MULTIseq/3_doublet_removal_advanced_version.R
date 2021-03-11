library(Seurat)
library(ggplot2)
library(sctransform)
library(affy)
library(dplyr)
library(stringr)
library(pheatmap)
library(umap)
library(DoubletFinder)

#setting global options is key for the PrepSCTIntegration step, the first number (20000) indicates im setting 
#it to 20Gb of space allotted for the calculation. 
options(future.globals.maxSize = 20000 * 1024^2)
#RUN1 

RUN1_import <- Read10X(data.dir = "MULTI-seq_directory/_Multiseq_3pr_v3/outs/filtered_feature_bc_matrix")
RUN1_SCT <- CreateSeuratObject(counts = RUN1_import, project = "Breast", names.field = 1, names.delim = "-", min.cells = 3)

RUN1_SCT$run <- "RUN1"

#this part is for MULTIseq
#import final_call_RUN1.RData file
#need to make final.calls into a character instead of a list
final.calls.RUN1 <- unlist(final.calls.RUN1)
#make a new metadata column in your seurat object for the final calls.
RUN1_SCT$barcode <- final.calls.RUN1

#Add new column with BRCA status that maps to the cells MULTIseq barcode 
RUN1_SCT$Genotype <- plyr::mapvalues(
  x = RUN1_SCT$barcode, 
  from = c("sample1", "sample2", "sample3", "sample4"), 
  to = c("RM", "RM", "mutneg", "BRCA1")) 

#first only look at cells with appropriate barcodes
#subset out just the negative cells to enable comparison between MULTIseq doublets and DoubletFinder doublets
RUN1_SCT <- subset(RUN1_SCT, subset = barcode %in% c("sample1", "sample2", "sample3", "sample4", "Doublet"))

#Now preprocess EACH Seurat dataset separately up to the SCTtransform step. 
#next, go through standard quality control metrics, firstly, add column with %mitochondrial RNA
RUN1_SCT[["percent.mt"]] <- PercentageFeatureSet(RUN1_SCT, pattern = "^MT-")

win.graph()
VlnPlot(RUN1_SCT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(RUN1_SCT, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(RUN1_SCT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#subset, vars.to.regress function for percent.mt
RUN1_SCT <- subset(RUN1_SCT, subset = nFeature_RNA > 200 & percent.mt < 15)
#-----------

###Note: Need to apply DoubletFinder to each individual run prior to aggregating for batch correction, only use SCTransform in this step, later analysis does not use SCTransform.
RUN1_SCT <- SCTransform(RUN1_SCT, verbose = FALSE, return.only.var.genes = FALSE)

RUN1_SCT <- RunPCA(RUN1_SCT)
RUN1_SCT <- RunUMAP(RUN1_SCT, dims = 1:30)

RUN1_SCT <- FindNeighbors(RUN1_SCT, dims = 1:30, verbose = FALSE)
RUN1_SCT <- FindClusters(RUN1_SCT, resolution = 0.8, verbose = FALSE)

win.graph()
DimPlot(RUN1_SCT, group.by= c("barcode", "seurat_clusters"), reduction="umap", pt.size=0.5)


##pK Identfication (no ground-truth)
sweep.res.list.RUN1 <- paramSweep_v3(RUN1_SCT, PCs = 1:30, sct = TRUE)
sweep.stats_RUN1 <- summarizeSweep(sweep.res.list.RUN1, GT = FALSE)
bcmvn_RUN1 <- find.pK(sweep.stats_RUN1)


#Homotypic doublet proportion estimation
homotypic.prop <- modelHomotypic(RUN1_SCT@meta.data$seurat_clusters)   ## ex: annotations <- RUN1_SCT@meta.data$ClusteringResults
nExp_poi <- round(0.134*length(RUN1_SCT@active.ident))  ## Assuming 13.4% doublet formation rate - this number came from Demultiplex doublet rate calculated from table(bar.tsne$Classification) doublet/total cells*100
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
RUN1_SCT <- doubletFinder_v3(RUN1_SCT, PCs = 1:30, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
#Run with adjusted nExp (the reuse.pANN column name is found in metadata)
RUN1_SCT <- doubletFinder_v3(RUN1_SCT, PCs = 1:30, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_1757", sct = TRUE)


DimPlot(RUN1_SCT, group.by= c("DF.classifications_0.25_0.005_1631", "barcode"), reduction="umap", pt.size=0.5)
win.graph()

DimPlot(RUN1_SCT, group.by= c("barcode", "seurat_clusters"), reduction="umap", pt.size=0.5)

#to visualize and compare the table of homotypic vs heterotypic doublets Identification (MULTIseq)
RUN1_SCT$singlets <- plyr::mapvalues(
  x = RUN1_SCT$barcode, 
  from = c("sample1", "sample2", "sample3", "sample4", "Doublet"), 
  to = c("singlet_homo", "singlet_homo", "singlet_homo", "singlet_homo", "Doublet_homo")) 

table(RUN1_SCT$DF.classifications_0.25_0.005_1631, RUN1_SCT$barcode)
table(RUN1_SCT$barcode)
#Now we have doublet classification from both MULTIseq demultiplex and DoubletFinder 
#Need to pull out cell ID's in sequence to avoid doublets identified from both methods, save this as a vector
#then pull out these cell IDs in the original Seurat Object prior to SCT transform in order to recluster minus the doublets. 

#Subset singlets from MULTIseq demultiplex
RUN1_singlets <- subset(RUN1_SCT, subset = barcode %in% c("sample1", "sample2", "sample3", "sample4"))
#Remove additional doublets found by DoubletFinder using the DF.classification
RUN1_singlets <- subset(RUN1_singlets, subset = DF.classifications_0.25_0.005_1631 %in% c("Singlet"))

table(RUN1_singlets$DF.classifications_0.25_0.005_1631, RUN1_singlets$barcode)

Idents(object=RUN1_singlets) <- "orig.ident"
RUN1_singlet_list <- WhichCells(RUN1_singlets, idents = 'Breast')

#-----------------Save workspace and then run code below and save final seurat object separately. 

#Now pull out CellIDs in order to finalize the RUN1_SCT file for downstream merging or integration (batch correction)
#Create clean Seurat object
RUN1_SCT_final <- CreateSeuratObject(counts = RUN1_import, project = "Breast", names.field = 1, names.delim = "-", min.cells = 3)

RUN1_SCT_final$run <- "RUN1"

final.calls.RUN1 <- unlist(final.calls.RUN1)
RUN1_SCT_final$barcode <- final.calls.RUN1

#Add new column with BRCA status that maps to the cells MULTIseq barcode as an example
RUN1_SCT_final$Genotype <- plyr::mapvalues(
  x = RUN1_SCT_final$barcode, 
  from = c("sample1", "sample2", "sample3", "sample4"), 
  to = c("RM", "RM", "BRCA1", "BRCA1")) 

RUN1_SCT_final[["percent.mt"]] <- PercentageFeatureSet(RUN1_SCT_final, pattern = "^MT-")


#This subsetting will account for mito percent and gene number correction (<15% mito, >200 nfeatures)
RUN1_SCT_final <- subset(RUN1_SCT_final, cells = RUN1_singlet_list)

save(RUN1_SCT_final, file = "RUN1_SCT_final.RData")
