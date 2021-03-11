
library(enrichR)
library(viridis)
library(genesorteR)
library(tidyr)
library(readxl)
library(tidyverse)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

#here we are loading the metabolism genelist to subset the main seurat object to just these genes of interest before
#performing differential expression testing (ie/ finding metabolic markers of each cell type/cluster)
metab_list <- read_xlsx("F:/MULTI-seq/possemato_metabolic_gene_list_nature10350-s2.xlsx", sheet = 2, col_names = F)
metab_list <- as.matrix(metab_list)
metab_list <- as.character(metab_list)

#subset based on the metabolism gene list
b2_met <- subset(b2, features = metab_list)

#to see how this looks, can make a dendrogram
Idents(b2_met) <- "celltype"
b2_met <- BuildClusterTree(b2_met)

# pull the tree
data.tree <- Tool(object = b2_met, slot = "BuildClusterTree")
# plot the tree
#extract order of tree in order to get colors. 
data.tree$tip.label

mycol <-  c("Basal1"="#ED1C24", "Basal2"="#ED1C24", "LP1" = "#00AEEF", "LP2" = "#00AEEF", "LP3" = "#00AEEF", "LP4"="#00AEEF", "ML1"="#2E3192", "ML2"="#2E3192", "ML3"="#2E3192", "ML4"="#2E3192",
            "ML5"="#2E3192")

ape::plot.phylo(x = data.tree, direction = "downwards",  show.node.label = FALSE,
               tip.color = mycol, show.tip.label = T, node.pos = 1, srt = 90, adj = 0.5,
                font = 2, edge.width = 6, cex = 4) 



#Using the GenesortR package and ClusterProfiler R package to get top differentially expressed genes per MEC cluster
#subset on each cell type and make corrplot?

Idents(b2_met) <- "celltype"

# median method enriches for genes with a higher detection rate than adaptiveMedian, but the latter may be cleaner - dont use naive
gs = sortGenes(b2_met@assays$RNA@data, Idents(b2_met), binarizeMethod = "median")

head(gs$specScore)

mm = getMarkers(gs, quant = 0.99)

pp = plotMarkerHeat(gs$inputMat, gs$inputClass, mm$markers, clusterGenes=TRUE, outs = TRUE)

plotTopMarkerHeat(gs, top_n = 5, outs = TRUE, plotheat = T)

pp = getPValues(gs) 

names(which(apply(pp$adjpval, 1, function(x) any(x < 0.05))))

rownames(pp$adjpval[which(pp$adjpval[,1] < 0.05),])

tab = getTable(gs, pp, islog = T, adjpval_cutoff = 0.025)

#diagnostic plot, logFC should correlate with spec score
plot(tab$Average.Log.Fold.Change, tab$Specificity.Score, col = as.factor(tab$Cluster), pch = 20)

unique(tab$Gene.Name)

#test clusters
getClassAUC(gs, plotCurves = T)


#plot corr of top 100 gene per cluster
corr = genesorteR::plotCorrelationHeat(gs, markers = unlist(plotTopMarkerHeat(gs, top_n = 100, plotheat = FALSE, outs = TRUE)), displayNumbers=FALSE, outs = TRUE)

#load GO database
msig = msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol)

#category = "C5", subcategory = "GO:BP"
#msig = data.frame(msig$gs_name[which(msig$gs_subcat == "CP:KEGG")], msig$gene_symbol[which(msig$gs_subcat == "CP:KEGG")])

topMarkerGenes = genesorteR::plotTopMarkerHeat(gs, top_n = 50, averageCells = 10^6, outs = TRUE, plotheat = F)

#To pull genelists for each cluster (for example, I took these lists for input into ssGSEA analysis for the metabric cohort)
LP1_enrich <- topMarkerGenes$LP1
LP2_enrich <- topMarkerGenes$LP2
LP3_enrich <- topMarkerGenes$LP3
LP4_enrich <- topMarkerGenes$LP4

ML1_enrich <- topMarkerGenes$ML1
ML2_enrich <- topMarkerGenes$ML2
ML3_enrich <- topMarkerGenes$ML3
ML4_enrich <- topMarkerGenes$ML4
ML5_enrich <- topMarkerGenes$ML5

basal1_enrich <- topMarkerGenes$Basal1
basal2_enrich <- topMarkerGenes$Basal2

#set options for the compareCluster call below
ont_option = "BP"
pvalueCutoff = 0.025
qvalueCutoff = 0.025
maxGSSize = 300
minGSSize = 3

#using enrichGO with C5 geneset
resGO <- compareCluster(geneClusters = topMarkerGenes, pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff,
                        qvalueCutoff = qvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize,
                        fun = "enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", readable = F, ont = ont_option)

res_simplifyGO <- clusterProfiler::simplify(x = resGO, cutoff = 0.65, by = "p.adjust")

clusterProfiler::dotplot(res_simplifyGO, showCategory = 5, font.size = 15, color = "p.adjust") + 
  theme(axis.text.y = element_text(size = 18, color = "black", family = "sans"),
        axis.text.x = element_text(size = 17, color = "black", family = "sans"),
        legend.text = element_text(size = 16, color = "black", family = "sans"),
        legend.title = element_text(size = 16, color = "black", family = "sans")) 


#plot custom dotplot

#pull out pvalue, geneRatio, to plot 
c <- res_simplifyGO@compareClusterResult
frac <- c$GeneRatio
frac <- sapply(frac, function(x) eval(parse(text=x)))
c$GeneRatio <- frac

#need to pull out top 5 terms 
  c_final <- c %>% group_by(Cluster) %>%  top_n(Count, n = 100) %>% slice(1:6) %>%
    mutate(Description = fct_reorder(Description, p.adjust))
  
  ggplot(c_final, mapping = aes(x= Cluster, y= Description, stroke = 2, color = Cluster, size = GeneRatio)) +
  geom_point()  +
  theme_bw()  +
  scale_shape_discrete(solid = T) +
  scale_color_manual(values = c("Basal1"="#ED1C24", "Basal2"="#ED1C24", "LP1" = "#00AEEF", 
                                "LP2" = "#00AEEF", "LP3" = "#00AEEF", "LP4"="#00AEEF",
                                "ML1"="#2E3192", "ML2"="#2E3192", "ML3"="#2E3192",
                                "ML4"="#2E3192", "ML5"="#2E3192")) + 
  theme(axis.text = element_text(size = 16, color = "black", family = "sans"),
        axis.title = element_text(size = 20, color = "black", family = "sans"),
        legend.text = element_text(size = 16, color = "black", family = "sans"),
        legend.title =element_text(size = 20, color = "black", family = "sans")) +
  guides(color = guide_legend(override.aes = list(size = 8))) +
  guides(color = guide_legend(override.aes = list(size = 8))) +
  theme(strip.text.x = element_text(size = 22, colour = "black", family = "sans"))


#to make heatmap
  
mycol <-  c("Basal1"="#ED1C24", "Basal2"="#ED1C24", "LP1" = "#00AEEF", "LP2" = "#00AEEF", "LP3" = "#00AEEF", "LP4"="#00AEEF", "ML1"="#2E3192", "ML2"="#2E3192", "ML3"="#2E3192", "ML4"="#2E3192",
              "ML5"="#2E3192")

b2_met$celltype <- factor(b2_met$celltype, levels = c("Basal1", "Basal2","LP1", "LP2", "LP3", "LP4", "ML1", "ML2", "ML3", "ML4", "ML5"))

#from genesorter
topMarkerGenes_heat = genesorteR::plotTopMarkerHeat(gs, top_n = 5, averageCells = 10^6, outs = TRUE, plotheat = F)
genes <- stack(topMarkerGenes_heat)

DoHeatmap(subset(b2_met, downsample = 50), features = genes$values, 
          group.by = "celltype",
          group.colors = mycol,
          size = 8, hjust = 0.5,
          raster = F,
          angle = 0, draw.lines = TRUE, lines.width = 4) + 
  NoLegend() +
  theme(axis.text.y = element_text(size = 16, color = "black", family = "sans")) +
  scale_fill_viridis(option = "viridis", na.value = "white")

