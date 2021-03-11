
library(readxl)
library(AUCell)
library(GSEABase)
library(mixtools)
library(Seurat)
library(SingleCellExperiment)
library(escape)
library(dittoSeq)
library(ggridges)
library(limma)
library(Matrix)
library(cowplot)
library(Nebulosa)
library(RColorBrewer)
#from escape package can use
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))

#load 50 Hallmark genesets within escape package 
GS <- getGeneSets(library = "H")

#load in signatures/custom genesets 
B.genes <- read.delim("F:/MULTI-seq/Current/MAT_B.txt", header = FALSE, sep = "\t" , stringsAsFactors = FALSE)
LP.genes <- read.delim("F:/MULTI-seq/Current/MAT_LP.txt", header = FALSE, sep = "\t" , stringsAsFactors = FALSE)
ML.genes <- read.delim("F:/MULTI-seq/Current/MAT_ML.txt", header = FALSE, sep = "\t" , stringsAsFactors = FALSE)

#convert to characters
B.genes <- as.list(B.genes$V1)
B.genes <- as.character(B.genes)
LP.genes <- as.list(LP.genes$V1)
LP.genes <- as.character(LP.genes)
ML.genes <- as.list(ML.genes$V1)
ML.genes <- as.character(ML.genes)

#create GeneSet object from the B and ML lists 
B_g <- GeneSet(B.genes, setIdentifier = "B", setName="B_genes")
LP_g <- GeneSet(LP.genes, setIdentifier = "LP", setName="LP_genes")
ML <- GeneSet(ML.genes, setIdentifier = "ML", setName="ML_genes")

geneSets <- GeneSetCollection(B_g, LP_g, ML)


#object b2
#enrich for the signatures on the raw counts for each group of genesets (could combine into a single call if you want)
ES <- enrichIt(obj = b2, gene.sets = geneSets, groups = 1000)
ES_base <- enrichIt(obj = b2, gene.sets = GS, groups = 1000) 

#add enrichment scores to seurat object metadata 
b2 <- AddMetaData(b2, ES)
b2 <- AddMetaData(b2, ES_base)

#visualize complex patterns using the ggridges R package 
Idents(b2) <- "celltype"
minusML <- subset(b2, idents =  c("LP1", "LP2", "LP3", "LP4","Basal1", "Basal2"))

#here is an example subsetting seurat object to just plot certain clusters on ridgeplot below
ES2 <- data.frame(minusML[[]], Idents(minusML))

#need to factor the cell type column in order to get a specified order in the ridgeplot
ES2$celltype <- factor(ES2$celltype,levels = rev(c("LP1", "LP2", "LP3", "LP4", "Basal1", "Basal2")))


#here are three ridgeplot for proteomes derived metabolic signatures from Mats paper
B_rid <- ridgeEnrichment(ES2, gene.set = "B_genes", group = "celltype", add.rug = TRUE,
                colors = c("LP1" = "#00AEEF", "LP2" = "#00AEEF", "LP3" = "#00AEEF", "LP4"="#00AEEF", "ML1"="#2E3192", "ML2"="#2E3192", "ML3"="#2E3192", "ML4"="#2E3192",
                           "ML5"="#2E3192", "Basal1"="#ED1C24", "Basal2"="#ED1C24")) +
  theme(axis.text.x= element_blank(),
        plot.title = element_text(lineheight=.8, face="bold"),
        axis.text = element_text(size = 35, color = "black", family = "sans"),
        axis.title.x=element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.title.y = element_blank(),
        legend.title=element_text(size=10,face = 'bold'),
        legend.text=element_text(size=10),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = 'white'),
        plot.background = element_blank(),
        legend.key = element_rect(color = "transparent", fill = "transparent")) +
  geom_density_ridges(quantile_lines = T, quantile_fun = function(x,...)mean(x), size = 1.25)
  

LP_rid <- ridgeEnrichment(ES2, gene.set = "LP_genes", group = "celltype", add.rug = TRUE,
                          colors = c("LP1" = "#00AEEF", "LP2" = "#00AEEF", "LP3" = "#00AEEF", "LP4"="#00AEEF", "ML1"="#2E3192", "ML2"="#2E3192", "ML3"="#2E3192", "ML4"="#2E3192",
                                     "ML5"="#2E3192", "Basal1"="#ED1C24", "Basal2"="#ED1C24")) +
  theme(axis.text.x= element_blank(),
        plot.title = element_text(lineheight=.8, face="bold"),
        axis.text = element_text(size = 35, color = "black", family = "sans"),
        axis.title.x=element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.title.y = element_blank(),
        legend.title=element_text(size=10,face = 'bold'),
        legend.text=element_text(size=10),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = 'white'),
        plot.background = element_blank(),
        legend.key = element_rect(color = "transparent", fill = "transparent")) +
  geom_density_ridges(quantile_lines = T, quantile_fun = function(x,...)mean(x), size = 1.25)

ML_rid <- ridgeEnrichment(ES2, gene.set = "ML_genes", group = "celltype", add.rug = TRUE,
                          colors = c("LP1" = "#00AEEF", "LP2" = "#00AEEF", "LP3" = "#00AEEF", "LP4"="#00AEEF", "ML1"="#2E3192", "ML2"="#2E3192", "ML3"="#2E3192", "ML4"="#2E3192",
                                     "ML5"="#2E3192", "Basal1"="#ED1C24", "Basal2"="#ED1C24")) +
  theme(axis.text.x= element_blank(),
        plot.title = element_text(lineheight=.8, face="bold"),
        axis.text = element_text(size = 35, color = "black", family = "sans"),
        axis.title.x=element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.title.y = element_blank(),
        legend.title=element_text(size=10,face = 'bold'),
        legend.text=element_text(size=10),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = 'white'),
        plot.background = element_blank(),
        legend.key = element_rect(color = "transparent", fill = "transparent")) +
  geom_density_ridges(quantile_lines = T, quantile_fun = function(x,...)mean(x), size = 1.25)

#combine all three plots into a nice even grid using cowplot package
plot_grid(B_rid, LP_rid, ML_rid, ncol = 3)

#for complementary UMAP plots of the ridgeplots:
Bsig <- FeaturePlot(b2, "B_genes", cols = c("yellow", "lightblue", "blue"), pt.size = 0.8) +
  theme(axis.text = element_text(size = 25, color = "black", family = "sans"),
        axis.title=element_text(size = 25, color = "black", family = "sans"),
        legend.text=element_text(size=20)) +
        ggtitle(NULL)
LPsig <- FeaturePlot(b2, "LP_genes", cols = c("yellow", "lightblue", "blue"), pt.size = 0.8)+
  theme(axis.text = element_text(size = 25, color = "black", family = "sans"),
        axis.title.x = element_text(size = 25, color = "black", family = "sans"),
        axis.title.y=element_blank(),
        legend.text=element_text(size=20)) +
  ggtitle(NULL)
MLsig <- FeaturePlot(b2, "ML_genes", cols = c("yellow", "lightblue", "blue"), pt.size = 0.8)+
  theme(axis.text = element_text(size = 25, color = "black", family = "sans"),
        axis.title.x = element_text(size = 25, color = "black", family = "sans"),
        axis.title.y=element_blank(),
        legend.text=element_text(size=20)) +
  ggtitle(NULL)

plot_grid(Bsig, LPsig, MLsig, ncol =3)


#get significance using the Escape package
output <- getSignificance(ES2, group = "celltype", fit = "ANOVA")


#One-way ANOVA with Tukeys post-test - this could is built from the getSignificance function above. 
library(dplyr)
  
  group2 <- ES2[, 'celltype']
  gr_names <- unique(group2)
  input <- select_if(ES2, is.numeric)
  output <- NULL
  
  out <- lapply(input, function(x) aov(x ~ group2))
  
  for (i in seq_along(out)) {
    df <- out[[i]]
    fval <- summary(df)[[1]]$"F value"[[1]]
    pval <- summary(df)[[1]]$"Pr(>F)"[[1]]
    output <- rbind(output, c(fval, pval))
  }
  
  output <- as.data.frame(output)
  colnames(output) <- c("f.value", "p.value")
  rownames(output) <- colnames(input)
  output$FDR <- p.adjust(output$p.value)
  
  #extract geneset aov of interest for Tukey
  
  
  #KAZEERA - if you could figure out how to change the below to perform the anova on each column of "out" 
  o2 <- out$HALLMARK_OXIDATIVE_PHOSPHORYLATION
  
  q <- TukeyHSD(o2)
  
  q$group2

write.csv(q$group2, file = "F:/MULTI-seq/OXPHOS_ridge_anova.csv")
