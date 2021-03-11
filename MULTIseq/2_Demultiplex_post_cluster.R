#load libraries
library(KernSmooth)
library(reshape2)
library(Rtsne)
library(stringdist)
library(ShortRead)
library(deMULTIplex)
library(ggplot2)
library(dplyr)
library(Seurat)

#------------------------------------
#This code was run on the veryhimem cluster using the shell script "MULTIseq_RUN1.sh"

# #-------------------------------------RUN NUMBER 1 -----------------------------------------------------
# ## Define vectors for reference barcode sequences and cell IDs - these are four example barcodes for the four patients multiplexed
# bar.ref1 <- c("GCACACGC", "AGAGAGAG", "CGAGATTC", "GTAGCACT")
# 
# #you need to move the barcode.tsv.gz and matched .fastq file into your working directory on the cluster before running this code
# #*****a laptop cannot handle this step*******
# cell.id.vec1 <- read.table("/cluster/home/barcodes1.tsv.gz")
# cell.id.vec1 <- gsub("-1", "", cell.id.vec1[,1])
# readTable1 <- MULTIseq.preProcess(R1 = '/cluster/home/Barcode_S14_L003_R1_001.fastq.gz', R2 = '/cluster/home/Barcode_S14_L003_R2_001.fastq.gz', cellIDs = cell.id.vec1, cell=c(1,16), umi=c(17,28), tag=c(1,8))
# 
# bar.table1 <- MULTIseq.align(readTable1, cell.id.vec1, bar.ref1)
# 
# rm(readTable1)
# 
# save.image("/cluster/home/MULTIseq_RUN1.RData")

#-----------------------------------

#Now import the MULTIseq_RUN1.RData file containing bar.table1
#Now rename columns from Bar#-> patient/sample ID
colnames(bar.table1)[1] <- "sample1"
colnames(bar.table1)[2] <- "sample2"
colnames(bar.table1)[3] <- "sample3"
colnames(bar.table1)[4] <- "sample4"

## Visualize barcode space calling only the columns with barcodes
bar.tsne <- barTSNE(bar.table1[,1:4]) 
## Note: Exclude columns 97:98 (assuming 96 barcodes were used) which provide total barcode UMI counts for each cell (ignore this if using small number of barcodes)

pdf("MULTI_run1_barcodes.pdf")
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "right") 
  print(g)
}
dev.off()

#this .pdf will appear in your working directory.


## Round 1 -----------------------------------------------------------------------------------------------------
## Perform Quantile Sweep
#Shave off the nUMI columns in bar.table1
bar.table1.full <- bar.table1
bar.table1 <- bar.table1[,1:4]

good.bars <- c("sample1", "sample2", "sample3", "sample4")
bar.table1 <- bar.table1[, good.bars]  # Remove missing bars and summary columns. 

#This step assigns cell ID's to quantiles based on the lognormalized barcode UMI counts.
bar.table1_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table1_sweep.list[[n]] <- classifyCells(bar.table1, q=q)
  names(bar.table1_sweep.list)[n] <- paste("q=",q,sep="")
}

## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 <- findThresh(call.list=bar.table1_sweep.list)
#the plot computed below shows the proportion of cells that fall within quantiles based on their log normalized barcode UMI counts. 
ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "right") +
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))

## Finalize round 1 classifications, remove negative cells
round1.calls <- classifyCells(bar.table1, q=findQ(threshold.results1$res, threshold.results1$extrema))
neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]
#now remove negative cells from bar.table1 (note, backup info stored in bar.table1.full)
bar.table1 <- bar.table1[-which(rownames(bar.table1) %in% neg.cells), ]


## the online tutorial says to repeat until all no negative cells remain (usually 3 rounds)
final.calls <- round1.calls
#to view distrubtion of demultiplex calls
table(final.calls)

#-------I've skipped this step since I've had few numbers of negative cells
## Perform semi-supervised negative cell reclassification
reclass.cells <- findReclassCells(bar.table1, names(final.calls)[which(final.calls=="Negative")])
reclass.res <- rescueCells(bar.table1, final.calls, reclass.cells)

## Finalize negative cell rescue results
final.calls.rescued <- final.calls
rescue.ind <- which(reclass.cells$ClassStability >= 15) ## Note: Value will be dataset-specific.
final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]
#-------

#take final.calls.rescued and group with tSNE coordinate to see what ended up getting called doublet, negative, or singlet
bar.tsne$Classification <- "Singlet"
bar.tsne$Classification[which(final.calls.rescued[rownames(bar.tsne)]=="Doublet")] <- "Doublet"
bar.tsne$Classification[which(final.calls.rescued[rownames(bar.tsne)] == "Negative")] <- "Negative"

#to determine how many cells are in each classification
table(bar.tsne$Classification)

#to view overall distrubution of calls
tsne_classification <-ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2)) + 
  geom_point(size = 0.25, aes(color = Classification)) + 
  theme_void()
ggsave(tsne_classification, file = "workingdirectory/RUN1_classification_tsne.png", 
       dpi = 600, width = 8, height = 5.5)

win.graph()
tsne_classification


#take unique samples that you have (number of barcodes), then plot tsne showing which cells have which barcode. 
samples <- unique(final.calls.rescued)
#use the below code to alphabetize or arrange all unique labels within final.calls
sorted_samples <- unique(final.calls.rescued)
sorted_samples <- unlist(sorted_samples)
sorted_samples <- sort(sorted_samples)
sorted_samples <- as.list(sorted_samples)


#plot tSNE barcode space
plotSampleTSNE <- function(sorted_samples){
  data <- bar.tsne
  data$Sample <- "Other"
  data$Sample[which(final.calls.rescued[rownames(bar.tsne)] == sorted_samples)] <- sorted_samples
  sample_plot <- ggplot(data, aes(x = TSNE1, y = TSNE2)) +
    geom_point(size = 0.25, alpha = 0.5, aes(color = Sample)) +
    scale_color_manual(values=c("red","lightgrey")) +
    theme_void() }
 

#the below call adds all plots together ie/ all unique call groups (bar1, bar2, bar3, singlet, doublet, negative)
barcode_tsne_plots <- lapply(sorted_samples, plotSampleTSNE)
barcode_tsne_plots_combined <- CombinePlots(plots = barcode_tsne_plots)
win.graph()
barcode_tsne_plots_combined

#Now save the full analysis as .RData file. These final calls will be added as a new metadata column in your Seurat object
final.calls.RUN1 <- final.calls.rescued
save(final.calls.RUN1, file = "final_calls_RUN1.RData")