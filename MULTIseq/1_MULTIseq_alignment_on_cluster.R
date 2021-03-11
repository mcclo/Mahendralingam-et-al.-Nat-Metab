#load libraries
library(KernSmooth)
library(reshape2)
library(Rtsne)
library(stringdist)
library(ShortRead)
library(deMULTIplex)
library(ggplot2)
library(dplyr)
#-------------------------------------RUN NUMBER 1 -----------------------------------------------------
## Define vectors for reference barcode sequences and cell IDs - these are four example barcodes for the four patients multiplexed
bar.ref1 <- c("GCACACGC", "AGAGAGAG", "CGAGATTC", "GTAGCACT")

#you need to move the barcode.tsv.gz and matched .fastq file into your working directory on the cluster before running this code
#*****a laptop cannot handle this step*******
cell.id.vec1 <- read.table("/cluster/home/barcodes1.tsv.gz")
cell.id.vec1 <- gsub("-1", "", cell.id.vec1[,1])
readTable1 <- MULTIseq.preProcess(R1 = '/cluster/home/Barcode_S14_L003_R1_001.fastq.gz', R2 = '/cluster/home/Barcode_S14_L003_R2_001.fastq.gz', cellIDs = cell.id.vec1, cell=c(1,16), umi=c(17,28), tag=c(1,8))

bar.table1 <- MULTIseq.align(readTable1, cell.id.vec1, bar.ref1)

rm(readTable1)

save.image("/cluster/home/MULTIseq_RUN1.RData")


