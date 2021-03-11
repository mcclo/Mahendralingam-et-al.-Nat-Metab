####
#"Plotting moving average CytoTRACE"
####
#this is example code for just the basal subset generated from "7_CytoTrace_differentiation_state_analysis.R"
#based on code generously provided by G. Gulati from the original CytoTRACE publication.

#Functions
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
MoveAve <- function(x, width) {
  as.vector(stats::filter(x, rep(1/width, width), sides=1));
}

census <- function(mat, counts) {
  
  xnl <- 2^data.matrix(mat)-1
  rs <- colSums(xnl)
  rnorm <- t(t(xnl)*counts/rs)
  A <- x <- log2(1+rnorm)
  
  return(A)
}


#Define important variables
pheno <- basal$celltype
cyto <- basal$cytoTRACE

#pull out data for genes of interest (log-normalized counts)
LBH <- as.matrix(GetAssayData(basal, slot = "data")["LBH", ])
KRT5 <- as.matrix(GetAssayData(basal, slot = "data")["KRT5", ])
ACTA2 <- as.matrix(GetAssayData(basal, slot = "data")["ACTA2", ])
TNFRSF11A <-  as.matrix(GetAssayData(basal, slot = "data")["TNFRSF11A", ])

#calculate RNAcontent per cell
RNAcontent <- range01(colSums(mat_basal))

#Order phenotypes by median CytoTRACE for boxplot
ord1 <- aggregate(cyto, list(pheno), median)
pheno <- factor(pheno)
pheno <- factor(pheno, levels(pheno)[order(ord1$x)])

#Subsample phenotypes for visualization
exm <- table(pheno)
set.seed(123)
for(i in 1:length(exm)){
  if(exm[i] > 50){
    s1 <- sample(which(pheno == names(exm[i])), 50)
    s2 <- pheno[s1]
  } else {
    s2 <- pheno[which(pheno == names(exm[i]))]
  }
  ifelse(i == 1, s3 <- as.character(unlist(s2)), s3 <- c(s3, as.character(unlist(s2))))
}

cyto2 <- cyto[which(pheno %in% s3)]
RNAcontent <- RNAcontent[which(pheno %in% s3)]
pheno2 <- pheno[which(pheno %in% s3)]

#colors matching from dimplot
cols <-  c("#E69F00","#56B4E9")
cols2 <- cols[which(colnames(mat_basal) %in% names(s3))]


##graph
layout(matrix(c(1:3), nrow = 3, byrow = TRUE), height = c(3.25, 3, 2))
par(oma = c(5,0,0,0))
par(mar = c(2.5, 10,0, 10))
par(xpd = T)
boxplot(range01(-cyto)~pheno, outline = F,las = 2,xaxt = "n", yaxt = "n",
        staplelwd = 2,
        medlwd = 2,
        whisklty = 1,
        border = cols,
        col = adjustcolor(cols, alpha.f = 0.25),
        ylab = "",
        xlab = "",
        cex.lab = 1.5,
        las = 1,
        frame.plot = F,
        ylim = c(0,1),
        xlim = c(0, length(levels(pheno))+1),
        cex.axis = 1.5,
        horizontal = TRUE)
title(main = "",font.main = 1, cex.main = 2)
mtext("Cell state",
      cex = 2, side = 2, line = 2.5)
axis(2, pos = 0, at = 1:length(levels(pheno)), labels = levels(pheno), cex.axis = 2.5)
segments(0,length(levels(pheno)) + 1,0,0)
text(-0.025, seq_along(1:length(levels(pheno))), cex =2,
     labels = rep("", length(levels(pheno))), srt = 0, adj = 1, xpd = TRUE)
axis(1, pos = 0, at = seq(0, 1, 0.1),
     labels = rev(format(seq(0,1,0.1), nsmall = 1)), cex.axis = 2)
segments(0,0,1,0)
n = 1
stripchart(range01(-cyto2)~pheno2, vertical = FALSE,
           method = "jitter", add = TRUE, pch = 16, col = cols, bg = 'white',
           lwd = 1, cex = 0.75)
par(mar = c(0, 5, 1, 5))
y0 <- 0
y1 <- 1
step <- 0.1
move <- MoveAve(RNAcontent[order(-cyto)],200)


#for gene expression
par(mar = c(2.5, 10,0, 10))
par(xpd = T)
plot(0,col = "white", type = "l",
     lwd =1, pch = 16, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1),
     cex.axis = 1, cex.lab = 2, frame = F, yaxt = "n", xaxt = "n")

m1 <- -cyto[order(-cyto)][!is.na(move)]
m2 <- MoveAve(LBH[order(-cyto)][!is.na(move)],200)
m1 <- m1[!is.na(m2)]
m2 <- m2[!is.na(m2)]
test <- smooth.spline(range01(m1),range01(m2), spar = 0.85)
x1 <- range01(as.numeric(unlist(predict(test, range01(m1))$x)))
y1 <- range01(as.numeric(unlist(predict(test, range01(m1))$y)))
lines(x1, y1, col = "firebrick", lwd = 4)

m1 <- -cyto[order(-cyto)][!is.na(move)]
m2 <- MoveAve(TNFRSF11A[order(-cyto)][!is.na(move)],200)
m1 <- m1[!is.na(m2)]
m2 <- m2[!is.na(m2)]
test <- smooth.spline(range01(m1),range01(m2), spar = 0.85)
x1 <- range01(as.numeric(unlist(predict(test, range01(m1))$x)))
y1 <- range01(as.numeric(unlist(predict(test, range01(m1))$y)))
lines(x1, y1, col = "orange", lwd = 4)


axis(2, pos = 0, at = seq(0, 1, 0.2), labels = format(seq(0,1,0.2), nsmall = 1),
     las = 2, cex.axis = 2)
segments(0, y0,0,y1)
mtext("Scaled log2\ngene expression",
      cex = 2, side = 2, line = 4)
axis(1, pos =y0, at = seq(0, 1, 0.1),
     labels = rev(format(seq(0,1,0.1), nsmall = 1)), cex.axis = 2)
legend(0.9, 0.8, c("LBH", "TNFRSF11A"), fill = c("firebrick", "orange"), cex = 2, border = NULL, bty = "n")


plot(0,col = "white", type = "l",
     lwd =1, pch = 16, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1),
     cex.axis = 1, cex.lab = 2, frame = F, yaxt = "n", xaxt = "n")

m1 <- -cyto[order(-cyto)][!is.na(move)]
m2 <- MoveAve(KRT5[order(-cyto)][!is.na(move)],200)
m1 <- m1[!is.na(m2)]
m2 <- m2[!is.na(m2)]
test <- smooth.spline(range01(m1),range01(m2), spar = 0.85)
x1 <- range01(as.numeric(unlist(predict(test, range01(m1))$x)))
y1 <- range01(as.numeric(unlist(predict(test, range01(m1))$y)))
lines(x1, y1, col = "dodgerblue", lwd = 4)

m1 <- -cyto[order(-cyto)][!is.na(move)]
m2 <- MoveAve(ACTA2[order(-cyto)][!is.na(move)],200)
m1 <- m1[!is.na(m2)]
m2 <- m2[!is.na(m2)]
test <- smooth.spline(range01(m1),range01(m2), spar = 0.85)
x1 <- range01(as.numeric(unlist(predict(test, range01(m1))$x)))
y1 <- range01(as.numeric(unlist(predict(test, range01(m1))$y)))
lines(x1, y1, col = "purple", lwd = 4)

axis(2, pos = 0, at = seq(0, 1, 0.2), labels = format(seq(0,1,0.2), nsmall = 1),
     las = 2, cex.axis = 2)
segments(0, y0,0,y1)
mtext("Scaled log2\ngene expression",
      cex = 2, side = 2, line = 4)
mtext("CytoTRACE",
      cex = 2, side = 1, line =2.5)
axis(1, pos =y0, at = seq(0, 1, 0.1),
     labels = rev(format(seq(0,1,0.1), nsmall = 1)), cex.axis = 2)
legend(1.00, 0.8, c("KRT5", "ACTA2"), fill = c("dodgerblue", "purple"), cex = 2, border = NULL, bty = "n")





