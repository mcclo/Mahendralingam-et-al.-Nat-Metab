#' This function checks whether an R package is installed and loads it into the system
#' list.of.packages is a vector of package names 
#' 
#' ** does not work with Bioconductor packages unless already installed #TODO 
load_packages <- function(list.of.packages){
  # Find packages that are not currently installed
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  # Install these packages
  if(length(new.packages)) 
    install.packages(new.packages)
  # Run require() on the list of packages, which is equivalent to library() but gives a return value
  lapply(list.of.packages, require, character.only = TRUE)
}

# Filtering and cleaning data
CleanRnaSeqFunc <- function(rnaSeq) {
  library(dplyr)
  # Remove duplicated gene symbols
  rnaSeq <- metabricExpression[!duplicated(metabricExpression$Symbol),]
  rownames(rnaSeq) <- rnaSeq$Symbol
  # Select columns
  rnaSeq <- rnaSeq %>% select(
    -CLID,
    -Symbol
  )
  # Remove duplicated columns
  rnaSeq <- rnaSeq[!duplicated(colnames(rnaSeq)),]
  # return
  as.matrix(as.data.frame(rnaSeq))
}

# Read signatures to a list
read_list_from_excel <- function(file){
  library(readxl)
  lapply(read_xlsx(file), function(x) x[!is.na(x)])
}

# Stats functions ====================================
# Performs t-test for multiple comparisons (ie. all subtypes vs each other)
ttestFunc <- function (ssGSEA.subtype, signature) {
  # Subset data columns to only signature of interest
  signature.values <- ssGSEA.subtype[,signature]
  subtype <- ssGSEA.subtype$subtype
  # Do pairwise t test
  test <- pairwise.t.test(x = signature.values, g = subtype, p.adjust.method = "fdr")
  # Print values to console
  print(signature)
  print(test$method)
  print(test$p.value)
}

# ANOVA Func
anovaTestFunc <- function(ssGSEA.subtype, signature){
  # run 1 way ANOVA across all subtypes
  res.aov <- aov(formula = ssGSEA.subtype[,signature] ~ ssGSEA.subtype$subtype)
  return(res.aov)
}

# Plotting functions =================================
# Violin Plot Func
# dfAux - data frame with "subtype" and "signature" columns, where signatures contain the ssGSEA values for the signature
drawViolinPlot <- function(dfAux){
  # Create violin plot
  g <- ggplot(dfAux, aes(x = subtype, y = signature)) + 
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill=NA)+#, outlier.shape = NA) +
    # stat_summary(fun.y=median, geom="line", size=2, color="red")+
    geom_jitter(width = 0.2)+
    # scale_y_continuous(name = "ssGSEA Enrichment Score", breaks = seq(-2, 2, .2))+
    ylim(-0.42, 0.69)+
    theme_classic()+
    # scale_fill_manual(values = unname(c(color = brewer.pal(length(levels(dfAux$subtype)), "Greys"))))++
    labs(xlab = "Breast Cancer Subtype",
         ylab = "ssGSEA Enrichment Score",
         title = signature)
    
  return(g)
}

# Box Plot Func
drawBoxPlot <- function(dfAux){
   # Make boxplot 
  g <- ggplot(dfAux, aes(x = subtype, y = signature))+ 
    geom_boxplot(data = dfAux) +
    # geom_jitter(width=0.1)+
    # scale_y_continuous(name = "ssGSEA Enrichment Score", breaks = seq(-2, 2, .2))+
    ylim(-0.42, 0.69)+
    # scale_x_discrete(labels=c("claudin-low"="CL", "Basal"="BC", "LumA"="LA", "LumB"="LB", "Her2"="H2", "NC"="NC", "Normal"="No"))+
    theme_classic() +
    labs(xlab = "Breast Cancer Subtype",
         ylab = "ssGSEA Enrichment Score",
         title = signature)
    
  return(g)
}
