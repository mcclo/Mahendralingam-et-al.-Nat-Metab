#' This script was used to find enrichment of signatures in PAM50+Claudin-low breast cancer (BC) subtypes via RNA expression
#' 
#' It combines 3 separate parts - 1) prepare metabric data, 2) compute ssGSEA, 3) make box plots (x = BC subtype, y = ssGSEA score)
#' RNA expression data is from the METABRIC study from RData file in the folder. 
#' 
#' This was originally Luis' and Miguel's code but I refactored and optimized it.

# ============================================
# Call libraries
# BiocManager::install("GSVA")
source("metabric_functions.R")
load_packages(c("dplyr", "GSVA", "readxl")) # function in metabric_functions.R

# 1. Prepare input for ssGSEA
# This script filters metabric breast cancer data to only subtypes of interest
# RNA seq data
load(file="metabricExpression.RData")
metabricExpression <- metabricExpression %>% select(Symbol, everything())

# Reformat data table 
rnaSeq.clean <- CleanRnaSeqFunc(metabricExpression) #function in functions.R

# Clinical annotations
load("metabricClinical.RData")

# Keep patient ids that are present in both 
keep_ids <- 
  colnames(rnaSeq.clean) %in% rownames(metabricClinical) %>%
  colnames(rnaSeq.clean)[.]

# Subset clinical annotations
clinAnn <- metabricClinical[keep_ids, ]

# Make a new column with only breast cancer subtypes of interest
x <- clinAnn$Pam50...Claudin.low.subtype %>% as.character 
x[! x %in% c("Basal", "claudin-low", "LumA", "LumB", "Her2")] <- NA
clinAnn$BC.Subtype <- x

# Save variables
save(rnaSeq.clean, clinAnn, file="clean_metabric.RData")

# ============================================
# 2. Create ssGSEAs (ie. Find Enrichment of Signatures in BC expression profiles)
# This script runs an ssGSEA that correlates signatures (gene sets) with Breast Cancer (BC) subtypes (Pam-50 + Claudin Low). 
source("metabric_functions.R")
load_packages(c("dplyr", "GSVA"))#function in functions.R
load("clean_metabric.RData")

# Get signatures # Excel or RDS file 
signature_file <- "..." # insert Excel file name from 4b here

# Read signatures to a list
if(grepl("xlsx", signature_file)){
  signatures <- read_list_from_excel(signature_file)  #function in functions.R
} else { # Rdata
  signatures <- readRDS(signature_file)
}

# Compute ssGSEA matrix
ssGSEA <- 
  gsva(expr = rnaSeq.clean, gset.idx.list = signatures, method="ssgsea", kcdf="Gaussian") %>% 
  t %>%
  data.frame
colnames(ssGSEA) <- paste(colnames(ssGSEA), "ssGSEA", sep="_")

# Subset data
ssGSEA_clean <- ssGSEA[rownames(clinAnn),]
all(clinAnn$Patient.ID == rownames(ssGSEA_clean)) # TRUE

save(ssGSEA_clean, clinAnn, file = "ssGSEA.RData")

# ============================================
# 3. Plot ssGSEA data 
analysis_name <- "Mahendralingam_metabolic_proteomic"

# Read in required data/libraries/functions
source("metabric_functions.R")
load("ssGSEA.RData")
load_packages(c("ggplot2", "ggpubr"))

# Create output folder
output_folder <- sprintf("Plots")
dir.create(output_folder)

# variable of interest - column in annotation table
varAnn <- "BC.Subtype"
keep_rows <- !is.na(clinAnn[,varAnn])

# Colors to annotate BC subtypes
PAM50_colors <- c(LumA="mediumpurple1", LumB="mediumpurple4", Basal="seagreen1", Her2="goldenrod1", "claudin-low"="gray")

# For stat testing - specify all comparisons
elements <- names(PAM50_colors)
comparisons <- gtools::combinations(n=length(elements),r=2,v=elements, repeats.allowed=F)
comparisons <- split(comparisons, seq(nrow(comparisons)))

# Font and line size for plot
font_size <- 11
line_size <- 1

# Create pdf file of all plots
pdf_filename <- sprintf("%s/METABRIC_%s_ssGSEA_indiv.pdf", output_folder, analysis_name)
pdf(pdf_filename, onefile=TRUE)

for(curr in colnames(ssGSEA_clean)){
  # Create data frame
  df <- data.frame(var = clinAnn[keep_rows,varAnn] %>% factor,
                   ssGSEA = ssGSEA_clean[keep_rows,curr])
  
  # Sort in descending order according to subtype
  x <- split(df$ssGSEA, df$var) %>% lapply(median) %>% unlist
  y <- sort(x, decreasing = TRUE, na.last = TRUE)
  df$var <- factor(df$var, levels = names(y), ordered = TRUE)
  
  # Make box plot
  p <- ggplot(df, aes(x = var, y = ssGSEA, fill=var))+
    geom_boxplot(aes(fill = var), outlier.shape = NA)+
    # geom_violin(aes(fill = var), draw_quantiles = c(0.25, 0.5, 0.75))+#, outlier.shape = NA) +
    # geom_point()+
    geom_jitter(size = 0.5, width = 0.1)+
    labs(title = sprintf("%s enrichment in PAM50+CL", curr),
         subtitle = sprintf("METABRIC, %s, t.test", analysis_name),
         x = varAnn, 
         y = "ssGSEA")+
    # theme_classic()+
    theme(panel.background = element_blank(),#remove background color and lines
          axis.line = element_line(colour = 'black', size = line_size), # increase the axis-line thickness and change the color to blac
          # Ticks
          axis.ticks = element_line(colour = "black", size = line_size),#increase the tick thickness)
          # axis.ticks.x = element_line(margin = margin(t = 20, r = 0, b = 0, l = 0)), #increase space between x axis title and labels
          # axis.ticks.y = element_line(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.ticks.length = unit(.25, "cm"),
          # Axes labels
          axis.text = element_text(colour = "black", size = font_size),
          axis.text.x = element_text(margin = margin(t = 7, r = 0, b = 0, l = 0),vjust = 1, hjust = 1, angle=45), #increase space between x axis title and labels
          axis.text.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 0)),
          # axes tick labels
          axis.title = element_text(colour = "black", size = font_size, face = "bold"), # axes title labels
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), #increase space between x axis title and labels
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          # axis.text.x = element_text(angle = 45, hjust = 1),
          # legend
          legend.text = element_text(colour = "black", size = font_size),
          legend.title = element_blank())+ 
    stat_compare_means(method="t.test", comparisons = comparisons)+
    scale_fill_manual(values = PAM50_colors)
  print(p)
}

dev.off()
