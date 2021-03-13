#' Helper functions -- 
#' 
#' This function returns a list of genes, where the position of the gene name matches the position in the original gene_list (in which gene names are concatenated by ";")
#' if there is a gene name, it is present in the signatures, if not then it is absent
#' 
#' e.g. INPUTS gene_list = "a;b;c", "c;d;f", "d", "l;d;a"  
#'             signatures = "l", "q", "b     
#'      OUTPUT "b", "", "", "l" 
get_present_gene_name <- function(gene_column, signatures){
  # Split the gene symbols delimited by ";" into seperate entities as elements of a list
  gene_column <- strsplit(gene_column,";")
  gene_column <- toupper(gene_column)
  signatures <- toupper(signatures)
  # Remove the ".number" that was appended to make the gene names unique
  # gene_column <- lapply(gene_column, function(x) gsub("\\..*","",x))
  # Look for any genes present in the signature in each split list vector 
  present_genes <- lapply(gene_column, function(x){
    if(any(x %in% signatures)){
      genes <- x[x %in% signatures]  #find whether genes are in signatures
      genes <- paste(genes, collapse=";") #concatenate using ";" if necessary
    } 
    else return("")
  })
  # unlist(sapply(1:6222, function(x){ if(length(present_genes[[x]]) > 1) return(paste(present_genes[[x]], collapse=";")) })
  # Unlist list to vector form
  present_genes <- unlist(present_genes)
  
  # Return a vector of genes that are present in the signatures ("" if absent), conserving the order of the original gene list
  return(present_genes)
}

# This function finds a subset of table (data frame with expression values) by matching its gene_column (vector/column in table) to the signatures (vector)
subset_proteome_table2 <- function(table, gene_column, signatures){
  # Find genes that match the genes of interest list. Non-matches become ""
  present_genes <- get_present_gene_name(gene_column, signatures)
  # Subset
  table <- table[present_genes != "",]
  rownames(table) <- make.unique(present_genes[present_genes != ""], sep = ".")
  # Returns subset of table 
  return (table)
}

#' This function plots a heatmap --
#' Arguments: 
#' full proteome - all protein expression data, matrix where row names correspond to proteins/gene symbols and columns correspond to sample info
#' signatures - list of genes to subset
#' title - name of heat map and filename id
#' celltypes - vector with (B, ML, LP, S cell types) according to order of full proteome
make_human_heatmap <- function(full_proteome, signatures, title=NA, filename=NA){
  # Call libraries (usually at beginning of script)
  library(pheatmap)
  library(RColorBrewer)
  library(cluster)
  library(scales)
  
  # Define variables for heat map ----
  # Make annotation column (sample information)
  celltypes <- c("BC", "ML", "LP")
  num_replicates <- 10  #9 for basal cell
  hormoneStatuses <- c("Follicular", "Luteal", "Post.Menopausal")
  
  # Subset proteome ----
  raw_value_mat <- subset_proteome_table2(full_proteome, rownames(full_proteome), signatures)
  
  # Compute z-scores to plot on heat map
  z_score_mat <- apply(raw_value_mat, 1, function(x) (x - mean(x)) / sd(x))
  exp_matrix <- t(z_score_mat)
  # # Rescale
  # exp_matrix <- rescale(exp_matrix, to = c(-2, 2))
  
  # Make hierarchical cluster objects (dendogram - like structures)----
  # Clustering method ----
  cluster_method <- "pearson" 
  
  # a) Cluster samples (genes)
  hc_samples <- as.dist(1-cor(exp_matrix, method =cluster_method))
  #use = "pairwise.complete.obs",
  hc_samples <- as.dendrogram(diana(hc_samples)) #DIvisive ANAlysis Clustering

  # b) Cluster rows (genes)
  hc_proteins <- as.dist(1-cor(t(exp_matrix), method =cluster_method))
                                      #use = "pairwise.complete.obs",
  hc_proteins <- as.dendrogram(diana(hc_proteins)) #DIvisive ANAlysis Clustering
  
  # Make annotations -----------------    
  #make annotation column (sample information)
  annotation_col <- data.frame(
    CellType = gsub("LC", "ML", pint.pheno$cell_type)
    # PatientAge = pint.pheno$patient_age_group,
    # HormoneStatus = pint.pheno$hormone_status
  )
  annotation_col <- data.frame(CellType = annotation_col)
  rownames(annotation_col) <- colnames(exp_matrix)
  
  # Make annotation row (gene information) -- not used right now
  annotation_row = NA
  # rownames(annotation_row) <- rownames(exp_matrix)
  
  # Define annotation colors 
  # a) Pick colors for the age groups
  newAgeCols <- colorRampPalette(c("white", "black")) #"yellowgreen
  myAgeColors <- newAgeCols(5)#patient samples are stratified into 5 patient age groups
  names(myAgeColors) <-c("<30", "31-40", "41-50", "51-60",  ">60")

  # b) Pick colors for the age groups
  newHormoneStatusCols <- colorRampPalette(c("lightgoldenrodyellow", "goldenrod2"))
  myHormoneStatusColors <- newHormoneStatusCols(length(hormoneStatuses))
  names(myHormoneStatusColors) <- hormoneStatuses

  # c) Pick colors for the cell-types
  myCellColors = c(BC="red", LP="skyblue", ML="blue")#, S="darkgray")
  
  # d) Make a list of the annotation colors
  annotation_colors <- list(
    CellType = myCellColors,
    HormoneStatus = myHormoneStatusColors,
    PatientAge = myAgeColors
  )
  
  # Make heat map ----------
  # a) Define title and filename if desired
  title_for_plot <- sprintf("%s (%s proteins)", title, nrow(exp_matrix))
 
  # b) pheatmap function call
  pheatmap(exp_matrix, 
           color = colorRampPalette(c("darkolivegreen", "darkgreen", "white", "darkorchid3", "darkorchid4"))(50),
           scale = "row",#"none", #already scaled
           annotation_col = annotation_col,
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           cluster_cols =  as.hclust(hc_samples),
           cluster_rows = as.hclust(hc_proteins),
           border_color = NA,
           cellheight = 6,
           # cutree_rows = 11,
           cellwidth = 10,
           # show_rownames = FALSE,
           show_colnames = FALSE,
           treeheight_col = 10,
           treeheight_row = 10,
           fontsize_row = 6,
           height = 6,
           main=NA, #title_for_plot
           filename = filename
  )
}
