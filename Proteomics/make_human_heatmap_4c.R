#' This function, plot h (human) heatmap2, takes in the following to create a heat map with an annotation bar for significant hits.
#' 
#' Arguments:
#' exp_matrix = numeric matrix of expression data with columns as samples and genes as rows
#' hc_samples, hc_proteins = hc_clust objects indicating hierarchal clustering 
#' num_replicates = number of replicates of each cell type (synonymous with number of patients)
#' genes_of_interest = vector of genes used to annotate rows
#' genes_of_interest_label = string indicating group name of these genes of interest
#' filename =  string - name of the file with directory and format as .pdf, .jpeg, etc
#' title = string - title of the heat map 
#' celltypes = vector of the cell type names eg. c("B", "LM", "LP", "S")

#' heat map colors:
#' purple - low 
#' green - high 

plot_h_heatmap2 <-function(exp_matrix, hc_samples, hc_proteins, filename=NA, title=NA, pint.pheno=NA, cluster_vector=NA){
  # Call libaries (usually at beginning of script)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  
  # Rescale expression matrix
  exp_matrix <- rescale(exp_matrix, to = c(-2, 2))

  # celltypes <- c("BC", "LM", "LP")
  # num_replicates <- 10  #9 for basal cell
  phases <- c("Follicular", "Luteal", "Post.Menopausal")
  
  # Make annotations -----------------    
  # Make annotation column (sample information)
  annotation_col <- data.frame(
    CellType = gsub("LC", "LM", pint.pheno$cell_type),
    # PatientAge = rep(c("31-40","31-40","31-40", "<30","<30","31-40", "51-60", "51-60",  ">60", ">60"), length(celltypes)),
    PatientAge = pint.pheno$patient_age_group,
    HormoneStatus = pint.pheno$hormone_status
  )
  rownames(annotation_col) <- colnames(exp_matrix)
  
  # Make annotation row (gene information)
  annotation_row = data.frame(
    SignificantHits = cluster_vector
  )
  rownames(annotation_row) <- rownames(exp_matrix)
  
  # Define annotation colors -----------------
  # For age variable
  newAgeCols <- colorRampPalette(c("white", "black")) #"yellowgreen
  myAgeColors <- newAgeCols(length(unique(annotation_col$PatientAge))+1)
  names(myAgeColors) <-c("<30", "31-40", "41-50", "51-60",  ">60")
  # For hormone status variable
  newPhaseCols <- colorRampPalette(c("lightgoldenrodyellow", "goldenrod2"))
  myPhaseColors <- newPhaseCols(length(phases))
  names(myPhaseColors) <- phases
  
  annotation_colors <- list(
    CellType = c(BC="red", LP="skyblue", LM="blue"),#, S="darkgray"),
    HormoneStatus = myPhaseColors,
    PatientAge = myAgeColors,
    SignificantHits =  c(BC="red", LM="blue", LP="skyblue", Unassigned="white")
    #Pathway = c(genes[genes%in%hdr_genes]="blue", genes[!genes%in%hdr_genes]="black")
  )
  
  # pheatmap function call -----------------    
  pheatmap(exp_matrix, 
           color = colorRampPalette(c("darkgreen", "forestgreen", "white", "darkorchid3", "darkorchid4"))(50),
           scale = "none", #"row",
           annotation_col = annotation_col,
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           cluster_cols = hc_samples,
           cluster_rows = hc_proteins,
           border_color = NA,
           # cutree_rows = 10,
           cellwidth = 15,
           show_rownames = FALSE,
           show_colnames = FALSE,
           treeheight_col = 10,
           treeheight_row = 50,
           # cellheight = 1,
           # fontsize_row = 1,
           # height = 10,
           main=title,
           filename = filename
  )
}
