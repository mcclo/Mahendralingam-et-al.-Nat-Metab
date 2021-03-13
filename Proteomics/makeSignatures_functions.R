# These functions are limited to the analysis of this data.
# Check out devtools::install_github("kazeera/kazDGE") for more generic functions

# This function averages and returns the columns of expression matrix according to celltype
average_hmpint_columns <- function(x){
  cell_columns <- list(BC = grep("BC", colnames(x)), 
                       ML = grep("LC", colnames(x)), 
                       LP = grep("LP", colnames(x)))
  
  df <-data.frame(BC=apply(x[,cell_columns$BC], 1, mean),
                  ML=apply(x[,cell_columns$ML], 1, mean),
                  LP=apply(x[,cell_columns$LP], 1, mean))
  rownames(df) <- rownames(x)
  return(df)    
}

# Main function to do ANOVA and Tukey HSD----------------------------
# This function takes in a protein expression table (where rownames are gene ids) and computes the ANOVA and Tukey test p-vals for each gene
compute_ANOVA_tukey_pvals <- function(exp_matrix){
  #indices/columns that match the cell type to proteome dataset
  # cell_indices <- list(BC = 1:9, ML = 10:19, LP = 20:29) 
  cell_indices <- list(BC = grep("BC", colnames(exp_matrix)), 
                       ML = grep("LC", colnames(exp_matrix)), 
                       LP = grep("LP", colnames(exp_matrix)))
  # run ANOVA and Tukey on each row (gene) provided in the data frame using apply
  #returns a list of dataframes
  ANOVA_tukey <- apply(exp_matrix, 1, function(x){
    # Make expression values into a 3 x 10 table with three columns corresponding to each cell type
    df <- data.frame(
      Expression= unlist(c(x[cell_indices$BC], x[cell_indices$ML], x[cell_indices$LP])),
      CellType =factor(rep(c("BC", "ML", "LP"), times=c(length(x[cell_indices$BC]), length(x[cell_indices$ML]), length(x[cell_indices$LP]))))
    )
    # run 1 way ANOVA across all values
    res.aov <- aov(formula = Expression ~ CellType, data=df)
    #reformat output of anova
    anova_output <- anova(res.aov)
    # extract p-value from output
    anova_pval <- anova_output$`Pr(>F)`[1]
    # run TukeyHSD test
    tukey_output <- TukeyHSD(res.aov)
    #extract p-values of all 3 comparisons (in column 4)
    tukey_pvals <- tukey_output$CellType[,4]
    # create table with p-value results
    significance_df <- data.frame(ANOVA_pval = anova_pval,
                                  Tukey_pval_LP.BC = unname(tukey_pvals[grep("LP-BC", names(tukey_pvals))]),
                                  Tukey_pval_ML.BC = unname(tukey_pvals[grep("ML-BC", names(tukey_pvals))]),
                                  Tukey_pval_ML.LP = unname(tukey_pvals[grep("ML-LP", names(tukey_pvals))])
    )
    return(significance_df)
  })
  # reformat list of dataframes to a single dataframe
  significance_df <- do.call(rbind, unname(Map(cbind, Gene = names(ANOVA_tukey), ANOVA_tukey)))
  #colnames "Gene, ANOVA_pval, Tukey_pval_LP.BC, Tukey_pval_ML.BC, Tukey_pval_ML.LP"
  return(significance_df)
}

# Main function to compute log2 FC (Fold Change)-----------------------
# Takes in log2 expression data and returns a matrix comparing the 3 cell types to each other
compute_logFC_3groups <- function(exp_matrix_avg){
  # Make expression values into a 3 x 10 table with three columns corresponding to each cell type
  # Data is already log2 transformed, so we can subtract them
  # log2FC = log2(x/y) = log2(x) - log2(y)
  log2FC_df2 <- data.frame(log2FC_BC.vs.ML = exp_matrix_avg$BC - exp_matrix_avg$ML, 
                           log2FC_BC.vs.LP = exp_matrix_avg$BC - exp_matrix_avg$LP, 
                           log2FC_LP.vs.BC = exp_matrix_avg$LP - exp_matrix_avg$BC, 
                           log2FC_LP.vs.ML = exp_matrix_avg$LP - exp_matrix_avg$ML, 
                           log2FC_ML.vs.BC = exp_matrix_avg$ML - exp_matrix_avg$BC, 
                           log2FC_ML.vs.LP = exp_matrix_avg$ML - exp_matrix_avg$LP
                           
  )
  return(log2FC_df2)
}

# Helper functions-------
subset_df <- function(significance_df, log2FC_df, ID) {
  log2FC_ID_columns <- grep(sprintf("log2FC_%s", ID), colnames(log2FC_df))
  Tukey_ID_columns <- which(grepl(sprintf("(?=.*%s)(?=.*Tukey_pval)", ID), colnames(significance_df), perl = TRUE))
  
  subset_df <- data.frame(Gene = significance_df$Gene,
                          ANOVA_pval = significance_df$ANOVA_pval,
                          significance_df[,Tukey_ID_columns],
                          log2FC_df[,log2FC_ID_columns])
  return(subset_df)
}

# This function takes in a data frame with Tukey p-values and log2FC values for a cell type 
# It returns the original table with t
add_sig_columns <- function(subset_df, pval_thres, log2FC_thres_up){
  Tukey_columns <- grep("Tukey_pval", colnames(subset_df))
  log2FC_columns <- grep("log2FC", colnames(subset_df))
  
  Tukey_Significant <- ifelse((subset_df[,Tukey_columns[1]] < pval_thres & subset_df[,Tukey_columns[2]] < pval_thres), "Yes", "")
  log2FC_Cutoff_Met  <- ifelse((subset_df[,log2FC_columns[1]] > log2FC_thres_up & subset_df[,log2FC_columns[2]] > log2FC_thres_up), "Yes", "")
  
  both_conditions_met <- ifelse((Tukey_Significant == "Yes" & log2FC_Cutoff_Met == "Yes"), "Yes", "")
  
  subset_df <- cbind(subset_df, 
                     Tukey_Significance_Met = Tukey_Significant,
                     log2FC_Cutoff_Met = log2FC_Cutoff_Met,
                     Both_Conditions_Met = both_conditions_met)
  return(subset_df)
}

# This function writes a list (named) with vectors of varying sizes to Excel
# It will make a worksheet in the workbook (wb) specified with the name of "current_sheet"
write_list_to_Excel <- function(current_sheet, my_list, wb){
  # # call required package to write to Excel
  require(openxlsx)
    # add a worksheet to the workbook with the name specified in the argument 
  addWorksheet(wb, current_sheet)
  col_num <- 1
  # In a loop, make a column for the signatures' ensembl ids
  for(i in 1:length(my_list)){
    # rename columns to specify cell type
    writeData(wb, sheet = current_sheet, x = names(my_list[i]), startCol = col_num, startRow = 1)
    # write data to worksheet
    writeData(wb, sheet = current_sheet, x = my_list[[i]], startCol = col_num, startRow = 2)
    # update counter
    col_num <- col_num + 1
  }
  return(wb)
}
