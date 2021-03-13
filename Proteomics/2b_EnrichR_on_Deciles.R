# This script can be used to recursively run pathway analysis (via enrichR) on genes in each decile for each cell type

# Load libraries
library(openxlsx) # install.packages("openxlsx")
library(enrichR) # install.packages("enrichR")

# This function reads a table from an Excel file into a list, where each column is a list element
read_list_from_excel <- function(file, sheet){
  file = openxlsx::read.xlsx(file, sheet=sheet)
  file_list = lapply(file, function(x)  x[!is.na(x)])
  return(file_list)
}

# Get name of Excel file produced in previous script
genes_per_decile_file <- list.files(pattern="genes_per_decile_in_each_cell.xlsx")

# Read lists of genes, each sheet was a different cell type
openxlsx::getSheetNames(genes_per_decile_file)
BC <- read_list_from_excel(genes_per_decile_file, sheet = "BC Genes in Each Decile")
ML <- read_list_from_excel(genes_per_decile_file, sheet = "ML Genes in Each Decile")
LP <- read_list_from_excel(genes_per_decile_file, sheet = "LP Genes in Each Decile")

# Define function to iteratively run enrichR on each column in excel sheet/list element
add_enrichr_table_sheet <- function(wb, sheet_name, list_of_genes_per_decile){
  # Initialize an empty data frame # TODO optimize with vectorization
  main_df <- data.frame(matrix(nrow = 0, ncol = 16))
  
  for(decile in 1:10){
    # Run enrichR pathway analysis on specific decile gene list
    enriched <- enrichr(list_of_genes_per_decile[[decile]],databases = "GO_Biological_Process_2018")
    # Get GO BP element
    df <- enriched[["GO_Biological_Process_2018"]]
    # Sort by adj p-value, ascending
    df <- df[order(df$Adjusted.P.value),]
    # Get top 5 of ordered table
    df <- df[1:5, ]
    # Add relevant columns to beginning of table
    df <- cbind(Decile = decile, 
                Pathway = gsub(" \\(.*", "", df$Term), 
                Adjusted.P.Value = df$Adjusted.P.value, 
                "",
                df)
    # Add to main data frame
    main_df <- rbind(main_df, df)
  }
  # Rename column names
  colnames(main_df) <- colnames(df)
  
  # Add as Excel sheet
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, x = main_df)
  
  return(wb)
}

# Create a new Excel workbook
wb <- createWorkbook()
# Add Enrichr info for each cell type into respective sheets
add_enrichr_table_sheet(wb, sheet_name = "BC EnrichR Output", BC)
add_enrichr_table_sheet(wb, sheet_name = "ML EnrichR Output", ML)
add_enrichr_table_sheet(wb, sheet_name = "LP EnrichR Output", LP)
# Save as file
saveWorkbook(wb, sprintf("%s_enrichr_per_decile_in_each_cell2.xlsx", format(Sys.Date(), format="%Y%m%d")), overwrite = T)
