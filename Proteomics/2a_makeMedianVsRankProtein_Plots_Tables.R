# This script can be used to make protein distribution plots that show the total number of proteins detected, stratified into deciles
# The plots show the log2 median intensity of a protein vs its rank (rank 1 = highest intensity)
# Make a plot for each of the cell types - total proteome of B, ML, LP, stroma

#' Input:
#' -total_ibaq values (non-imputed, raw expression, averaged across cell types)

#' Output:
#' -distribution plots (.png)
#' -distribution plots (.svg) - no labels or axes, just points and axis ticks
#' -Excel file with "decile table", including protein rank and number in decile
#' -Excel file with tabs for gene symbols of proteins in each decile

#' Notes: 
#' - The script is set up to create plots without any axis labels/lines/legend, but can be undone by removing theme() arguments in ggplot
#' - Refer to Ankit Sinha's Cancer Cell 2019 paper figure 1B
#' - Stroma was included but is exclusive from other cell types for this type of analysis

# Load required libraries
library(openxlsx) # for saving worksheets in Excel file
library(ggplot2) # plotting function
library(StatMeasures) # for decile() 
library(RColorBrewer) # for making color pallette
library(svglite) # for saving as svg

# Define plotting function --------------------------------------------
# This function makes the plot, that shows the log2 median intensity of a protein vs its rank.
# It takes in non-imputed expression values (i,e. with zeros) as dataframe x, plot (and filename) title.
make_distribution_plot <- function(avg_x, combat_x, non_imputed_x, plot_title) {
  # Keep rows that are non-0s ie. take out rows that are 0
  non_zero_indices <- which(avg_x != 0)
  combat_x <- combat_x[non_zero_indices,]
  non_imputed_x <- non_imputed_x[non_zero_indices,] 
 
  # Find the number of patients protein is found by calculating sum 
  num_patient_observed <- apply(non_imputed_x, 1, function(x) sum(x!=0))
  
  # Replace original 0s in non-imputed matrix with NA in combat matrix
  combat_x[non_imputed_x == 0] <- NA
  
  # log2 transform the matrix and apply median computation to all rows 
  log2_med_int <- apply(combat_x, 1, function(x) median(x, na.rm = TRUE))
  
  # Make data frame with num patients observation and log2 median density, as well as decile
  df <- data.frame(Decile = factor(decile(log2_med_int, decreasing = TRUE)),
                   Gene = names(log2_med_int),
                   log2_Med_Intensity = log2_med_int,
                   Rank = rank(-log2_med_int),
                   No.Patients_Observed_In = num_patient_observed#factor(num_patient_observed, levels = rev(levels(factor(num_patient_observed))))
  )
  
  # Sort based on Rank (should order accordingly into deciles 1 to 10 as well)
  df <- df[order(df$Rank),]
  
  # Have two files, one with labels(png), one for importing into adobe illustrator (svg)
  filename <- sprintf("%s_prot_distribution_plot.png", plot_title)
  filename2 <- sprintf("%s_prot_distribution_plot.svg", plot_title)

  plot_title <- sprintf("%s distribution (%d proteins)", plot_title, length(non_zero_indices))
  p <- ggplot(df, aes(x=Rank, y=log2_Med_Intensity, label=Gene, color = Decile))+
    geom_point(shape=21, size=3)+
    # geom_label_repel(aes(label = GeneName))+
    scale_color_manual(values=RColorBrewer::brewer.pal(n = 10, name = "Spectral"))+
    # geom_text(data=subset(df, #logical#), aes(10_Obs, log2_Med_Intensity,label=name), size=1.5) +
    xlab("Protein Rank") +
    ylab("Log2 Median Intensity") +
    xlim(0, 6150)+
    ylim(0, 35)+
    theme(panel.grid.major = element_blank(), #remove grid in plot
          panel.grid.minor = element_blank(),
          legend.key=element_blank(),
          panel.background = element_blank(), #clear background
          axis.line = element_line(colour = "black"))+
    ggtitle(plot_title)

  # save map
  ggsave(filename, p)

  # Plot "clean" map (ie no labels/axis titles/text)
  p2 <- ggplot(df, aes(x=Rank, y=log2_Med_Intensity,label=Gene, color = Decile))+
    geom_point(shape=21, size=3)+
    # geom_label_repel(aes(label = GeneName))+
    scale_color_manual(values=RColorBrewer::brewer.pal(n = 10, name = "Spectral"))+
    # geom_text(data=subset(df, #logical#), aes(10_Obs, log2_Med_Intensity,label=name), size=1.5) +
    xlab("Protein Rank") +
    ylab("Log2 Median Intensity") +
    xlim(0, 6150)+
    ylim(0, 35)+
    # scale_x_discrete(limits = c(1:10)) +
    theme(panel.grid.major = element_blank(), #remove grid in plot
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), #clear background
          # axis.line = element_line(colour = "black"),
          axis.line = element_blank(),
          legend.key=element_blank(),
          axis.title = element_blank(),
          title = element_blank(),
          axis.text = element_blank(),
          legend.position = "none")+# coord_fixed(ratio = 1.001) +
    ggtitle(plot_title)
  # Save plot
  ggsave(filename2, p2)
  
  return(df)
}

# Helper function
# Define function to write a list (named) with vectors of varying sizes to Excel
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

# Run analysis ------------------------------------------------------
# Load data 
load("h.total_pint.nonimputed.RData")
load("h.total_proteome_mathepan.RData")
h.pint.nonimputed.avg <- h.pint.nonimputed.avg[order(row.names(h.pint.nonimputed.avg)),]
h.pint <- h.pint[order(h.pint$Gene.names), ]

# Check row order is the same
all(rownames(h.pint.combat_matrix) == rownames(h.pint.nonimputed.avg))
all(h.pint$Gene.names == rownames(h.pint.nonimputed.avg))
# Check column order is the same
all(colnames(h.pint.combat_matrix) == colnames(h.pint[,1:29]))


# Defne the 3 main groups #indices/columns that match the cell type to proteome dataset
cell_columns <- list(BC = grep("BC", colnames(h.pint)), 
                     ML = grep("LC", colnames(h.pint)), 
                     LP = grep("LP", colnames(h.pint)))

plot_titles <- c("BC_population", "ML_population", "LP_population") #names of titles/file names


# Initialize list to hold data frames (ie. the decile tables)
decile_tables <- list(BC=NULL, ML=NULL, LP=NULL)

# Create Excel workbook
wb <-  createWorkbook(sprintf("%s_decileTable", format(Sys.Date(), "%Y%m%d")))

# Run the plotting function on each of the cell types 
for(i in 1:3){ #TODO try sapply
  # Make table that shows deciles in each cell type
  columns = cell_columns[[i]]
  df <- make_distribution_plot(h.pint.nonimputed.avg[,i], h.pint.combat_matrix[,columns], h.pint[,columns],plot_titles[i]) 
  # Save data frame to list for further analysis
  decile_tables[[i]] <- df
  
  # Write data to new sheet in Excel
  # a) define name of sheet
  current_sheet <- paste(plot_titles[i], "Decile Table", sep=" ")
  # b) add worksheet to workbook
  addWorksheet(wb, sheetName = current_sheet)
  # c) write data to Excel file
  writeData(wb, sheet = current_sheet, x = df)
}

# Save Excel workbook to file
saveWorkbook(wb, file = sprintf("%s_TotalProteomeRankedDecileTable.xlsx", format(Sys.Date(), "%Y%m%d")), overwrite = TRUE)

# Get genes per decile for pathway enrichment
# Create a workbook to write data to Excel
wb <- createWorkbook("GenesInEachDecile")

#In a loop, create a 
for (i in 1:3){
  # Define current table (cell type)
  current_table <- decile_tables[[i]]
  
  # Unfactor gene names (messes with gProfiler output)
  current_table$Gene <- as.character(current_table$Gene)
  
  # Initialize list
  list_decile_genes <-vector("list",10)
  
  # Run apply to create list of genes in each decile
  list_decile_genes <- sapply(1:10, function(j){
    current_indices <- which(current_table$Decile == j)
    gene_list <- as.character(current_table$Gene[current_indices])
    gene_list <- vapply(strsplit(gene_list,";"), `[`, 1, FUN.VALUE=character(1))
    return(gene_list)
  })
  # Rename list elements as Decile Number
  names(list_decile_genes) <- paste("Decile", 1:10, sep="_")
  
  # Add worksheet to workbook 
  # Define worksheet name
  current_sheet <- paste(names(decile_tables[i]), "Genes in Each Decile", sep = " ")
  # Write data to worksheet
  addWorksheet(wb, current_sheet)
  col_num <- 1
  # In a loop, make a column for the signatures' ensembl ids
  for(i in 1:10){
    # rename columns to specify cell type
    writeData(wb, sheet = current_sheet, x = names(list_decile_genes[i]), startCol = col_num, startRow = 1)
    # write data to worksheet
    writeData(wb, sheet = current_sheet, x = list_decile_genes[[i]], startCol = col_num, startRow = 2)
    # update counter
    col_num <- col_num + 1
  }
}
# Save Excel workbook to file
filename <- paste(format(Sys.Date(), format="%Y%m%d"), "genes_per_decile_in_each_cell.xlsx", sep="_")
saveWorkbook(wb, file = filename, overwrite = TRUE)  