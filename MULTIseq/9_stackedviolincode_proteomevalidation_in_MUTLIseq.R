

library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(ggsignif)

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab("Expression level") + ggtitle(feature) + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 90), 
          axis.text.y = element_text(size = rel(1), angle = 0), 
          plot.margin = plot.margin ) +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "crossbar", width = 0.5, show.legend = F) +
    geom_signif(comparisons = my_comparisons, test = "wilcox.test",
                map_signif_level = FALSE, textsize=4, y_position = c(3.5,4.5,5.5)) +
    ylim(NA, 7.0)
    
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 5)
  return(p)
}


fun_mean <- function(x){
  return(round(data.frame(y=mean(x),label=mean(x,na.rm=T)), digits = 4)) }


mycol <-  c("Luminal progenitor" = "#00AEEF", "Mature luminal"="#2E3192", "Basal"="#ED1C24")
my_comparisons <- list( c("Luminal progenitor", "Mature luminal"), c("Luminal progenitor", "Basal"), c("Mature luminal", "Basal"))

Idents(b2) <- "broadcelltype"



#combined for paper 
b2$broadcelltype <- factor(b2$broadcelltype, levels = c("Basal", "Luminal progenitor", "Mature luminal"))

features<- c("PKM", "GAPDH", "IDH3A",  "FBP1", "ACLY", "ACO1", "IDH2", "MDH2",
             "COQ5", "COQ9", "SDHC", "SUCLA2",
             "NDUFA2", "NDUFA10", "NDUFA12", "NDUFB5","NDUFS1", "NDUFS2", "NDUFS3", "NDUFS5", "NDUFS7", "NDUFS8",
  "NDUFV1", "NDUFV3", "MT-CO2", "MT-ND4", "MT-ND5",
  "UQCRC1", "UQCRC2", "COX5B", "ATP5F1A", "ATP5PD", "ATP5PO")
StackedVlnPlot(obj = b2, features = features, col = mycol) 



