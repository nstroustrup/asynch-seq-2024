plotExampleGeneDVA <- function(norm, 
                               annots,
                               color,
                               genesymbol, 
                               geneid, 
                               legend_colors = NULL,
                               legend_linetype = NULL,
                               legend_colors_name = NULL, 
                               legend_linetype_name = NULL, 
                               bins = 25,
                               linetype = NULL,
                               type = c("hist", "boxplot", "ecdf")) {
  
  type <- match.arg(type)
  
  genename <- geneid[genesymbol, on = "GeneSymbol"]$GeneName[1]
  
  ggdf <- cbind(Norm = norm[genename, ], annots)
  
  gp <- ggplot(ggdf) +
    scale_x_log10() +
    theme(legend.position = "top") +
    ggtitle(genesymbol)
  
  if (type == "hist") {
    gp <- gp + 
      geom_histogram(aes_string("Norm + 1", fill = color), bins = bins, position = "identity", alpha = .7) +
      xlab("Normalized Read Counts") +
      ylab("Number of Individual Nematodes")
    
  } else if (type == "boxplot") {
    gp <- gp + 
      geom_jitter(width = .4) + 
      geom_boxplot(aes_string(color, "Norm + 1", color = color), outlier.colour = NA) +
      ylab("Normalized Read Counts")
    
  } else if (type == "ecdf") {
    gp <- gp + 
      stat_ecdf(aes_string("Norm + 1", color = color, linetype = linetype)) +
      xlab("Normalized Read Counts") +
      ylab("ECDF")
  }
  
  if (! is.null(legend_colors)) {
    gp <- gp + 
      scale_fill_manual(values = legend_colors,  name = legend_colors_name) + 
      scale_color_manual(values = legend_colors,  name = legend_colors_name)
  }
  
  if (! is.null(legend_linetype)) {
    gp <- gp + 
      scale_linetype_manual(values = legend_linetype,  name = legend_linetype_name)
  }
  
  gp 
}