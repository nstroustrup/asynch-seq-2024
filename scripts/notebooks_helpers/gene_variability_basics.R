formatDA <- function(basics_da, geneid) {
  mean_df <- as.data.table(basics_da@Results$Mean@Table)
  disp_df <- as.data.table(basics_da@Results$Disp@Table)
  resdisp_df <- as.data.table(basics_da@Results$ResDisp@Table)
  resdisp_df[, MeanOverall := NULL]
  
  da <- Reduce(
    function(x, y) merge(x, y, by = "GeneName", all = TRUE),
    list(mean_df, disp_df, resdisp_df))
  da <- merge(geneid[, .(GeneName, GeneSymbol)], da, all.x = FALSE, all.y = TRUE, by = "GeneName")
}

compareConditionsBaSiCS <- function(data, condition1, condition2, geneid, colors, parameter = c("mean", "disp", "resdisp"), ngenes = 4L) {
  
  
  data <- copy(data)
  
  parameter <- match.arg(parameter, c("mean", "disp", "resdisp"))
  
  
  coefficient_colname <- switch(parameter, mean = "MeanLog2FC", disp = "DispLog2FC", resdisp = "ResDispDistance")
  diff_result_colname <- switch(parameter, mean = "ResultDiffMean", disp = "ResultDiffDisp", resdisp = "ResultDiffResDisp")
  
  # each column correspond to which condition?
  value_colname <- switch(parameter, mean = "Mean", disp = "Disp", resdisp = "ResDisp")
  if (data[get(diff_result_colname) == condition1][[coefficient_colname]][1] > 0) {
    value1_colname <- paste0(value_colname, 1)
    value2_colname <- paste0(value_colname, 2)
  } else {
    value1_colname <- paste0(value_colname, 2)
    value2_colname <- paste0(value_colname, 1)
    data.table::set(data, j = coefficient_colname, value = - data[[coefficient_colname]])
  }
  
  # find genes we will label
  condition1_labeled_genes <- data[get(diff_result_colname) == condition1] %>%
    .[order(abs(get(coefficient_colname)), decreasing = TRUE)] %>%
    .[seq_len(ngenes)]
  
  condition2_labeled_genes <- data[get(diff_result_colname) == condition2] %>%
    .[order(abs(get(coefficient_colname)), decreasing = TRUE)] %>%
    .[seq_len(ngenes)]
  
  # nudges for ggrepel
  nudge <- c(x = 0.3, y = -0.3)
  
  condition1_nudge <- nudge
  condition2_nudge <- nudge
  
  if (condition1_labeled_genes[[coefficient_colname]][1] < 0) {
    condition1_nudge <- - condition1_nudge
  }
  
  if (condition2_labeled_genes[[coefficient_colname]][1] < 0) {
    condition2_nudge <- - condition2_nudge
  }
  
  # reformat colors
  colors <- c(colors, "#9d9683")
  names(colors) <- c(condition1, condition2, "NoDiff")
  
  # plot
  gp <- ggplot(data, aes_string(value1_colname, value2_colname, color = diff_result_colname)) +
    
    geom_point(alpha = .7, size = .8) +
    
    geom_text_repel(aes(label = GeneSymbol),
                    box.padding = 1, 
                    direction="both",
                    nudge_x = condition1_nudge["x"],
                    nudge_y = condition1_nudge["y"],
                    
                    condition1_labeled_genes,
                    alpha = 0.8,
                    color = colors[1],
                    min.segment.length = unit(0, 'lines'),
                    size = 4.5,
                    arrow = arrow(length = unit(0.005, "npc"))) +
    
    geom_text_repel(aes(label = GeneSymbol), 
                    box.padding = 1, 
                    direction="both",
                    nudge_x = condition2_nudge["x"],
                    nudge_y = condition2_nudge["y"],
                    
                    condition2_labeled_genes, 
                    alpha = 0.8,
                    color = colors[2], 
                    min.segment.length = unit(0, 'lines'),
                    size = 4.5,
                    arrow = arrow(length = unit(0.005, "npc"))) +
    
    geom_abline() +
    scale_color_manual(values = colors) +
    theme(legend.position = "none")
  
  if (parameter == "mean" | parameter == "disp") {
    gp <- gp + scale_x_log10() + scale_y_log10()
  }
  
  
  gp
}


volcanoPlot <- function(data, condition1, condition2, geneid, colors, parameter = c("mean", "disp", "resdisp")) {
  
  coefficient_colname <- switch(parameter, mean = "MeanLog2FC", disp = "DispLog2FC", resdisp = "ResDispDistance")
  prob_colname <- switch(parameter, mean = "ProbDiffMean", disp = "ProbDiffDisp", resdisp = "ProbDiffResDisp")
  diff_result_colname <- switch(parameter, mean = "ResultDiffMean", disp = "ResultDiffDisp", resdisp = "ResultDiffResDisp")
  
  colors <- c(colors, "#9d9683")
  names(colors) <- c(condition1, condition2, "NoDiff")
  
  gp <- ggplot(data, aes_string(coefficient_colname, prob_colname, color = diff_result_colname)) +
    geom_point(alpha = .7, size = .8) +
    scale_color_manual(values = colors) +
    theme(legend.position = "none")
  
  
  gp
}

compateChangeMeanAndVariability <- function(data, condition1, 
                                            condition2,
                                            geneid, colors,
                                            ngenes = 3L, 
                                            variability = c("disp", "resdisp"),
                                            flip = TRUE) {
  
  
  
  data <- copy(data)
  variability <- match.arg(variability, choices = c("disp", "resdisp"))
  
  diff_col <- ifelse(variability == "disp", "ResultDiffDisp", "ResultDiffResDisp")
  effect_col <- ifelse(variability == "disp", "DispLog2FC", "ResDispDistance")
  
  # print(head(data[[effect_col]]))
  
  if (flip) {
    data.table::set(x = data, j = "MeanLog2FC", value = - data[["MeanLog2FC"]])
    data.table::set(x = data, j = effect_col, value = - data[[effect_col]])
  }
  
  # print(head(data[[effect_col]]))
  
  condition11_labeled_genes <- data[data[[diff_col]] == condition1][MeanLog2FC > 0.4] %>%
    .[order(abs(MeanLog2FC), decreasing = TRUE)] %>%
    .[seq_len(ngenes)]
  
  condition12_labeled_genes <- data[data[[diff_col]] == condition1][MeanLog2FC < -0.4] %>%
    .[order(abs(MeanLog2FC), decreasing = TRUE)] %>%
    .[seq_len(ngenes)]
  
  condition13_labeled_genes <- data[data[[diff_col]] == condition1][abs(MeanLog2FC) < 0.4] %>%
    .[order(abs(get(effect_col)), decreasing = TRUE)] %>%
    .[seq_len(ngenes)]
  
  condition21_labeled_genes <- data[data[[diff_col]] == condition2][MeanLog2FC > 0.4] %>%
    .[order(abs(MeanLog2FC), decreasing = TRUE)] %>%
    .[seq_len(ngenes)]
  
  condition22_labeled_genes <- data[data[[diff_col]] == condition2][MeanLog2FC < -0.4] %>%
    .[order(abs(MeanLog2FC), decreasing = TRUE)] %>%
    .[seq_len(ngenes)]
  
  condition23_labeled_genes <- data[data[[diff_col]] == condition2][abs(MeanLog2FC) < 0.4] %>%
    .[order(abs(get(effect_col)), decreasing = TRUE)] %>%
    .[seq_len(ngenes)]
  
  colors <- c(colors, "#9d9683")
  names(colors) <- c(condition1, condition2, "NoDiff")
  
  gp <- ggplot(data, aes_string("MeanLog2FC", effect_col)) +
    geom_smooth(method = "lm", color = "purple", se = FALSE, formula = y ~ x, alpha = .7) +
    
    geom_point(aes_string(color = diff_col), alpha = .7, size = .8) +
    xlab("Change in Mean Gene Expression") +
    
    geom_text_repel(aes(label = GeneSymbol), 
                    condition11_labeled_genes,
                    max.overlaps = 9999,
                    alpha = 0.8,
                    color = colors[1], 
                    min.segment.length = unit(0, 'lines'),
                    size = 4.5,
                    direction = "both") +
    
    geom_text_repel(aes(label = GeneSymbol), 
                    condition12_labeled_genes,
                    max.overlaps = 9999,
                    alpha = 0.8,
                    color = colors[1], 
                    min.segment.length = unit(0, 'lines'),
                    size = 4.5,
                    direction = "both") +
    
    geom_text_repel(aes(label = GeneSymbol), 
                    condition21_labeled_genes,
                    max.overlaps = 9999,
                    size = 4.5,
                    alpha = 0.8,
                    color = colors[2], 
                    min.segment.length = unit(0, 'lines'),
                    direction = "both") +
    
    geom_text_repel(aes(label = GeneSymbol), 
                    condition22_labeled_genes,
                    max.overlaps = 9999,
                    size = 4.5,
                    alpha = 0.8,
                    color = colors[2], 
                    min.segment.length = unit(0, 'lines'),
                    direction = "both") +
    
    geom_text_repel(aes(label = GeneSymbol), 
                    condition13_labeled_genes,
                    max.overlaps = 9999,
                    size = 4.5,
                    alpha = 0.8,
                    color = colors[1], 
                    min.segment.length = unit(0, 'lines'),
                    direction = "both") +
    
    geom_text_repel(aes(label = GeneSymbol), 
                    condition23_labeled_genes,
                    max.overlaps = 9999,
                    size = 4.5,
                    alpha = 0.8,
                    color = colors[2], 
                    min.segment.length = unit(0, 'lines'),
                    direction = "both") +
    
    
    scale_color_manual(values = colors) +
    theme(legend.position = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", size = .5) +
    geom_hline(yintercept = 0, linetype = "dashed", size = .5)
  
  if (variability == "disp") {
    gp <- gp + 
      ylab("Change in Overdispersion")
  } else {
    gp <- gp + 
      ylab("Change in Residual Overdispersion")
  }
  
  gp
  
}


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
