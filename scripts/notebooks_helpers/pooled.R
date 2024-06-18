plotCompareSinglePooled <- function(counts, annots) {
  
  x <- counts[, annots$Sample]
  # x <- filterGenes(x, min_count = 1, min_prop = .5)
  x <- normalizeCounts(x, nf = normalizationFactor(x, groups = (annots$NumberOfWorms != 1)))
  select_single <- annots[NumberOfWorms == 1]$Sample
  select_pooled <- annots[NumberOfWorms != 1]$Sample
  
  mean_df <- data.table(GeneName = rownames(x),
                        Single = rowMeans(x[, select_single]), 
                        Pooled = rowMeans(x[, select_pooled]))
  
  ggplot(mean_df, aes(Pooled, Single)) + 
    geom_point(size = .9) + 
    geom_abline() + 
    scale_x_log10(limit = c(0.01, max(mean_df$Pooled)), labels = scales::comma_format(), breaks = 10^seq(-2, 10, by = 2)) + 
    scale_y_log10(limit = c(0.01, max(mean_df$Single)), labels = scales::comma_format(), breaks = 10^seq(-2, 10, by = 2)) +
    annotate("label", x = 6, y = 5e4, 
             label = paste("Squared Correlation\n", 
                           round(100*cor((mean_df$Single), 
                                         (mean_df$Pooled), 
                                         method = "s")**2, 3), 
                           "%")) +
    xlab("Average Counts\nPooled Worms Lysate") +
    ylab("Average Counts\nSingle Nematodes") +
    theme(legend.position = "top")
}

computeGeneCorrelation <- function(x, annots) {
  
  x <- x[, annots$Sample]
  
  keep <- rowMeans(x >= 10) > 0.9
  x <- x[keep, ]
  
  mu <- rowMeans(x)
  r2 <- apply(x, 2, cor, y = mu, method = "s")**2
  
  results <- data.table(Sample = colnames(x), 
                        Correlation = r2, 
                        SampleSize = length(mu))
  results <- merge(results, 
                   annots[, .(Sample, NumberOfWorms)], 
                   by = "Sample")
  results
}

plotGeneCorrelation <- function(biological, technical) {
  
  ggdf <- rbind(cbind(biological, Sequence = "mRNA"),
                cbind(technical, Sequence = "ERCC Synthetic RNA"))
  ggdf[, Sequence := factor(Sequence)]
  ggdf[, Sequence := relevel(Sequence, "mRNA")]
  ggdf[, NumberOfWorms := as.integer(levels(NumberOfWorms)[as.integer(NumberOfWorms)])]
  
  summ_ggdf <- ggdf[, .(Correlation = mean(Correlation),
                        SD = sd(Correlation),
                        Q25 = quantile(Correlation, 0.25),
                        Q75 = quantile(Correlation, 0.75)), 
                    by = .(NumberOfWorms, Sequence)]
  
  ggplot(ggdf, aes(NumberOfWorms, 1-Correlation, color = Sequence)) +
    geom_jitter(width = 0.85) +
    expand_limits(y = 0) + 
    # geom_line(data = summ_ggdf, size = 1) +
    scale_color_manual(name = "", values = c("gray0", "darkorange")) +
    xlab("Number Of Pooled Worms in RNA-seq Sample") +
    ylab("Distance of Sample\nto Mean Gene Expression") +
    theme(legend.position = "top")
}

volcanoPlot <- function(da, colors, ngenes = 7) {
  
  down <- da[Strain_QZ120_vs_QZ0_log2FoldChange < 0]
  down <- down[order(Strain_QZ120_vs_QZ0_padj)]$GeneName[1:ngenes]
  
  up <- da[Strain_QZ120_vs_QZ0_log2FoldChange > 0]
  up <- up[order(Strain_QZ120_vs_QZ0_padj)]$GeneName[1:ngenes]
  
  ggplot2::ggplot(da, aes(Strain_QZ120_vs_QZ0_log2FoldChange, 
                          -log10(Strain_QZ120_vs_QZ0_padj + min(Strain_QZ120_vs_QZ0_padj[Strain_QZ120_vs_QZ0_padj!=0])),
                          color = Strain_QZ120_vs_QZ0_log2FoldChange < 0)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(name = "",
                                values =  setNames(rev(colors), NULL),
                                guide = guide_legend(reverse = TRUE),
                                labels = c("Up in daf-2(e1368)",
                                           "Down in daf-2(e1368)")) +
    ggplot2::xlab("log2 Fold Change\ndaf-2(e1368) vs N2") +
    ggplot2::ylab("- log10 Adjusted p-value") +
    ggrepel::geom_text_repel(data = da[GeneName %in% up],
                             min.segment.length = unit(0, 'lines'),
                             aes(label = GeneSymbol), 
                             show.legend = FALSE) + 
    # vjust = 1,
    # hjust = 1, nudge_x = Inf, nudge_y = 0, direction = "y") +
    ggrepel::geom_text_repel(data = da[GeneName %in% down], 
                             min.segment.length = unit(0, 'lines'),
                             aes(label = GeneSymbol), 
                             show.legend = FALSE) +
    # vjust = 1,
    # hjust = 1, nudge_x = -Inf, nudge_y = 0, direction = "y") +
    ggplot2::theme(legend.position = "top")
}

plotHeatmap <- function(counts, genes, annots, row_colors, column_colors, column_title = c("daf-2(e1368)", "N2")) {
  
  # handle counts
  genes <- genes[genes %in% rownames(counts)]
  counts <- filterGenes(counts[genes, annots$Sample], droplevels(annots$Group))
  nf <- normalizationFactor(counts, annots$Group)
  norm <- normalizeCounts(counts, nf)
  log_norm <- log(norm+1)
  scaled_log_norm <- t(scale(t(log_norm)))
  
  # cluster
  row_hc <- hclust(as.dist(1 - cor(t(norm), method = "s")), method = "ward.D2")
  col_hc <- hclust(as.dist(1 - cor(norm, method = "s")), method = "ward.D2")
  
  row_dend <- as.dendrogram(row_hc)
  row_dend <- dendextend::color_branches(row_dend, k = length(row_colors), col = row_colors)
  
  col_dend <- as.dendrogram(col_hc)
  col_dend <- dendextend::color_branches(col_dend, k = length(column_colors), col = column_colors)
  
  plot(col_dend)
  
  # make heatmap
  col_fun <- circlize::colorRamp2(c(min(scaled_log_norm), 0, max(scaled_log_norm)), c("blue", "white", "red"))
  
  ComplexHeatmap::Heatmap(scaled_log_norm, 
                          use_raster = TRUE,
                          name = " ",
                          col = col_fun, 
                          
                          row_split = length(row_colors), 
                          column_split = length(column_colors), 
                          
                          row_title_gp = grid::gpar(col = row_colors), 
                          
                          column_title = names(column_colors),
                          column_title_gp = grid::gpar(col = column_colors), 
                          
                          cluster_rows = row_dend,
                          cluster_columns = col_dend,
                          show_row_names = FALSE, 
                          show_column_names = FALSE)
  
}
