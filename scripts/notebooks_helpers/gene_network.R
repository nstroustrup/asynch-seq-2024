plotMeanCor <- function(mu, rho) {
  ggdf <- reshape2::melt(rho)
  ggdf <- as.data.table(ggdf)
  colnames(ggdf) <- c("GeneName", "GeneName2", "Correlation")
  ggdf <- ggdf[GeneName != GeneName2]
  ggdf[, Mean := mu[GeneName]]
  ggdf[, AbsCorrelation := abs(Correlation)]
  
  thetas <- unique(round(10^seq(0, log10(200), length.out = 101)))
  cor <- sapply(thetas, function(x) {
    with(ggdf[Mean >= x], cor(log(Mean), Correlation^2, method = "p"))
  })
  mean_cor_df <- as.data.table(cbind(Mean = thetas, Correlation = cor))
  
  list(ggplot(ggdf, aes(Mean, Correlation^2)) +
         geom_bin2d(bins = 80) +
         scale_fill_viridis_c() +
         scale_x_log10() +
         geom_vline(xintercept = 30, linetype = "dashed", color = "darkred") +
         xlab("Mean Gene Expression") +
         ylab("Squared Correlation") +
         theme(legend.position = "none"),
       
       ggplot(mean_cor_df, aes(Mean, Correlation^2)) +
         geom_point() +
         geom_line() +
         geom_vline(xintercept = 30, linetype = "dashed", color = "darkred") +
         xlab("Mean Expression Threshold") +
         ylab("Squared Correlation Between\nMean Gene Expression and Correlation"))
}

bootstrapR2 <- function(norm1, clusters, norm2 = NULL, nboot = 500) {
  
  .formatBootECDF <- function(results, clusters, ecdf_x, ecdf_y, which) {
    df <- rbindlist(lapply(seq_along(clusters), function(i) {
      
      out <- results[[i]][[which]]
      out <- as.data.table(out)
      
      colnames(out) <- stringr::str_remove_all(colnames(out), "%")
      colnames(out) <- stringr::str_pad(colnames(out), width = 2, pad = "0")
      colnames(out) <- paste0("CI", colnames(out))
      
      out[, ECDF := ecdf_y[[i]]]
      out[, Community := i]
      out[, Rho := ecdf_x]
      out
    }))
    
    if (! is.null(names(clusters))) {
      df[, Community := names(clusters)[as.integer(Community)]]
    }
    
    df
  }
  
  .quantileFast <- function(x, probs = c(0.1, 0.9), na.rm = F) {
    
    # taken from https://gist.github.com/sikli/f1775feb9736073cefee97ec81f6b193
    
    if (na.rm) x <- x[!is.na(x)]
    
    n <- length(x)
    index <- 1 + (n - 1) * probs
    
    lo    <- floor(index)
    hi    <- ceiling(index)
    
    x  <- sort(x, partial = unique(c(lo, hi)))
    qs <- x[lo]
    
    i     <- 1:length(probs) 
    h     <- index - lo
    qs    <- (1 - h) * qs + h * x[hi]
    names(qs) <- paste0(100*probs, "%")
    qs
    
  }
  
  computeCorrelation <- function(x) {
    # cor(log(t(x)), method = "p")^2
    cor(t(x), method = "s")^2
  }
  
  # define x-axis for ECDF
  ecdf_x <- seq(0, 1, length.out = 1001)
  
  # what type of test are we doing?
  if (is.null(norm2)) {
    one_sample <- TRUE
  } else {
    one_sample <- FALSE
  }
  
  # make sure genes are everywhere
  if (one_sample) {
    all_genes <- rownames(norm1)
  } else {
    all_genes <- intersect(rownames(norm1), rownames(norm2))
    norm1 <- norm1[all_genes, ]
    norm2 <- norm2[all_genes, ]
  }
  clusters <- lapply(clusters, function(x) x[x %in% all_genes])
  
  # all genes ECDF for one-sample case
  if (one_sample) {
    r <- computeCorrelation(norm1)
    ecdf2 <- ecdf(as.numeric(r[lower.tri(r)]))(ecdf_x)
    rm(r)
  }
  
  # iterate over clusters
  results <- pblapply(seq_along(clusters), function(i) {
    
    cluster_genes <- clusters[[i]]
    
    cluster_norm1 <- norm1[cluster_genes, ]
    rho1 <- computeCorrelation(cluster_norm1)
    ecdf1 <- ecdf(as.numeric(rho1[lower.tri(rho1)]))(ecdf_x)
    
    if (! one_sample) {
      cluster_norm2 <- norm2[cluster_genes, ]
      rho2 <- computeCorrelation(cluster_norm2)
      ecdf2 <- ecdf(as.numeric(rho2[lower.tri(rho2)]))(ecdf_x)
    }
    
    # bootstrap
    boot_results <- lapply(seq_len(nboot), function(j) {
      
      boot_samples1 <- sample(ncol(cluster_norm1), replace = TRUE)
      boot_cluster_norm1 <- cluster_norm1[, boot_samples1]
      
      if (one_sample) {
        boot_genes <- sample(nrow(norm1), length(cluster_genes))
        boot_samples2 <- sample(ncol(cluster_norm1), replace = TRUE)
        boot_cluster_norm2 <- norm1[boot_genes, ][, boot_samples2]
      } else {
        boot_samples2 <- sample(ncol(cluster_norm2), replace = TRUE)
        boot_cluster_norm2 <- cluster_norm2[, boot_samples2]
      }
      
      boot_rho1 <- computeCorrelation(boot_cluster_norm1)
      boot_rho2 <- computeCorrelation(boot_cluster_norm2)
      
      boot_rho_vector1 <- as.numeric(boot_rho1[lower.tri(boot_rho1)])
      boot_rho_vector2 <- as.numeric(boot_rho2[lower.tri(boot_rho2)])
      
      boot_ecdf1 <- ecdf(boot_rho_vector1)(ecdf_x)
      boot_ecdf2 <- ecdf(boot_rho_vector2)(ecdf_x)
      
      # KS
      w <- c(boot_rho_vector1, boot_rho_vector2)
      z <- cumsum(ifelse(order(w) <= length(boot_rho_vector1), -1/length(boot_rho_vector1), 1/length(boot_rho_vector2)))
      statistic <- z[which.max(abs(z))]
      
      # Median
      # statistic <- median(boot_rho_vector2) - median(boot_rho_vector1)
      
      # Mean distance
      # w <- c(boot_rho_vector1, boot_rho_vector2)
      # z <- cumsum(ifelse(order(w) <= length(boot_rho_vector1), -1/length(boot_rho_vector1), 1/length(boot_rho_vector2)))
      # statistic <- mean(z)
      
      list(boot_ecdf1 = boot_ecdf1, boot_ecdf2 = boot_ecdf2, statistic = statistic)
      
    })
    
    statistics <- sapply(boot_results, `[[`, "statistic")
    boot_ecdf1 <- sapply(boot_results, `[[`, "boot_ecdf1")
    boot_ecdf2 <- sapply(boot_results, `[[`, "boot_ecdf2")
    
    pvalue <- mean(statistics > 0)
    pvalue <- 2*min(pvalue, 1-pvalue)
    
    boot_ecdf_quantiles1 <- t(apply(boot_ecdf1, 1, .quantileFast, c(0.01, 0.99)))
    boot_ecdf_quantiles2 <- t(apply(boot_ecdf2, 1, .quantileFast, c(0.01, 0.99)))
    
    list(ecdf1 = ecdf1,
         ecdf2 = ecdf2,
         boot_ecdf_quantiles1 = boot_ecdf_quantiles1,
         boot_ecdf_quantiles2 = boot_ecdf_quantiles2,
         pvalue = pvalue)
  })
  
  ecdf1 <- lapply(results, `[[`, "ecdf1")
  ecdf2 <- lapply(results, `[[`, "ecdf2")
  pvalue <- sapply(results, `[[`, "pvalue")
  names(pvalue) <- names(clusters)
  
  boot_ecdf_quantiles_df1 <- .formatBootECDF(results, clusters, ecdf_x, ecdf1, "boot_ecdf_quantiles1")
  boot_ecdf_quantiles_df2 <- .formatBootECDF(results, clusters, ecdf_x, ecdf2, "boot_ecdf_quantiles2")
  
  list(boot_ecdf_quantiles_df1 = boot_ecdf_quantiles_df1, 
       boot_ecdf_quantiles_df2 = boot_ecdf_quantiles_df2,
       pvalue = pvalue) 
}

plotDiffAna <- function(da, 
                        clusters, 
                        colors, 
                        alpha = 1e-3, 
                        nboot = 10000) {
  
  # format data
  da <- copy(da)
  
  clusters_vector <- rbindlist(lapply(names(clusters), function(i) 
    data.table(Community = stringr::str_pad(gsub("Community", "", i), width = 2, pad = "0"), 
               GeneName = clusters[[i]])))
  clusters_vector <- pull(clusters_vector, "Community", "GeneName")
  
  da[, Community := clusters_vector[as.character(GeneName)]]
  da[is.na(Community), Community := "None"]
  da[, Community := factor(Community)]
  
  da[, Significant := ifelse(padj < alpha, "True", "False")]
  da[, Significant := factor(Significant, c("True", "False"))]
  
  # bootstrap  
  boot <- da[, {
    
    out <- sapply(seq_len(nboot), function(i) {
      sel <- sample(.N, .N, replace = TRUE)
      mean(rnorm(.N, log2FoldChange[sel], lfcSE[sel]))
    })
    mu <- mean(out)
    sdev <- sd(out)
    pvalue <- pnorm(0, mu, sdev)
    pvalue <- 2 * min(c(pvalue, 1-pvalue))
    
    list(Mean = mu, SE = sdev, PValue = pvalue)
    
  }, by = .(Community)]
  
  # t_test <- da[, {
  #   out <- t.test(log2FoldChange)
  #   list(Mean = out$estimate, SE = out$stderr, PValue = out$p.value)
  #   
  # }, by = .(Community)]
  
  # plot
  ggplot(mapping = aes(x = Community)) + 
    
    geom_jitter(aes(y = log2FoldChange, color = Significant), 
                data = da,
                width = .15,
                size = .7, 
                alpha = .7) + 
    
    geom_errorbar(aes(ymin = Mean-SE, ymax = Mean+SE), 
                  data = boot, 
                  width = .2) + 
    
    
    geom_point(aes(y = Mean), 
               data = boot, 
               size = 1) + 
    
    geom_text(aes(x = Community, 
                  y = min(da$log2FoldChange) - .2, 
                  label = paste0(format(PValue, digits = 1))), 
              boot) +
    
    scale_color_manual(values = rev(setNames(colors, NULL)), 
                       name = paste0("p-value < ", alpha)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("log2 Fold Change")
}

compareGraphs <- function(graphs, 
                          genes, 
                          geneid,
                          communities = NULL,
                          edge_width = 0.2,
                          edge_color = c("black", "red", "blue"),
                          plot_vertex_names = TRUE, 
                          layout_algorithm = "nicely",
                          ...) {
  
  # make subgraph containing only specified genes
  subgraphs <- lapply(graphs, function(x) induced_subgraph(x, genes))
  
  # delete nodes with no degree in any graph
  delete <- table(unlist(lapply(subgraphs, function(g) {
    x <- degree(g)==0
    names(x)[x]
  })))
  delete <- names(delete)[delete > 1]
  for (i in 1:length(subgraphs)) {
    subgraphs[[i]] <- delete_vertices(subgraphs[[i]], delete)
  }
  rm(delete)
  
  # if graph contains no nodes, return nothing
  if (any(sapply(subgraphs, function(x) length(V(x))) == 0)) {
    # browser()
    return(NULL)
  }
  
  # match vertices between graphs
  if (length(subgraphs) > 1) {
    subgraphs[2:length(subgraphs)] <- lapply(subgraphs[2:length(subgraphs)], function(x) {
      idx <- match(V(x)$name, V(subgraphs[[1]])$name)
      permute.vertices(x, idx)
    })
  }
  
  # are edges shared or specific?
  inter_graph <- Reduce(intersection, subgraphs)
  E(inter_graph)$type <- "Shared"
  inter_names <- apply(as_edgelist(inter_graph), 1, paste, collapse="|")
  for (i in 1:length(subgraphs)) {
    
    g_names <- apply(as_edgelist(subgraphs[[i]]), 1, paste, collapse="|")
    shared_edge <- g_names %in% inter_names
    
    E(subgraphs[[i]])$type <- paste("Specific to", names(subgraphs)[i])
    E(subgraphs[[i]])$type[shared_edge] <- "Shared"
  }
  
  # union graph
  union_graph <- Reduce(union, subgraphs)
  type <- sapply(1:length(subgraphs), function(i) 
    get.edge.attribute(union_graph, paste0("type_", i)))
  type <- apply(type, 1, function(x) unique(x[! is.na(x)]))
  E(union_graph)$type <- type
  
  # centrality (out degree)
  for (i in 1:length(subgraphs)) {
    cent <- degree(subgraphs[[i]], mode = "out")
    vertex_attr(graph = subgraphs[[i]], name = "centrality") <- cent
  }
  
  # find communities in individual graphs
  # individual_comm <- lapply(subgraphs, cluster_walktrap, steps = 3)
  # for (i in 1:length(subgraphs)) {
  #   x <- letters[individual_comm[[i]]$membership]
  #   vertex_attr(graph = subgraphs[[i]], name = "separate_community") <- as.character(x)
  # }
  
  # find limits of centrality
  centrality_limits <- lapply(subgraphs, function(x) 
    vertex_attr(graph = x, name = "centrality"))
  centrality_limits <- range(unlist(centrality_limits))
  
  vertex_attr(graph = union_graph, name = "centrality") <- mean(centrality_limits)
  vertex_attr(graph = inter_graph, name = "centrality") <- mean(centrality_limits)
  
  # make list with all graphs
  all_subgraphs <- c(list(Union = union_graph, 
                          Intersection = inter_graph), 
                     subgraphs)
  
  # add community information
  if (! is.null(communities)) {
    for (i in 1:length(all_subgraphs)) {
      x <- communities[genes %in% V(all_subgraphs[[i]])$name]
      vertex_attr(graph = all_subgraphs[[i]], name = "community") <- as.character(x)
    }
  } else {
    for (i in 1:length(all_subgraphs)) {
      vertex_attr(graph = all_subgraphs[[i]], name = "community") <- "1"
    }
  }
  
  # create tidy graph
  .toTblGraph <- function(graph, geneid) {
    genesymbol <- geneid[match(names(V(graph)), GeneName)]$GeneSymbol
    graph <- set.vertex.attribute(graph, "name", value = genesymbol)
    tbl <- as_tbl_graph(graph)
    tbl
  }
  tidygraphs <- lapply(all_subgraphs, .toTblGraph, geneid = geneid)
  
  # layout is computed from union graph
  union_layout <- create_layout(tidygraphs$Union, 
                                layout = 'igraph',
                                algorithm = layout_algorithm)
  type_levels <- sort(unique(E(union_graph)$type))
  names(edge_color) <- type_levels
  
  # actual plot
  .plotGraph <- function(plot_gg, 
                         layout_gg, 
                         centrality_limits,
                         centrality_name = "centrality",
                         community_name = "community",
                         plot_vertex_names = TRUE,
                         layout_algorithm = "nicely") {
    
    set.seed(63)
    plot_gglayout <- create_layout(plot_gg, 
                                   layout = 'igraph', 
                                   algorithm = layout_algorithm)
    
    set.seed(63)
    layout_gglayout <- create_layout(layout_gg, 
                                     layout = 'igraph',
                                     algorithm = layout_algorithm)
    
    centrality_nbins <- 3
    centrality_breaks <- seq(0, max(centrality_limits), length.out = centrality_nbins+1)
    centrality_breaks <- unique(round(centrality_breaks))
    centrality_bins <- cut(plot_gglayout[, centrality_name], 
                           breaks = centrality_breaks, 
                           include.lowest = TRUE)
    
    layout_gglayout[, "centrality_bins"] <- centrality_bins
    layout_gglayout[, "community"] <- plot_gglayout[, community_name]
    attributes(layout_gglayout)$graph <- plot_gg
    
    gp <- ggraph(layout_gglayout) +
      coord_fixed() +
      geom_edge_link(aes(color = factor(type, levels = type_levels)),
                     edge_width = edge_width,
                     alpha = 0.8) +
      scale_edge_color_manual(name = "Edge", 
                              labels = type_levels,
                              values = edge_color, 
                              drop = FALSE)
    
    if (plot_vertex_names) {
      gp <- gp + 
        suppressWarnings(geom_node_label(aes(label = name, 
                                             fill = community),#, 
                                         # size = centrality_bins),
                                         alpha = 0.9)) +
        # scale_fill_stata(name = "Community") +
        guides(fill = FALSE) +
        suppressWarnings(scale_size_discrete(name = "Degree", 
                                             drop = FALSE,
                                             range = c(2, 4)))
      
    } else {
      gp <- gp +
        geom_node_point(aes(color = community), 
                        #size = centrality_bins),
                        alpha = 0.9) +
        # scale_color_stata(name = "Community") +
        guides(color = FALSE)# +
      # suppressWarnings(scale_size_discrete(name = "Degree", 
      #                                      drop = FALSE,
      #                                      range = c(2, 4)))
    }
    
    gp <- gp +
      theme(panel.border = element_rect(colour = "black", fill = NA))
    
    gp
  }
  
  list(
    N2 = .plotGraph(tidygraphs$N2, 
                    tidygraphs$Union, 
                    centrality_limits = centrality_limits, 
                    community_name = "community",
                    plot_vertex_names = plot_vertex_names,
                    layout_algorithm = layout_algorithm) +
      ggtitle("N2"),
    
    `daf-2` = .plotGraph(tidygraphs$`daf-2`, 
                         tidygraphs$Union, 
                         centrality_limits = centrality_limits, 
                         community_name = "community",
                         plot_vertex_names = plot_vertex_names,
                         layout_algorithm = layout_algorithm) +
      ggtitle("daf-2(e1368)"),
    
    `Union` = .plotGraph(tidygraphs$Union, 
                         tidygraphs$Union, 
                         centrality_limits = centrality_limits,
                         community_name = "community",
                         plot_vertex_names = plot_vertex_names,
                         layout_algorithm = layout_algorithm) +
      ggtitle("Union Graph"),
    
    
    `Intersection` = .plotGraph(tidygraphs$Intersection, 
                                tidygraphs$Union, 
                                centrality_limits = centrality_limits,
                                community_name = "community",
                                
                                plot_vertex_names = plot_vertex_names,
                                layout_algorithm = layout_algorithm) +
      ggtitle("Intersection Graph")
  )
  
}

compareGraphsList <- function(graphs, 
                              genesets, 
                              geneid, 
                              communities = NULL, 
                              light_graphs = NULL,
                              return_list = FALSE,
                              min_n_genes_for_names = 100,
                              ...) {
  
  if (! is.list(genesets))
    genesets <- list("NoName" = genesets)
  
  dense_graphs <- graphs
  
  plots <- pbapply::pblapply(names(genesets), function(i) {
    
    cat(paste(i, "\n"))
    
    genes <- genesets[[i]]
    if (! is.null(communities))
      comm <- communities[[i]]
    else 
      comm <- NULL
    
    if (length(genes) < min_n_genes_for_names) {
      
      gp_graph <- compareGraphs(dense_graphs, 
                                genes = genes, 
                                geneid = geneid,
                                communities = comm,
                                plot_vertex_names = TRUE, 
                                edge_width = 0.3,
                                ...)
    } else {
      
      
      if(! is.null(light_graphs)) {
        graphs <- light_graphs
      } else {
        graphs <- dense_graphs
      }
      
      gp_graph <- compareGraphs(graphs, 
                                genes = genes, 
                                geneid = geneid,
                                communities = comm,
                                plot_vertex_names = FALSE,
                                edge_width = 0.2,
                                ...)
    }
    
    if (is.null(gp_graph)) {
      return(NULL)
    }
    
    if (return_list) {
      return(gp_graph)
    } else {
      gp <- ggarrange(plotlist = gp_graph, 
                      common.legend = TRUE, 
                      align = "hv",
                      nrow = 2,
                      ncol = 2)
      gp <- annotate_figure(gp, top = gsub("Community", " Community ", i))
      return(gp)
    }
  })
  
  names(plots) <- names(genesets)
  plots <- plots[! sapply(plots, is.null)]
  plots
}

coordsWholeGraph <- function(graph, 
                             communities, 
                             radius = 30,
                             target_area = 50, 
                             layout_function = layout_randomly) {
  
  intersect_genes <- intersect(unlist(communities), V(graph)$name)
  communities <- lapply(communities, function(x) x[x %in% intersect_genes])
  
  gr <- induced_subgraph(graph, unlist(communities))
  
  coords <- lapply(seq_along(communities), function(i) {
    
    genes <- communities[[i]]
    genes <- intersect(genes, V(gr)$name)
    sg <- induced_subgraph(gr, genes)
    
    # coords <- layout_nicely(sg)
    coords <- layout_function(sg)
    rownames(coords) <- genes
    
    if (nrow(coords) > 1)
      coords <- scale(coords)
    
    area <- max(apply(coords, 1, function(x) sum(x^2)))
    coords <- sqrt(target_area) * coords / sqrt(area)
    
    alpha <-  2*pi * (i-1) / length(communities)
    coords[, 1] <- coords[, 1] + radius * sin(alpha)
    coords[, 2] <- coords[, 2] + radius * cos(alpha)
    
    coords
  })
  coords <- Reduce(rbind, coords)
  coords <- coords[V(gr)$name, ]
  coords
}

plotWholeGraph <- function(graph, 
                           rho,
                           communities, 
                           community_colors,
                           radius = 30,
                           target_area = 50,
                           coords = NULL) {
  
  intersect_genes <- intersect(unlist(communities), V(graph)$name)
  communities <- lapply(communities, function(x) x[x %in% intersect_genes])
  
  communities_vector <- unlist(lapply(seq_along(communities), function(i) {
    x <- communities[[i]]
    setNames(rep_len(i, length(x)), x)
  }))
  
  
  gr <- induced_subgraph(graph, unlist(communities))
  comm <- cluster_walktrap(gr)
  comm$membership <- communities_vector[comm$names]
  
  if (is.null(coords)) {
    coords <- coordsWholeGraph(graph, 
                               communities, 
                               radius,
                               target_area)
  }
  
  edge_color <- as.data.table(as_data_frame(gr))
  edge_color[, from_comm := communities_vector[from]]
  edge_color[, to_comm := communities_vector[to]]
  x <- apply(edge_color[, .(from, to)], 1, function(x) rho[x[1], x[2]])
  edge_color[, correlation := x]
  edge_color[, col := cut(correlation, seq(-1, 1, length.out = 301))]
  pal <- colorspace::diverge_hcl(length(levels(edge_color$col)), rev = TRUE)
  levels(edge_color$col) <- rev(pal)
  edge_color[, col := as.character(col)]
  
  edge_color[, type :=  ifelse(from_comm == to_comm, 1, 0)]
  edge_color <- edge_color[order(type, abs(correlation), decreasing = FALSE),]
  gr2 <- igraph::graph_from_data_frame(edge_color[, .(from, to)], 
                                       directed = FALSE,
                                       vertices = V(gr)$name)
  
  plot(comm, 
       gr2, 
       col = "cornsilk",
       edge.color = edge_color$col,
       edge.arrow.size = .05,
       vertex.label = NA,
       vertex.size = 2, 
       layout = coords,
       mark.border = "black", 
       mark.col = community_colors)
}

corMatToDataTable <- function(n2, daf2, communities, geneid) {
  
  communities_vector <- unlist(lapply(names(communities), function(i) {
    comm <- as.integer(gsub("[^0-9]*", "", i))
    out <- rep_len(comm, length(communities[[i]]))
    names(out) <- communities[[i]]
    out
  }))
  
  rho_n2_df <- as.data.table(reshape2::melt(n2))
  colnames(rho_n2_df) <- c("GeneName1", "GeneName2", "Correlation_N2")
  
  rho_daf2_df <- as.data.table(reshape2::melt(daf2))
  colnames(rho_daf2_df) <- c("GeneName1", "GeneName2", "Correlation_daf2")
  
  rho_df <- merge(rho_n2_df, rho_daf2_df, by = c("GeneName1", "GeneName2"))
  
  rho_df <- merge(rho_df, geneid[, .(GeneSymbol2 = GeneSymbol, GeneName2 = GeneName)], by = "GeneName2")
  rho_df <- merge(rho_df, geneid[, .(GeneSymbol1 = GeneSymbol, GeneName1 = GeneName)], by = "GeneName1")
  
  rho_df[, Community1 := factor(communities_vector[GeneName1])]
  rho_df[, Community2 := factor(communities_vector[GeneName2])]
  rho_df
}

cross <- function(communities, graph) {
  m <- communities
  el <- as_edgelist(graph, names = FALSE)
  m1 <- m[el[, 1]]
  m2 <- m[el[, 2]]
  res <- m1 != m2
  if (!is.null(names(m1))) {
    names(res) <- paste(names(m1), names(m2), sep = "|")
  }
  res
}

plotHeatmapCommunitiesAndTissues <- function(rho, clusters, tissues, beta = 2) {
  
  distfun <- function(x) {
    as.dist(1-x[rownames(x), rownames(x)])
  }
  
  intersect_genes <- intersect(rownames(rho), unlist(clusters))
  
  rho <- rho[intersect_genes, intersect_genes]
  clusters <- lapply(clusters, function(x) x[x %in% intersect_genes])
  
  clusters_vector <- unlist(lapply(names(clusters), function(i) {
    comm <- i# as.integer(gsub("[^0-9]*", "", i))
    out <- rep_len(comm, length(clusters[[i]]))
    names(out) <- clusters[[i]]
    out
  }))
  clusters_vector <- clusters_vector[rownames(rho)]
  
  rho <- rho[names(clusters_vector), names(clusters_vector)]
  rho <- sign(rho) * abs(rho)^beta
  
  split <- Hmisc::capitalize(gsub("_", " ", tissues[names(clusters_vector)]))
  
  comm_color <- colorspace::rainbow_hcl(n_distinct(clusters_vector))
  names(comm_color) <- paste0("Community", 1:length(comm_color))
  
  heatmap <- ComplexHeatmap::Heatmap(rho, 
                                     name = "Correlation",
                                     use_raster = TRUE,
                                     
                                     top_annotation = HeatmapAnnotation(
                                       Community = anno_simple(clusters_vector, col = comm_color), 
                                       annotation_name_side = "left"),
                                     
                                     left_annotation =  rowAnnotation(
                                       Community = anno_simple(clusters_vector, col = comm_color), 
                                       show_annotation_name = FALSE),
                                     
                                     col = circlize::colorRamp2(c(-1, 0, 1), 
                                                                c("blue", "white", "red")),
                                     
                                     clustering_method_rows = "ward.D2", 
                                     clustering_method_columns = "ward.D2",
                                     clustering_distance_rows = distfun,
                                     clustering_distance_columns = distfun,
                                     
                                     row_split = split,
                                     column_split = split, 
                                     
                                     cluster_row_slices = FALSE, 
                                     cluster_column_slices = FALSE,
                                     
                                     show_row_names = FALSE, 
                                     show_column_names = FALSE,
                                     
                                     border = TRUE)
  
  legend <- list(
    Legend(labels = 1:length(comm_color), 
           title = "Community", 
           type = "grid",
           legend_gp = gpar(fill = comm_color)))
  
  # ComplexHeatmap::draw(gp, annotation_legend_list = legend)
  list(heatmap = heatmap, legend = legend)
}

printEnrichmentLatexTable <- function(enrich, 
                                      file = NULL, 
                                      log10_transform = list(TRUE)[rep(1, length(enrich))], 
                                      capitalize = list(FALSE)[rep(1, length(enrich))], 
                                      alpha = list(.05)[rep(1, length(enrich))]) {
  
  enrich_tab <- rbindlist(lapply(seq_along(enrich), function(i) {
    
    rbindlist(lapply(seq_len(ncol(enrich[[i]]$qvalue)), function(j) {
      
      qv <- enrich[[i]]$qvalue[, j]
      sig_qv <- qv[qv < alpha[[i]]]
      
      if (length(sig_qv) == 0) {
        out <- "No significant enrichment"
      } else {
        
        if (log10_transform[[i]]) {
          sig_qv <- round(- log10(sig_qv), 1)
        }
        sig_qv <- sort(sig_qv, decreasing = TRUE)
        
        out <- paste0(names(sig_qv), " (", sig_qv, ")")
        if (capitalize[[i]]) {
          cap_grep <- ! grepl(pattern = "^[a-z]RNA|^mTOR", out)
          out[cap_grep] <- Hmisc::capitalize(out[cap_grep])
        }
        out <- paste(out, collapse = " \\\\ ")
        out <- paste0("\\makecell[l]{", out, "}")
      }
      
      data.table(Community = j, Dataset = names(enrich)[i], Enrichment = out)
      
    }))
    
  }))
  enrich_tab[, Dataset := factor(Dataset, names(enrich))]
  enrich_tab <- enrich_tab[order(Community, Dataset)]
  
  
  add_to_row_pos <- lapply(seq_along(enrich), function(i) {
    if (i == 1) {
      seq(length(enrich), nrow(enrich_tab)-1, by = length(enrich))
    } else {
      seq(i-1, nrow(enrich_tab), by = length(enrich))
    }
  })
  add_to_row_command <- sapply(seq_along(enrich), function(i) {
    if (i == 1) {
      "\\hline"
    } else {
      "\\cline{2-3}"
    }
  })
  
  options(xtable.sanitize.text.function=identity)
  if (is.null(file)) {
    
    print(xtable::xtable(enrich_tab, align = "l|c|l|l|"), 
          tabular.environment = "longtable",
          floating = FALSE,
          include.rownames=FALSE, 
          add.to.row = list(pos = add_to_row_pos,
                            command = add_to_row_command))
    
  } else {
    
    print(xtable::xtable(enrich_tab, align = "l|c|l|l|"), 
          tabular.environment = "longtable",
          floating = FALSE,
          include.rownames=FALSE, 
          add.to.row = list(pos = add_to_row_pos,
                            command = add_to_row_command),
          file = file)
  }
  
  invisible(NULL)
}

cat2 = function(val){
  writeLines(val,con=outf,sep="")
}
catg = function(symb,gene_name){
  if (grepl("-",symb)){
    cat2(paste0("\\textit{",symb,"}"))
  }else cat2(paste0(symb));
  if (gene_name %in% orthlist$WormBase.ID){
    cat2(paste0(" (",olist[gene_name,"HGNC.Symbol"],")"))
  }
}
gene_count_str = function(num){
  if (num == 1){ return("1 gene");
  }else return(paste0(num, " genes"))
}
print_pretty_table = function(d8_tissue_communities, enrich,orthlist,many_genes=F){
  olist = aggregate(HGNC.Symbol~WormBase.ID,data=orthlist,FUN=function(x)paste0(unique(x)[1],collapse="/"))
  
  
  tissue_tab <- sapply(names(d8_tissue_communities), function(i) {
    
    genes <- d8_tissue_communities[[i]]
    comm_tissues <- tissues[names(tissues) %in% genes]
    
    tab_comm_tissues <- table(factor(comm_tissues, levels = sort(unique(tissues))))
    names(tab_comm_tissues) <- gsub("_", " ", names(tab_comm_tissues))
    names(tab_comm_tissues) <- Hmisc::capitalize(names(tab_comm_tissues))
    tab_comm_tissues
  })
  
  rownames(olist) = olist$WormBase.ID;
  
  cat2("\\documentclass{article}\n\\usepackage{booktabs}\n\\usepackage{longtable}\n\\usepackage[margin=.8in]{geometry}\n\\begin{document}\n");

    cat2("\\begin{longtable}{|p{\\linewidth}|}\n\\toprule\n\\endhead\n")
  
  for (i in 1:nrow(d8_tissue_communities)){
    all_genes <- d8_tissue_communities[[i]]
    cat(paste0("**Co-Expression Group ",i)," (",length(all_genes),")**\n")
    cat2(paste0("\\textbf{Co-Expression Group ",i,"} (",length(all_genes)," Genes)"))
    cat2("\\\\*\n");
    cat2("{\\small ");
    if (many_genes){
      to_write = c("enrichment","genes") 
    }
    else to_write = c("genes","enrichment")  
    for (write_order in to_write){
      if (write_order == "genes"){
        for (t in c(order(tissue_tab[,i],decreasing=T))){
          t_name =stringr::str_replace(tolower(rownames(tissue_tab)[t])," ","_")
          t_label = tolower(rownames(tissue_tab)[t]);
          if (t_label == "intestine")
            t_label = "intestine/pharynx"
          cur_t = tissues[tissues == tolower(t_name)];
          tissue_specific_genes = names(cur_t)[names(cur_t) %in% all_genes]
          rownames(geneid) = geneid$GeneName;
          gene_sym = geneid[tissue_specific_genes,]$GeneSymbol
          if (length(gene_sym) > 0){
            cat(paste0("*",t_label,"*", "(",tissue_tab[t,i],"):"))
            cat2(paste0("\\textbf{\\textit{",t_label,"}}", "(",gene_count_str(tissue_tab[t,i]),"): "))
            if(length(gene_sym)>1)
              for (j in 1:(length(gene_sym)-1)){
                cat(paste0(gene_sym[j],", "))
                catg(gene_sym[j],tissue_specific_genes[j])
                cat2(", ")
              }
            cat(paste(gene_sym[[length(gene_sym)]], " "))
            catg(gene_sym[[length(gene_sym)]],tissue_specific_genes[length(gene_sym)])
            cat2(" ")
          }
        }
        
        non_tissue_specific_genes = all_genes[!(all_genes %in% names(tissues))];
        gene_sym = geneid[non_tissue_specific_genes,]$GeneSymbol
        if (length(gene_sym) > 0){
          cat2(paste0("\\textbf{\\textit{Multi-tissue}}(",length(non_tissue_specific_genes),"): "))
          if(length(gene_sym)>1)
            for (j in 1:(length(gene_sym)-1)){
              cat(paste0(gene_sym[j],", "))
              catg(gene_sym[j],non_tissue_specific_genes[j])
              cat2(", ")
              if (j%%48==0){
                cat2("} \\\\\n{\\small ");
              }
            }
          cat(paste(gene_sym[[length(gene_sym)]], " "))
          catg(gene_sym[[length(gene_sym)]],non_tissue_specific_genes[length(gene_sym)])
          cat2(" ")
        }
        #cat2("}")
      }
      else{#write enrichment
        en = enrich[["Wormcat"]][["qvalue"]][,i]
        en_s = en[en<.01]
        en_s = en_s[order(en_s,decreasing = F)]
        if (length(en_s) > 0){
          cat2("\\textit{\\textbf{Functional Enrichment}}: ")
          cat("Functional Enrichment:")
          if (length(en_s)>1){
            for (j in 1:(length(en_s)-1)){
              cat(paste0(names(en_s)[j],", "))
              if(is.na(log10(en_s[j])))
                browser()
              cat2(paste0(names(en_s)[j]," p=10\\textsuperscript{",round(log10(en_s[j]),1),"}, "))
            }
          }
          cat(paste0(names(en_s)[length(en_s)],". "))
          cat2(paste0(names(en_s)[length(en_s)]," p=10\\textsuperscript{",round(log10(en_s[length(en_s)]),1),"}. "))
        }else {cat2("\\textit{(No Functional Enrichment)}")
          cat("No Functional Enrichment")}
        cat("\n")
        cat2("\n")
        en = enrich[["TF"]][["qvalue"]][,i]
        en_s = en[en<.01]
        en_s = en_s[order(en_s,decreasing = F)]
        p_val_suffix = " p="
        if (length(en_s) > 10)
           p_val_suffix = " "
        if (length(en_s) > 0){
          cat("*Transcription Factors:*")
          cat2("\\textit{\\textbf{Transcription Factors}}: ")
          if (length(en_s)>1){
            for (j in 1:(length(en_s)-1)){
              cat(paste0(names(en_s)[j],","))
              cat2(paste0("\\textit{",names(en_s)[j],"}",p_val_suffix,"10\\textsuperscript{",round(log10(en_s[j]),0),"}, "))
            }
          }
          cat(paste0(names(en_s)[length(en_s)]," "))
          cat2(paste0("\\textit{",names(en_s)[length(en_s)],"}",p_val_suffix,"10\\textsuperscript{",round(log10(en_s[length(en_s)]),1),"}. "))
        }else {cat("No TF Enrichment"); 
          cat2("\\textit{(No TF Enrichment)}")}
        cat("\n")
        cat2("\n")
        en = enrich[["KEGG"]][["qvalue"]][,i]
        en_s = en[en<.01]
        en_s = en_s[order(en_s,decreasing = F)]
        if (length(en_s) > 0){
          cat("*KEGG pathway:*")
          cat2("\\textit{\\textbf{KEGG pathway}}: ")
          if (length(en_s)>1){
            for (j in 1:(length(en_s)-1)){
              cat(paste0(names(en_s)[j],", "))
              cat2(paste0(names(en_s)[j]," p=10\\textsuperscript{",round(log10(en_s[j]),1),"}, "))
            }
          }
          cat(paste0(names(en_s)[length(en_s)]," "))
          cat2(paste0(names(en_s)[length(en_s)]," p=10\\textsuperscript{",round(log10(en_s[length(en_s)]),1),"}. "))
        }else {cat("No kegg enrichment");
          cat2("\\textit{(No KEGG enrichment)}");
        }
      }
      cat("\n")
      #output spacer in between genes and enrichment
      if (which(write_order %in% to_write) == 1)
        cat2("}\\hspace{.25in} {\\small")
    }
    cat("\n")
    #if (sep_each_community){
   #   cat2("}");
   # }
    cat2("}\\\\\n\\midrule\n");
    #if (sep_each_community){
    #  cat2("\\end{tabular}\\end{tabular}\\\\\n\\begin{tabular}{|p{\\linewidth}|}\n\\begin{tabular}{p{.98\\linewidth}}\n\\toprule\n")
    #}
  }
  #if(sep_each_community){
    cat2("\\end{longtable}");
  #}else{
   #cat2("\\end{tabular}\\end{tabular}\\end{table}");
 # }
  cat2("\\end{document}");
}


