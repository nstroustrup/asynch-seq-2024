fisherEnrichment <- function(annotations,
                             genesets,
                             key,
                             column = "GeneName",
                             apply_function = pbapply::pbsapply,
                             universe = NULL) {
    
  if (is.null(universe))
    universe <- sort(unique(annotations[[column]]))
  # universe <- unlist(genesets, use.names = FALSE)
  annotations <- annotations[annotations[[column]] %in% universe]
  
  # remove annotations with not enough genes
  tab <- table(annotations[[key]])
  uniq_annots <- sort(names(tab)[tab >= 5])
  annotations <- annotations[annotations[[key]] %in% uniq_annots]
  #browser()
  setkeyv(annotations, key)
  pvalue <- apply_function(uniq_annots, function(x) {
    annot_genes <- annotations[x][[column]]
    sapply(names(genesets), function(y) {
      annot   <- universe %in% annot_genes
      hit <- universe %in% genesets[[y]]
      
      a <- sum(annot & hit)
      b <- sum(annot & ! hit)
      c <- sum(! annot & hit)
      d <- sum(! annot & ! hit)
      m <- matrix(c(a,b,c,d), 2, 2)
      test <- fisher.test(m, alternative = "greater")
     # browser()
      # lfc <- log10(as.numeric(test$estimate))
      pv <- as.numeric(test$p.value)
    })
  })
  if (length(genesets) > 1) pvalue <- t(pvalue)
  else {
    pvalue <- matrix(pvalue)
    colnames(pvalue) <- names(genesets)
  }
  rownames(pvalue) <- uniq_annots
  qvalue <- apply(pvalue, 2, p.adjust, method = "BH")
  
  if (! is.matrix(qvalue)) {
    qvalue <- matrix(qvalue, nrow = nrow(pvalue), ncol = ncol(pvalue))
    rownames(qvalue) <- rownames(pvalue)
    colnames(qvalue) <- colnames(pvalue)
  }
  
  list(pvalue = pvalue, qvalue = qvalue)
}

plotEnrichment <- function(enrich, alpha = 0.05, nmax = 30, adjusted = TRUE,text_size=NULL) {
  if (adjusted == TRUE) pv <- enrich$qvalue
  else pv <- enrich$pvalue
  
  gp <- lapply(colnames(pv), function(comm) {
    ggdf <- data.table(Annotation = rownames(pv), 
                       QValue = pv[, comm] + .Machine$double.xmin)
    ggdf <- ggdf[QValue < alpha][order(QValue, decreasing = FALSE)]
    if (nrow(ggdf) == 0) return(NULL)
    
    ggdf <- ggdf[1:min(c(nrow(ggdf), nmax))]
    ggdf[, Annotation := factor(Annotation, rev(unique(Annotation)))]
    
    # levels(ggdf$Annotation) <- stringr::str_wrap(levels(ggdf$Annotation), 50)
    
    gp = ggplot(ggdf, aes(-log10(QValue), Annotation),cex=cex) +#, fill = -log10(QValue))) +
      geom_point(size = 2) + 
      expand_limits(x = 0) +
      # geom_vline(xintercept = 1) + 
      # geom_bar(stat = "identity") +
      # scale_fill_viridis_c() +
      xlab("-log10 adjusted p-value") +
      ylab("") +
      ggtitle(comm)
      if (!is.null(text_size)){
      	gp = gp + theme(text = element_text(size = text_size),axis.text = element_text(size = text_size),legend.position = "none") 
      }else  gp = gp + theme(legend.position = "none")
      gp
  })
  names(gp) <- colnames(pv)
  gp <- gp[! sapply(gp, is.null)]
  gp
}
