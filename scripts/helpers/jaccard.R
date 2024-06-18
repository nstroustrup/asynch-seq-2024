jaccardIndex <- function(x, y) {
  length(intersect(x, y)) / length(union(x, y))
}

computeJaccardMatrix <- function(genesets) {
  
  jaccard <- matrix(NA_real_, length(genesets), length(genesets))
  rownames(jaccard) <- names(genesets)
  colnames(jaccard) <- names(genesets)
  
  diag(jaccard) <- 1
  for (i in seq_len(ncol(jaccard)-1)) {
    for (j in seq_len(ncol(jaccard))) {
      jaccard[i, j] <- jaccardIndex(genesets[[i]], genesets[[j]])
    }
  }
  
  jaccard[lower.tri(jaccard)] = t(jaccard)[lower.tri(jaccard)]
  jaccard
}
