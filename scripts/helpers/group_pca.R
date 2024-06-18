groupPCA <- function(x, groups, min_genes = 0L) {
  
  x <- x[, unlist(groups)]
  
  out <- lapply(names(groups), function(group) {
    
    genes <- groups[[group]]
    
    if (length(genes) < min_genes) {
      return(list(scores = NULL,
                  loadings = NULL, 
                  prop_var_exp = NULL))
    }
    
    subx <- scale(x[, genes], center = TRUE, scale = TRUE)
    
    s <- svd(subx, nu = 1, nv = 1)
    u <- s$u
    v <- s$v
    prop_var_exp = s$d[1]^2/sum(s$d^2)
    
    rownames(u) <- rownames(subx)
    colnames(u) <- group
    colnames(v) <- colnames(u)
    
    if (mean(apply(subx, 2, function(xx) cor(xx, u, method = "s")) < 0) < 0.5) {
      u <- - u
      v <- - v
    }
    
    list(scores = u,
         loadings = v, 
         prop_var_exp = prop_var_exp,
         genes = genes)
  })
  names(out) <- names(groups)
  out
}