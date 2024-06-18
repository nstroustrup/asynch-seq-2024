doSim <- function(nsamples = 1000, 
                  ngenes = 100, 
                  Sigma_omega = 1, 
                  Sigma_epsilon = diag(ngenes),
                  mu_mean = 10,
                  mu_sd = 1) {
  
  
  # proportion of each tissue
  p <- sort(rep_len(0:1, ngenes), decreasing = TRUE)
  
  # simulate
  mu <- rnorm(ngenes, mean = mu_mean, sd = mu_sd)
  omega <- MASS::mvrnorm(n = nsamples, mu = rep_len(0, 2), Sigma = Sigma_omega)
  epsilon <- t(MASS::mvrnorm(n = nsamples, mu = rep_len(0, ngenes), Sigma = Sigma_epsilon))
  
  log_y <- sapply(seq_len(nsamples), function(j) { # genes are rows, samples are columns
    sapply(seq_len(ngenes), function(i) {
      mu[i] + (p[i] * omega[j, 1] + (1-p[i]) * omega[j, 2]) + epsilon[i, j]
    })
  })
  y <- round(exp(log_y))
  rownames(y) <- paste0("Gene", 1:nrow(y))
  colnames(y) <- paste0("Sample", 1:ncol(y))
  
  # k_hat1 <- DESeq2::estimateSizeFactorsForMatrix(y+1)
  # k_hat2 <-cbind(DESeq2::estimateSizeFactorsForMatrix(y[p == 1, ]+1),
  #                DESeq2::estimateSizeFactorsForMatrix(y[p == 0, ]+1))
  k_hat1 <- scran::calculateSumFactors(y)
  k_hat2 <-cbind(scran::calculateSumFactors(y[p == 1, ]),
                 scran::calculateSumFactors(y[p == 0, ]))
  
  # normalize
  norm1 <- sweep(y, 2, k_hat1, "/")
  norm2 <- t(sapply(seq_len(ngenes), function(i) {
    y[i, ] / (p[i] * k_hat2[, 1] + (1-p[i]) * k_hat2[, 2])
  }))
  rownames(norm2) <- paste0("Gene", 1:nrow(y))
  
  list(y = y, 
       norm1 = norm1,
       norm2 = norm2,
       mu = mu,
       omega = omega, 
       epsilon = epsilon,
       k_hat1 = k_hat1, 
       k_hat2 = k_hat2,
       nsamples = nsamples, 
       ngenes = ngenes, 
       p = p,
       Sigma_omega = Sigma_omega, 
       Sigma_epsilon = Sigma_epsilon,
       mu_mean = mu_mean,
       mu_sd = mu_sd)
}
