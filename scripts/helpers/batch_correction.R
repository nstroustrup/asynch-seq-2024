matchQuantiles <- function(counts, old_mu, new_mu, disp) {
  size <- 1/disp
  out <- sapply(1:ncol(counts), function(i) {
    tmp_p <- pnbinom(counts[, i], mu = old_mu[, i], size = size[i], log.p = TRUE)
    qnbinom(tmp_p, mu = new_mu[, i], size = size[i], log.p = TRUE)
  })
  rownames(out) <- rownames(counts)
  colnames(out) <- colnames(counts)
  out
}

computeRegressionMean <- function(counts, annots, nf, coefficients, dispersions, full_design, biological_design) {
  
  # model matrices
  x <-  model.matrix(full_design, annots)
  bio_x <-  model.matrix(biological_design, annots)
  tech_x <- x[, ! colnames(x) %in% colnames(bio_x), drop = FALSE]
 # browser()
  # regression coefficients
  colnames(coefficients) <- colnames(x)
  bio_coefficients <- coefficients[, colnames(bio_x), drop = FALSE]
  tech_coefficients <- coefficients[, colnames(tech_x), drop = FALSE]
  
  # compute log means without normalization factor
  full_mu <- x %*% t(coefficients)
  biological_mu <- bio_x %*% t(bio_coefficients)
  tech_mu <- tech_x %*% t(tech_coefficients)
  
  # add in normalization factor
  if (is.matrix(nf)) {
    log_nf <- log(nf)
  } else if (is.vector(nf)) {
    log_nf <- matrix(log(nf), ncol(counts), 1) %*% matrix(1, 1, nrow(counts))
  } else {
    stop("Invalid nf argument.")
  }

  full_mu <- full_mu + log_nf
  biological_mu <- biological_mu + log_nf
  tech_mu <- tech_mu + log_nf
  
  # exponentiation
  full_mu <- exp(full_mu)
  biological_mu <- exp(biological_mu)
  tech_mu <- exp(tech_mu)
  
  # print covariates
  cat("All covariates\n")
  print(colnames(coefficients))
  
  cat("Biological effect covariates\n")
  print(colnames(bio_coefficients))
  
  # output
  list(full = full_mu, biological = biological_mu, technical = tech_mu)
}


correctCountsQuantileMatrix <- function(counts, annots, nf, coefficients, dispersions, full_design, biological_design) {
  
  
  means <- computeRegressionMean(counts, annots, nf, coefficients, dispersions, full_design, biological_design)
  
  # match quantiles to get corrected counts
  out <- t(matchQuantiles(counts = t(counts), 
                          old_mu = means$full, 
                          disp = dispersions, 
                          new_mu = means$biological))
  
  # handle infinity
  isinf <- is.infinite(out)
  if (any(isinf)) {
    warnings(paste0(sum(isinf), " infinite values\n"))
    out_scale <- correctCountsScaleMatrix(counts, annots, nf, coefficients, dispersions, full_design, biological_design)
    out[isinf] <- floor(out_scale[isinf])
  }
  
  out 
}

correctCountsScaleMatrix <- function(counts, annots, nf, coefficients, dispersions, full_design, biological_design) {
  
  means <- computeRegressionMean(counts, annots, nf, coefficients, dispersions, full_design, biological_design)
  
  # scale counts by batch effect
  out <- sapply(1:ncol(counts), function(j) {
    as.double(counts[, j]) / as.double(means$technical[j, ])
  })
  rownames(out) <- rownames(counts)
  colnames(out) <- colnames(counts)
  
  out
}


# for DESeq2
correctCountsQuantile <- function(dds, biological_design) {
  correctCountsQuantileMatrix(counts = DESeq2::counts(dds), 
                              annots = as.data.table(colData(dds)), 
                              nf = setNames(DESeq2::sizeFactors(dds), colnames(dds)), 
                              coefficients = coef(dds) * log(2), 
                              dispersions = DESeq2::dispersions(dds),
                              full_design = dds@design, 
                              biological_design = biological_design)
}

correctCountsScale <- function(dds, samples, biological_design) {
  correctCountsScaleMatrix(counts = DESeq2::counts(dds), 
                           annots = as.data.table(colData(dds)), 
                           nf = setNames(DESeq2::sizeFactors(dds), colnames(dds)), 
                           coefficients = coef(dds) * log(2), 
                           dispersions = DESeq2::dispersions(dds),
                           full_design = dds@design, 
                           biological_design = biological_design)
}
