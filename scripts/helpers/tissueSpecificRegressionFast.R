# source("./scripts/helpers/preprocess.R")

tissueSpecificRegression <- function(counts, 
                                     ercc, 
                                     annots, 
                                     tissues, 
                                     mean_design = ~ 1, 
                                     normalization_groups = NULL, 
                                     max_iterations = 10L,
                                     wald_test = TRUE,
                                     cores = parallel::detectCores()-1) {
  
  # compute ercc technical factor
  tech_factors <- normalizationFactor(ercc, 
                                      groups = normalization_groups)
  
  # compute normalization factors for each tissue
  tissue_norm_factors <- normalizationFactorTissueSpecific(counts, 
                                                           tissues = tissues, 
                                                           groups = normalization_groups)
  tissue_norm_factors <- simplify2array(tissue_norm_factors)
  
  # compute tissue size factor
  tissue_size_factors <- sweep(tissue_norm_factors, 1, tech_factors, "/")
  
  # create model matrix
  X <- model.matrix(mean_design, annots)
  colnames(X) <- stringr::str_remove_all(colnames(X), "\\(|\\)")
  X <- X[, ! apply(X == 0, 2, all), drop = FALSE]
  
  # select apply function
  if (cores <= 1) {
    applyfun <- function(X, FUN) pbapply::pblapply(X = X, FUN = FUN)
  } else {
    applyfun <- function(X, FUN) parallel::mclapply(mc.cores = cores, X = X, FUN = FUN)
  }
  
  # fit model 
  fits <-applyfun(seq_len(nrow(counts)), function(i) {
    
    fitModel(y = counts[i, ],
             X = X,
             tech_factors = tech_factors,
             tissue_size_factors = tissue_size_factors,
             max_iterations = max_iterations,
             wald_test = wald_test)
    
  })
  
  # extract results
  convergence <- sapply(fits, function(x) x$convergence)
  names(convergence) <- rownames(counts)
  
  coefficients <- sapply(fits, function(x) x$beta)
  if (! is.matrix(coefficients)) {
    coefficients <- matrix(coefficients)
    colnames(coefficients) <- colnames(X)
  } else {
    coefficients <- t(coefficients)
  }
  rownames(coefficients) <- rownames(counts)
  
  proportions <- sapply(fits, function(x) x$prop)
  if (! is.matrix(proportions)) {
    proportions <- matrix(proportions)
    colnames(proportions) <- colnames(tissue_size_factors)
  } else {
    proportions <- t(proportions)
  }
  rownames(proportions) <- rownames(counts)
  
  delta <- sapply(fits, function(x) x$delta)
  names(delta) <- rownames(counts)
  
  if (wald_test == TRUE) {
    
    coefficients_se <- sapply(fits, function(x) x$beta_se)
    if (is.matrix(coefficients_se)) {
      coefficients_se <- t(coefficients_se)
    } else {
      coefficients_se <- matrix(coefficients_se)
    }
    rownames(coefficients_se) <- rownames(counts)
    colnames(coefficients_se) <- colnames(X)
    
    pvalue <- sapply(fits, function(x) x$pvalue)
    if (is.matrix(pvalue)) {
      pvalue <- t(pvalue)
    } else {
      pvalue <- matrix(pvalue)
    }
    rownames(pvalue) <- rownames(counts)
    colnames(pvalue) <- colnames(X)
    
    padj <- apply(pvalue, 2, function(x) p.adjust(x, "fdr"))
    rownames(padj) <- rownames(counts)
    colnames(padj) <- colnames(X)
    
  }
  
  # final normalization factors
  final_norm_factors <- sweep(proportions %*% t(tissue_size_factors), 2, tech_factors, "*")
  norm <- counts / final_norm_factors
  
  out <- list(tissue_norm_factors = tissue_norm_factors,
              tech_factors = tech_factors,
              tissue_size_factors = tissue_size_factors, 
              final_norm_factors = final_norm_factors,
              coefficients = coefficients,
              dispersions = delta,
              proportions = proportions,
              convergence = convergence,
              counts = counts,
              norm = norm)
  
  if (wald_test == TRUE) {
    out <- c(out, list(
      coefficients_se = coefficients_se,
      pvalue  = pvalue,
      padj = padj))
  }
  
  out
}

computeEta <- function(params, 
                       X,
                       log_tech_factors, 
                       tissue_size_factors) {
  
  index_beta  <- seq_len(ncol(X))
  index_q <- (ncol(X)+1):length(params)
  
  beta <- params[index_beta]
  prop <- exp(params[index_q])
  prop <- prop/sum(prop)
  
  eta <- X %*% beta + log_tech_factors + log(tissue_size_factors %*% prop)
  eta
}

computeNegativeLogLikelihood <- function(params, 
                                         y, 
                                         X,
                                         log_tech_factors, 
                                         tissue_size_factors, 
                                         theta,
                                         model) {
  
  
  eta <- computeEta(params, 
                    X,
                    log_tech_factors, 
                    tissue_size_factors)
  
  if (model == "nb") {
    ll <- dnbinom(y, mu = exp(eta), size = theta, log = TRUE)
    out <- - sum(ll)
    
  } else if (model == "poisson") {
    ll <- dpois(y, lambda = exp(eta), log = TRUE)
    out <- - sum(ll)
  }
  
  out
}

computeNegativeLogLikelihoodGradient <- function(params, 
                                                 y, 
                                                 X,
                                                 log_tech_factors, 
                                                 tissue_size_factors, 
                                                 theta,
                                                 model) {
  
  
  index_beta  <- seq_len(ncol(X))
  index_q <- (ncol(X)+1):length(params)
  
  beta <- params[index_beta]
  exp_q <- exp(params[index_q])
  prop <- exp_q/sum(exp_q)
  
  eta <- computeEta(params, 
                    X,
                    log_tech_factors, 
                    tissue_size_factors)
  exp_eta <- exp(eta)
  
  if (model == "nb") {
    dnll_deta <- (exp_eta - y) / (1 + exp_eta/theta)
  } else if (model == "poisson") {
    dnll_deta <- exp_eta - y
  }
  
  # gradient beta
  deta_dbeta <- X
  dnll_dbeta <- deta_dbeta * matrix(dnll_deta, nrow(deta_dbeta), ncol(deta_dbeta))
  grad_beta  <- matrixStats::colSums2(dnll_dbeta)
  
  # gradient q
  sum_phi_exp_q <- tissue_size_factors %*% exp_q
  deta_dq <- tissue_size_factors / matrix(sum_phi_exp_q, length(sum_phi_exp_q), ncol(tissue_size_factors))
  deta_dq <- deta_dq - 1/sum(exp_q)
  deta_dq <- deta_dq * matrix(exp_q, nrow(deta_dq), ncol(deta_dq), byrow = TRUE)
  # colnames(deta_dq) <- colnames(tissue_size_factors)
  
  # dnll_dq <- sweep(deta_dq, 1, dnll_deta, "*")
  dnll_dq <- deta_dq * matrix(dnll_deta, nrow(deta_dq), ncol(deta_dq))
  grad_q <- matrixStats::colSums2(dnll_dq)
  
  # return results
  out <- c(grad_beta, grad_q)
  out
}

initParams <- function(X, tissue_size_factors) {
  
  beta_init <- setNames(rnorm(ncol(X)), colnames(X))
  
  q_init <- log(runif(ncol(tissue_size_factors), min = .3, max = .7))
  q_init <- setNames(q_init, colnames(tissue_size_factors))
  
  c(beta_init, q_init)
}

#Fits the model described in supplementary note 3, decomposing the total expression of a gene across an animal into four parts
#1) v_s: the effect of sample size variation (that effects both mRNAs and spike-ins) 
#2) phi_st: the relative size of each tissue in individual s compared to the population mean
#3) q_st: the relative expression of gene s in tissue t, as a proportional contributor to the whole-animal expression of that gene
#4) beta: some covariate (eg environmental condition, mutation) that changes the average abundance of a gene across all tissues
#eta_g describes the mean expression (across the population of samples) of the whole-animal abundance measured for gene g.  
#eta_g is a function of v_s, phi_st, q_st and beta.  v_s and phi_st we can estimate directly from comparing mRNA to spike-in abundance.
#sigma is the overdispersion term--population variation relative to the mean--wich depends on v_s, phi_st, and q_st.
#Since our empiric data shows populations distributed according to a negative binomial distribution
#our regression identifies the values of q_st, beta, and theta that specify the negative binomial sample distribution
#Y = NB(exp(eta),sigma-1) that best match the data for each geen.

fitModel <- function(y, 
                     X, 
                     tech_factors, 
                     tissue_size_factors, 
                     max_iterations = 10,
                     wald_test = TRUE) {
  
  
  # initialize parameters
  iter <- 1
  init_converged <- FALSE
  
  while (init_converged == FALSE & iter <= max_iterations) {
    params_init <- initParams(X, tissue_size_factors)
    
    init_fit <- lbfgs::lbfgs(vars = params_init,
                             call_eval = computeNegativeLogLikelihood,
                             call_grad = computeNegativeLogLikelihoodGradient,
                             y = y,
                             X = X,
                             tissue_size_factors = tissue_size_factors,
                             log_tech_factors = log(tech_factors),
                             model = "poisson",
                             invisible = 1)
    init_converged <- init_fit$convergence >= 0
    iter <- iter + 1
  }
  
  # fit negative binomial
  theta <- theta_init <- rgamma(1, 1, 1)
  params_last <- init_fit$par
  
  epsilon <- 1e-9
  iter <- 1
  change <- Inf
  
  while(epsilon < change & iter <= max_iterations) {
    
    fit <- lbfgs::lbfgs(vars = params_last,
                        call_eval = computeNegativeLogLikelihood,
                        call_grad = computeNegativeLogLikelihoodGradient,
                        y = y,
                        X = X,
                        tissue_size_factors = tissue_size_factors,
                        log_tech_factors = log(tech_factors),
                        theta = theta,
                        model = "nb",
                        invisible = 1)
    
    eta <- computeEta(fit$par, 
                      X = X, 
                      log_tech_factors = log(tech_factors), 
                      tissue_size_factors = tissue_size_factors)
    theta <- MASS::theta.ml(y = y, mu = exp(eta))
    change <- sum( (params_last - fit$par)^2 )
    
    params_last <- fit$par
    iter <- iter + 1
  }
  if (iter > max_iterations) {
    fit$convergence <- -99
  }
  
  # format results
  fit$beta <- fit$par[seq_len(ncol(X))]
  names(fit$beta) <- colnames(X)
  
  fit$prop <- exp(fit$par[seq.int(ncol(X)+1, length(fit$par))])
  fit$prop <- fit$prop/sum(fit$prop)
  names(fit$prop) <- colnames(tissue_size_factors)
  
  fit$delta <- as.numeric(1/theta)
  fit$theta <- as.numeric(theta)
  fit$theta_se <- attr(theta, "SE")
  
  # compute hessian to get standard-error and p-values
  if (wald_test == TRUE) {
    gradient_f <- function(p) {
      computeNegativeLogLikelihood(params = p, 
                                   y = y,
                                   log_tech_factors = log(tech_factors),
                                   tissue_size_factors = tissue_size_factors,
                                   model = "nb", 
                                   theta = theta, 
                                   X = X)
    }
    
    hessian <- numDeriv::hessian(gradient_f, x = fit$par)
    
    fit$beta_se <- tryCatch(
      HelpersMG::SEfromHessian(hessian, silent = TRUE)[seq_len(ncol(X))],
      error = function(e) rep_len(NA_real_, ncol(X)))
    
    fit$pvalue <- pchisq(fit$beta^2/fit$beta_se^2, df = 1, lower.tail = FALSE)
  }
  
  fit
}
