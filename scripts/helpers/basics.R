# can handle chains with different number of genes
BASiCS_CorrectOffset2 <- function(Chain, ChainRef, min.mean = 1) {
  
  genes <- intersect(colnames(Chain@parameters$mu), colnames(ChainRef@parameters$mu))
  
  mu1 <- Chain@parameters$mu[, genes]
  mu2 <- ChainRef@parameters$mu[, genes]
  
  OffsetChain0 <- matrixStats::rowMedians(mu1)/matrixStats::rowMedians(mu2)
  OffsetEst0 <- median(OffsetChain0)
  OffsetRatio <- (colMedians(mu1)/OffsetEst0 + colMedians(mu2))/2
  
  include <- which(OffsetRatio >= min.mean)
  
  OffsetChain <- matrixStats::rowMedians(mu1[, include])/matrixStats::rowMedians(mu2[, include])
  OffsetEst <- median(OffsetChain)
  
  Chain_offset <- Chain
  Chain_offset@parameters$mu <- Chain@parameters$mu/OffsetEst
  
  if ("phi" %in% names(Chain@parameters)) {
    Chain_offset@parameters$phi <- Chain@parameters$phi * OffsetEst
  }
  else {
    Chain_offset@parameters$s <- Chain@parameters$s * OffsetEst
  }
  
  list(Chain = Chain_offset, Offset = OffsetEst, OffsetChain = OffsetChain)
}

# don't exclude genes for dispersion variability analysis
BASiCS_TestDE2 <- function(Chain1,
                           Chain2,
                           EpsilonM = log2(1.5),
                           EpsilonD = log2(1.5),
                           EpsilonR = log2(1.5) / log2(exp(1)),
                           ProbThresholdM = 2 / 3,
                           ProbThresholdD = 2 / 3,
                           ProbThresholdR = 2 / 3,
                           OrderVariable = c("GeneIndex", "GeneName", "Mu"),
                           GroupLabel1 = "Group1",
                           GroupLabel2 = "Group2",
                           Plot = TRUE,
                           PlotOffset = TRUE,
                           PlotOffsetType = c(
                             "offset estimate", 
                             "before-after",
                             "MA plot"
                           ),
                           Offset = TRUE,
                           EFDR_M = 0.05,
                           EFDR_D = 0.05,
                           EFDR_R = 0.05,
                           GenesSelect = rep(
                             TRUE,
                             ncol(Chain1@parameters[["mu"]])
                           ),
                           min.mean = 1,
                           MinESS = 100,
                           ...) {
  
  OrderVariable <- match.arg(OrderVariable)
  
  BASiCS:::HiddenHeaderTest_DE(
    Chain1 = Chain1,
    Chain2 = Chain2,
    EpsilonM = EpsilonM,
    EpsilonD = EpsilonD,
    EpsilonR = EpsilonR,
    EFDR_M = EFDR_M,
    EFDR_D = EFDR_D,
    EFDR_R = EFDR_R,
    ProbThresholdM = ProbThresholdM,
    ProbThresholdD = ProbThresholdD,
    ProbThresholdR = ProbThresholdR,
    OrderVariable = OrderVariable,
    GroupLabel1 = GroupLabel1,
    GroupLabel2 = GroupLabel2,
    GenesSelect = GenesSelect,
    Plot = Plot,
    PlotOffset = PlotOffset,
    Offset = Offset
  )
  
  GeneName <- colnames(Chain1@parameters$mu)
  GeneIndex <- seq_len(length(GeneName))
  
  message(
    "-------------------------------------------------------------\n",
    "Log-fold change thresholds are now set in a log2 scale. \n",
    "Original BASiCS release used a natural logarithm scale."
  )
  
  if (xor(
    is.null(Chain1@parameters[["epsilon"]]),
    is.null(Chain2@parameters[["epsilon"]]))
  ) {
    
    stop("Both chains should be run with the same setting for Regression.")
  }
  
  n1 <- ncol(Chain1@parameters$nu)
  n2 <- ncol(Chain2@parameters$nu)
  n <- n1 + n2
  
  IncludeEpsilon <- !is.null(Chain1@parameters[["epsilon"]])
  # With offset correction
  if (Offset) {
    A <- BASiCS_CorrectOffset(Chain1, Chain2, min.mean = min.mean)
    OffsetEst <- A$Offset
    Chain1_offset <- A$Chain
    Chain2_offset <- Chain2  # Chain2 requires no change
    
    # Post-offset correction gene-specific estimates
    Mu1 <- matrixStats::colMedians(Chain1_offset@parameters$mu)
    Mu2 <- matrixStats::colMedians(Chain2_offset@parameters$mu)
    Delta1 <- matrixStats::colMedians(Chain1_offset@parameters$delta)
    Delta2 <- matrixStats::colMedians(Chain2_offset@parameters$delta)
    
    # Pre-offset correction LFC estimates
    Mu1_old <- matrixStats::colMedians(Chain1@parameters$mu)
    MuBase_old <- (Mu1_old * n1 + Mu2 * n2) / n
    ChainTau_old <- log2(Chain1@parameters$mu / Chain2@parameters$mu)
    MedianTau_old <- matrixStats::colMedians(ChainTau_old)
    
    # Offset corrected LFC estimates
    MuBase <- (Mu1 * n1 + Mu2 * n2) / n
    ChainTau <- log2(Chain1_offset@parameters$mu / Chain2_offset@parameters$mu)
    MedianTau <- matrixStats::colMedians(ChainTau)
    
    if (!PlotOffset) {
      message(
        "-------------------------------------------------------------\n",
        "Offset estimate: ", round(OffsetEst, 4), "\n",
        "(ratio ", GroupLabel1, " vs ", GroupLabel2, ").\n",
        "To visualise its effect, please use 'PlotOffset = TRUE'.\n",
        "-------------------------------------------------------------\n"
      )
    } else {
      message(
        "-------------------------------------------------------------\n",
        "Offset estimate: ", round(OffsetEst, 4), "\n",
        "(ratio ", GroupLabel1, " vs ", GroupLabel2, ").\n",
        "-------------------------------------------------------------\n"
      )
    }
    
  } else {
    message(
      "-------------------------------------------------------------\n",
      "It is recomended to perform a global offset correction \n",
      "to remove global changes between the two groups of cells \n",
      "Default offset value set equal to 1.\n",
      "To perform offset correction, please set 'Offset = TRUE'. \n",
      "-------------------------------------------------------------\n"
    )
    Chain1_offset <- Chain1
    Chain2_offset <- Chain2
    
    # Default values when no offset correction is applied
    OffsetEst <- 1
  }
  
  Mu1 <- matrixStats::colMedians(Chain1_offset@parameters$mu)
  Mu2 <- matrixStats::colMedians(Chain2_offset@parameters$mu)
  Delta1 <- matrixStats::colMedians(Chain1_offset@parameters$delta)
  Delta2 <- matrixStats::colMedians(Chain2_offset@parameters$delta)
  ChainTau <- log2(Chain1_offset@parameters$mu / Chain2_offset@parameters$mu)
  
  MuBase <- (Mu1 * n1 + Mu2 * n2) / n
  
  orderVar <- switch(
    OrderVariable,
    "GeneIndex" = order(GeneIndex, decreasing = FALSE),
    "GeneName" = order(GeneName, decreasing = TRUE),
    "Mu" = order(as.numeric(MuBase), decreasing = TRUE)
  )
  q <- length(GenesSelect)
  
  GoodESS <- BASiCS:::.CheckESS(Chain1, Chain2, MinESS, "mu", q)
  ResM <- BASiCS:::.RunTest(
    Chain = ChainTau,
    Epsilon = EpsilonM,
    ProbThreshold = ProbThresholdM,
    EFDR = EFDR_M,
    Task = "Differential mean",
    Suffix = "M",
    GroupLabel1 = GroupLabel1,
    GroupLabel2 = GroupLabel2,
    GenesSelect = GenesSelect,
    Param1 = Mu1,
    Param2 = Mu2,
    n1 = n1,
    n2 = n2,
    GeneName = GeneName,
    GoodESS = GoodESS,
    Measure = "Mean"
  )
  
  # Genes with no change in mean expression
  # DE <- !(ResM@Table$ResultDiffMean %in% c("NoDiff", "ExcludedByUser"))
  
  GoodESS <- BASiCS:::.CheckESS(Chain1, Chain2, MinESS, "delta", q)
  
  ChainOmega <- log2(Chain1@parameters$delta / Chain2@parameters$delta)
  ResD <- BASiCS:::.RunTest(
    Chain = ChainOmega,
    Epsilon = EpsilonD,
    ProbThreshold = ProbThresholdD,
    EFDR = EFDR_D,
    Task = "Differential dispersion",
    Suffix = "D",
    GroupLabel1 = GroupLabel1,
    GroupLabel2 = GroupLabel2,
    GenesSelect = GenesSelect,
    Param1 = Delta1,
    Param2 = Delta2,
    n1 = n1,
    n2 = n2,
    GeneName = GeneName,
    Measure = "Disp",
    GoodESS = GoodESS
    # Excluded = DE
  )
  
  # Changes in residual over-dispersion - if regression approach was used
  if (IncludeEpsilon) {
    
    Excluded <- is.na(Chain1@parameters$epsilon[1, ]) |
      is.na(Chain2@parameters$epsilon[1, ])
    
    ChainPsi <- Chain1@parameters$epsilon - Chain2@parameters$epsilon
    Epsilon1 <- matrixStats::colMedians(Chain1@parameters$epsilon)
    Epsilon2 <- matrixStats::colMedians(Chain2@parameters$epsilon)
    
    GoodESS <- BASiCS:::.CheckESS(Chain1, Chain2, MinESS, "epsilon", q)
    ResR <- BASiCS:::.RunTest(
      Chain = ChainPsi,
      Epsilon = EpsilonR,
      ProbThreshold = ProbThresholdR,
      EFDR = EFDR_R,
      Task = "Differential residual dispersion",
      Suffix = "R",
      GroupLabel1 = GroupLabel1,
      GroupLabel2 = GroupLabel2,
      GenesSelect = GenesSelect,
      Param1 = Epsilon1,
      Param2 = Epsilon2,
      n1 = n1,
      n2 = n2,
      GeneName = GeneName,
      Measure = "ResDisp",
      GoodESS = GoodESS,
      Excluded = Excluded
    )
    ResR@Table[, "MeanOverall"] <- ResM@Table[, "MeanOverall"]
    ResR <- ResR[orderVar, ]
  }
  ResM <- ResM[orderVar, ]
  ResD <- ResD[orderVar, ]
  
  Results <- list(
    Mean = ResM,
    Disp = ResD
  )
  
  if (IncludeEpsilon) {
    Results <- c(
      Results,
      ResDisp = ResR
    )
  }
  
  Out <- new("BASiCS_ResultsDE",
             Results = Results,
             Chain1 = Chain1_offset,
             Chain2 = Chain2_offset,
             GroupLabel1 = GroupLabel1,
             GroupLabel2 = GroupLabel2,
             Offset = OffsetEst,
             RowData = DataFrame(GeneName = GeneName),
             Extras = list()
  )
  if (Plot) {
    Out@Extras[["Plots"]] <- BASiCS_PlotDE(Out)
  }
  Out
}


.DiffRes <- function(ResultDE) {
  ResultDE[.WhichDiffRes(ResultDE), ]
}

.WhichDiffRes <- function(ResultDE) {
  !ResultDE@Table$Result %in% c(
    "NoDiff",
    "ExcludedLowESS",
    "ExcludedFromTesting",
    "ExcludedByUser"
  )
}

summarizeChain <- function(chain) {
  data.table::data.table(GeneName = colnames(chain@parameters$mu),
                         Mean    = matrixStats::colMedians(chain@parameters$mu),
                         Disp    = matrixStats::colMedians(chain@parameters$delta),
                         ResDisp = matrixStats::colMedians(chain@parameters$epsilon))
}

additiveModel <- function(chain0, chain1, chain2) {
  
  # browser()
  
  effect0 <- list(mu      = log(chain0@parameters$mu),
                  delta   = log(chain0@parameters$delta),
                  epsilon = chain0@parameters$epsilon)
  
  effect1 <- list(mu      = log(chain1@parameters$mu) - log(chain0@parameters$mu),
                  delta   = log(chain1@parameters$delta) - log(chain0@parameters$delta),
                  epsilon = chain1@parameters$epsilon - chain0@parameters$epsilon)
  
  effect2 <- list(mu      = log(chain2@parameters$mu) - log(chain0@parameters$mu),
                  delta   = log(chain2@parameters$delta) - log(chain0@parameters$delta),
                  epsilon = chain2@parameters$epsilon - chain0@parameters$epsilon)
  
  add_model <- copy(chain0)
  
  add_model@parameters$mu <- exp(effect0$mu + effect1$mu + effect2$mu)
  add_model@parameters$delta <- exp(effect0$delta + effect1$delta + effect2$delta)
  add_model@parameters$epsilon <- effect0$epsilon + effect1$epsilon + effect2$epsilon
  add_model
}
