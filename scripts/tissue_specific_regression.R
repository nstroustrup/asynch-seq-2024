# fit a model that estimates the effect of some covariate on whole-animal gene expression
# assuming that this expression is influenced by variability in tissue sizes across individuals
# that generates variability in the expression of all genes expressed in these tissues.
# Important: the covariate is assumed not to effect the variance of gene expression
# so DO NOT use age as a covariate as age will almost certainly break this assumption.
# Usefully, the theta (overdispersion) parameter estimated by this model will be 
# the "underlying per-gene population variability" left over after accounting for 
# tissue-size differences across the population.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(scran)
  library(stringr)
  library(optparse)
  library(pbapply)
  library(lbfgs)
})

run_defaults = F

if (!run_defaults){
setwd("/users/nstroustrup/nstroustrup/projects/asynch-seq-2022")
}

source("./scripts/helpers/preprocess.R")
source("./scripts/helpers/tissueSpecificRegressionFast.R")
source("./scripts/helpers/batch_correction.R")
pbo = pboptions(type="txt")

opt <- list(input = "./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds",
            name = "single_and_pooled_worms_2",
            output_dir = "./data/tissue_specific_regression/tables/single_and_pooled_worms/",
            models = "./data/tissue_specific_regression/models/single_and_pooled_worms.csv",
            specific_model = c("n2_d8_sw"),
            core_count = 7,
	    exclude_ts = T)

if(1)
opt <- list(input = "./data/formated_counts/counts_and_annots_net_validation.rds",
            name = "net_validation",
	    specific_model="",
            output_dir = "./data/tissue_specific_regression/tables/net_validation_2/",
            models = "./data/tissue_specific_regression/models/net_validation.csv",
	    exclude_ts=T,
	    core_count = 1)

# handle arguments
option_list = list(
  
  make_option(c("-n", "--name"), 
              type="character",
              default=NULL, 
              help="Name of experiment."),
  
  make_option(c("-i", "--input"), 
              type="character",
              default=NULL, 
              help="Input formated counts RDS file."),
  
  make_option(c("-m", "--models"), 
              type="character",
              default=NULL, 
              help="Input models CSV file."),

make_option(c("-S", "--output_suffix"),
              type="character",
              default="",
              help="suffix_for_output_files"),              

  make_option(c("-s", "--specific_model"), 
		    type="character",
		    default=NULL, 
              help="Run only a specific model"),
  
  make_option(c("-o", "--output_dir"), 
              type="character",
              default=NULL,
              help="Output directory for result tables."),
  
  make_option(c("-r", "--output_rds"), 
              type="character",
              default=NULL,
              help="Output file for list RDS object."),

make_option(c("-e","--exclude_ts"),
	type="logical",
	default=F,
help="Exclude models that include only tissue-specific genes"),
              
    make_option(c("-c", "--core_count"), 
                type="character",
                default=parallel::detectCores()-1,
              help="How many cores should be used for parallel computation"))
if (!run_defaults){
opt_parser <- OptionParser(option_list = option_list, prog = "Differential Analysis")
opt <- parse_args(opt_parser)
}
if (is.null(opt$name)) {
  stop("Missing name of experiment (--name)")
}

if (is.null(opt$input)) {
  stop("Missing formated counts RDS file (--input)")
}

if (is.null(opt$models)) {
  stop("Missing models CSV file (--models)")
}

# load data
cat("Loading data\n")

data <- readRDS(opt$input)

models <- fread(opt$models)
models <- models[Run == TRUE]
if (opt$exclude_ts == T){
   if (!"TissueSpecificNormalization" %in% names(models)){
   stop("The column TissueSpecificNormalization is not present in the models file.")
   }
   cat("Excluding tissue-specific gene models\n")
   models = models[TissueSpecificNormalization == F,]
}
if (!is.null(opt$specific_model) & opt$specific_model != ""){
	if (!opt$specific_model %in% models$Name)
		stop(paste("Specified model was not found:",opt$specific_model))
	models <- models[Name == opt$specific_model]
}

if (any(duplicated(models$Name))) {
  stop("Duplicated names in model file!")
}

tissues_df <- fread("./data/tissue_unique_genes/genes_unique_to_tissues.csv.gz")
tissues <- pull(tissues_df, "Tissue", "GeneName")
tissues[tissues != "germ_line"] <- "soma"
uniq_tissues <- sort(unique(tissues))

# create directories
if (! is.null(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (! is.null(opt$output_rds)) {
  dir.create(dirname(opt$output_rds), recursive = TRUE, showWarnings = FALSE)
}
cat(paste(nrow(models),"models\n"))

# run analysis
cat("Started analyses\n")
results <- lapply(seq_len(nrow(models)), function(i) {
  
  # model
  model <- models[i, ]

  cat("-------\n")
  cat(paste("#### Running", i, "out of", nrow(models)))
  cat("\nModel\n")
  print(model)
  cat("\n")
  
  # get annotations
  annots <- copy(data$annots[eval(parse(text = model$Subset))])
  
  # drop levels
  for (j in colnames(annots)) {
    value <- annots[[j]]
    if (is.factor(value)) {
      set(annots, j = j, value = droplevels(value))
    }
  }
  
  if (nrow(annots) == 0) return(NULL)
  
  # handle specific case
  if (opt$name == "tissue_specific_rnai" | opt$name == "net_validation") {
    annots[, RNAi := relevel(RNAi, "EV")]
  }
  if (opt$name == "glp1_glp4") {
    levels(annots$Genotype) <- c("EV", "glp1", "glp4")
  }
  if (opt$name == "single_and_pooled_worms") {
    annots[, Pooled := as.integer(NumberOfWorms != 1)]
  }
  
  # make formula
  if (model$BiologicalCovariates != "NULL") {
    full_design <- model$BiologicalCovariates
    biological_design <- as.formula(paste0("~", model$BiologicalCovariates))
  } else {
    full_design <- NULL
    biological_design <- NULL
  }
  
  if (model$TechnicalCovariates != "NULL") {
    technical_design <- as.formula(paste0("~", model$TechnicalCovariates))
    
    if (! is.null(full_design)) {
      full_design <- paste(full_design, "+", model$TechnicalCovariates)
    } else {
      full_design <- model$TechnicalCovariates
    }
  } else {
    technical_design <- NULL
  }
  
  if (! is.null(full_design)) {
    full_design <- as.formula(paste("~", full_design))
  } else {
    full_design <- ~1
  }
  
  cat("\nFormula:\n")
  print(full_design)
  
  # covariates
  covariates <- labels(terms(full_design))
  covariates <- colnames(annots)[colnames(annots) %in% covariates]
 # if (nrow(unique(annots[, ..covariates])) < 2) return(NULL)
  
  
  # format boolean variables and factors
  for (j in covariates) {
    value <- annots[[j]]
    
    if (is.character(value)) {
      new_value <- factor(value)
      set(annots, j = j, value = new_value)
      
    } else if (is.logical(value)) {
      new_value <- as.integer(value)
      set(annots, j = j, value = new_value)
      
    } else if (is.factor(value)) {
      new_lvls <- str_remove_all(levels(value), "\\(|\\)|-")
      levels(annots[[j]]) <- new_lvls
    }
    
  }
  
  cat("\nAnnotations:\n")
  print(unique(annots[, ! c("Sample", "Basename", "Directory")]))
  
  # filter counts
  counts <- data$counts[, annots$Sample]
  counts <- filterGenes(counts,annots$Group)
  ercc <- filterGenes(data$ercc[, annots$Sample], min_count = 5, min_prop = .8,annots$Group)
  
  wald_test <- as.logical(model$WaldTest)
  if (wald_test) {
    cat("Will perform Wald test on regression coefficients.\n")
  }
  
  # selgenes <- sample(names(tissues)[names(tissues) %in% rownames(counts)], 100)
  reg <- tissueSpecificRegression(counts = counts, # [selgenes, ], 
                                  ercc = ercc, 
                                  annots = annots, 
                                  tissues = tissues,
                                  mean_design = full_design,  
                                  normalization_groups = annots$Group, 
                                  wald_test = wald_test,
                                  max_iterations = 50,
                                  cores = opt$core_count)
  
  if (is.null(biological_design) | ! is.null(technical_design)) {
    batch_corrected_counts <- counts# [selgenes, ]
  } else {
    batch_corrected_counts <- correctCountsQuantileMatrix(counts = counts,# [selgenes, ], 
                                                          annots = annots, 
                                                          nf = t(reg$final_norm_factors), # rep_len(1, ncol(counts)), 
                                                          coefficients = reg$coefficients,
                                                          dispersions = reg$dispersions, 
                                                          full_design = full_design, 
                                                          biological_design = biological_design)
  }
  

  # output results to CSV
  if (! is.null(opt$output_dir)) {
    
    coefficients <- reg$coefficients / log(2)
    colnames(coefficients) <- paste0(colnames(coefficients), "_log2FoldChange")
    
    prop <- reg$proportions
    
    df <- data.table(GeneName = rownames(coefficients), 
                     GeneSymbol = data$geneid[rownames(coefficients)]$GeneSymbol,
                     Tissue = tissues[rownames(coefficients)],
                     Converged = reg$convergence >= 0,
                     ConvergenceCode = reg$convergence,
                     log2Disp=log2(reg$dispersions),
                     prop,
                     coefficients)
                     
        
    
    if (wald_test == TRUE) {
      se <- reg$coefficients_se
      colnames(se) <- paste0(colnames(se), "_lfcSE")
      
      pvalue <- reg$pvalue
      colnames(pvalue) <- paste0(colnames(pvalue), "_pvalue")
      
      padj <- reg$padj
      colnames(padj) <- paste0(colnames(padj), "_padj")
      
      df <- cbind(df, data.table(coefficients,
                                 se,
                                 pvalue, 
                                 padj))
    }
    
    cat("Saving table\n")
    fwrite(df, paste0(opt$output_dir, "/", model$Name,opt$output_suffix, ".csv.gz"))
  }
  
  # output results
  out <- c(reg, 
           list(batch_corrected_counts = batch_corrected_counts, 
                ercc = ercc, 
                annots = annots, 
                tissues = tissues))
})
names(results) <- models$Name

# save results
if (! is.null(opt$output_rds)) {
  cat("Saving RDS object\n")
  
  saveRDS(results, opt$output_rds)
}

cat("Done\n")

