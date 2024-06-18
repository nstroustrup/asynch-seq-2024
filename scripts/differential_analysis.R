# differential analysis
cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(dplyr)
  library(DESeq2)
  library(data.table)
  library(scran)
  library(stringr)
  library(optparse)
})

use_defaults = F
if (!use_defaults)
	setwd("/users/nstroustrup/nstroustrup/projects/asynch-seq-2022")

source("./scripts/helpers/preprocess.R")
source("./scripts/helpers/differential_analysis.R")
 opt <- list(name = "single_and_pooled_worms",
             input = "./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds", 
             output_dir = "./data/differential_analysis/tables/single_and_pooled_worms/",
             models = "./data/differential_analysis/models/single_and_pooled_worms.csv",
             output_dds_dir = "./data/differential_analysis/deseq/single_and_pooled_worms",
             node_id = 2
           
             )
  opt = list(
	  name ="single_worms_bootstrap",
	  input ="./data/gene_network/temp/counts_and_annots_tst.rds",
	  models ="./data/differential_analysis/models/single_worms_bootstrap.csv",
	  output_dir ="./data/gene_network/temp/differential_analysis/tables/",
	  output_suffix ="_tst" ,
	  output_dds_dir = "./data/gene_network/temp/differential_analysis/deseq",
	  node_id=-1
	  )
 if (0){
 opt <- list(name = "net_validation",
              input = "./data/formated_counts/counts_and_annots_net_validation.rds", 
              output_dir = "./data/differential_analysis/tables/net_validation/",
              models = "./data/differential_analysis/models/net_validation.csv",
              output_dds_dir = "./data/differential_analysis/deseq/net_validation",
              add_specific_model = "\"N2_nlp28_Set3_Day8\""
             )
 }


# handle arguments
option_list = list(
  
  make_option(c("-n", "--name"), 
              type="character",
              default=NULL, 
              help="Name of experiment."),
  
  make_option(c("-N", "--node_id"), 
                type="integer",
                default=-1, 
              help="Which node of the parallel processing should be run?"),
           
  make_option(c("-i", "--input"), 
              type="character",
              default=NULL, 
              help="Input formated counts RDS file."),
  
  make_option(c("-m", "--models"), 
              type="character",
              default=NULL, 
              help="Input models CSV file."),
  
  make_option(c("-o", "--output_dir"), 
              type="character",
              default=NULL,
              help="Output directory for result tables."),
  
  make_option(c("-r", "--output_dds_dir"), 
              type="character",
              default=NULL,
              help="Output directory for DESeq2 RDS objects."),

  make_option(c("-S","--output_suffix"),
    type="character",
    default="",
  help="suffix to add to output files"),
        	         
  make_option(c("-s", "--add_specific_model"), 
              type="character",
              default=NULL,
              help="add a specific model to the existing file"))

if (!use_defaults){
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
add_specific_model =F
specific_model_name = "";
if (!is.null(opt$add_specific_model) && opt$add_specific_model != ""){
	add_specific_model = T
	specific_model_name = eval(parse(text=opt$add_specific_model))
}

# load data
cat("Loading Differential Analysis Model Specifications\n")

models <- fread(opt$models)

if (opt$node_id != -1){
	if (opt$node_id > nrow(models))
		stop(paste("Invalid node ID specified:",opt$node_id))
	add_specific_model = T
	specific_model_name = models$Name[opt$node_id]
	cat(paste("Node",opt$node_id,"processing",specific_model_name,"\n"))
}

models <- models[Run == TRUE]
if (any(duplicated(models$Name))) {
  print(table(models$Name))
  stop("Duplicated names in model file!")
}

if (any(as.logical(models$TissueSpecificNormalization))) {
  tissues_df <- fread("./data/tissue_unique_genes/genes_unique_to_tissues.csv.gz")
  tissues <- pull(tissues_df, "Tissue", "GeneName")
  uniq_tissues <- sort(unique(tissues))
}

# create directories
if (! is.null(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (! is.null(opt$output_dds_dir)) {
  dir.create(opt$output_dds_dir, recursive = TRUE, showWarnings = FALSE)
}else stop("No output file specified");


run_differential_analysis_internal = function(model,annots,only_validate_metadata){

  if (model$use_optional_strain_prefix==T){
  	print("Using extra Strain Suffixes")  #allows us to analyze the same data in different ways in different models
  	annots$Strain = str_trim(apply(annots,1,function(x)paste(x[["Strain"]],x[["optional_strain_prefix"]])))
  	annots$Group = str_trim(apply(annots,1,function(x)paste0(x[["Group"]],".",x[["optional_strain_prefix"]])))
  }
  # drop levels
  for (j in colnames(annots)) {
    value <- annots[[j]]
    if (is.factor(value)) {
      set(annots, j = j, value = droplevels(value))
    }
  }

  if (nrow(annots) == 0){
	stop("Data selection returned no annotations");
	return(NULL)
  }

  # handle specific case
  if (opt$name == "tissue_specific_rnai" | opt$name == "net_validation") {
    annots[, RNAi := relevel(RNAi, "EV")]
  }
  if (opt$name == "glp1_glp4") {
    levels(annots$Genotype) <- c("EV", "glp1", "glp4")
  }

  annots[, Pooled := ifelse(NumberOfWorms != 1,"pooled","single")]


  # make formula
  if (!(is.null(model$BiologicalCovariates) || model$BiologicalCovariates %in% c("","NULL"))){
    form <- model$BiologicalCovariates
    annots = set_reference_groups_from_model(model,annots);
  } else {
    form <- NULL
  }
  if (model$TechnicalCovariates != "NULL") {
    if (! is.null(form)) {
      form <- paste(form, "+", model$TechnicalCovariates)
    } else {
      form <- model$TechnicalCovariates
    }
  }

  form <- as.formula(paste("~", form))
	
  cat("\nFormula:\n")
  print(form)

  # covariates
  covariates <- labels(terms(form))
  covariates <- colnames(annots)[colnames(annots) %in% covariates]



  # format boolean variables and factors
  for (j in covariates) {
    value <- annots[[j]]

    if (is.character(value)) {
      new_value <- factor(value)	#this code should run only for technical factors,
					#because biological covariates are handled by set_reference_groups_from_model() above.
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
if (any(!is.na(annots$optional_strain_prefix))){  
tab = aggregate(as.formula(paste0("Sample ~ ",paste(names(annots)[! names(annots) %in% c("Sample", "Basename", "Directory")],collapse="+"))),annots,FUN=length);
}else{
tab = aggregate(as.formula(paste0("Sample ~ ",paste(names(annots)[! names(annots) %in% c("Sample", "Basename", "Directory","optional_strain_prefix")],collapse="+"))),annots,FUN=length);
}
  names(tab)[names(tab)=="Sample"] = "Sample Count"
  #print(unique(annots[, ! c("Sample", "Basename", "Directory")]))
  print(tab)
  if(only_validate_metadata) return(NULL)

  # filter counts
  counts <- data$counts[, annots$Sample]
  if (!any(annots$LenientCountCriteria)){
  	print("Removing low quality samples")
  	counts <- filterGenes(counts, annots$Group)
  }
  
  # normalize
  if (as.logical(model$TissueSpecificNormalization)) {
    cat("Running tissue-specific normalization\n")
    nf <- normalizationFactorTissueSpecific(counts, tissues, annots$Group)
  } else {
    cat("Running standard normalization\n")
    nf <- normalizationFactor(counts, annots$Group)
  }
  if (nrow(unique(annots[, ..covariates])) < 2){
    	cat("Less than two levels were found among samples on which to run differential analysis\n")
    	return(list(counts=counts,nf=nf));
  }

  # run diff. ana.
  cat("Running differential analysis\n")
  if (as.logical(model$TissueSpecificNormalization)) {
    dds <- lapply(uniq_tissues, function(tissue) {
      genes <- names(tissues)[tissues == tissue]
      tissue_counts <- counts[rownames(counts) %in% genes, ]
      tryCatch((runDESeq(counts = tissue_counts, annots = annots, nf = nf[[tissue]], form)),
	       error = function(e) {
		 print(paste("Error running Tissue-specific DESeq::",e))
		 NULL
	       })
    })
    names(dds) <- uniq_tissues

  } else {
    dds <- tryCatch((runDESeq(counts = counts, annots = annots, nf = nf, form)),
		    error = function(e) {
		      print(paste("Error running DESeq::",e))
		      NULL
		    })

  }

  if (! is.null(opt$output_dir)) {
    cat("Saving table\n")
    if (as.logical(model$TissueSpecificNormalization)) {
      da <- rbindlist(lapply(uniq_tissues, function(tissue) {
	out <- formatDESeq(dds[[tissue]], geneid = data$geneid)
	out[, Tissue := tissue]
	out
      }))
    } else {
      da <- formatDESeq(dds, data$geneid)
    }
    fwrite(da, paste0(opt$output_dir, "/", model$Name, opt$output_suffix,".csv.gz"))
  }

  cat("-------\n\n")

  dds
}

# run differential analysis
cat("Started differential analyses\n")
run_differential_analysis = function(models,only_validate_metadata=F){
	lapply(seq_len(nrow(models)), function(i) {
	 	model <- models[i, ]
		  cat("-------\n")
		  if (!only_validate_metadata){
		  	cat(paste("#### Running", i, "out of", nrow(models)))
			  cat("\nModel\n")
			  print(model)
		  }else cat(paste("#### Validating metadata for ",model$Name))
	  	cat("\n")
	  	# get annotations
		annots <- copy(data$annots[eval(parse(text = model$Subset))])

	  	if (model$BiologicalCovariates == "NULL" || model$BatchCorrectionType == "standard"){
	  		res = run_differential_analysis_internal(model,annots,only_validate_metadata);
	  	}else if (model$BatchCorrectionType == "stratify_by_bio_covariate"){
	  		res = list();
	  		cat("Stratified batch correction; running biological covariate model\n")
	  		#first, we compare biological covariates without accounting for technical covariates
	  		#This is the differential analysis used to understand the effect of the biological covariate
	  		new_model = model;
	  		new_model$TechnicalCovariates = "NULL";
	  		res[["unstratified"]] = run_differential_analysis_internal(new_model,annots,only_validate_metadata);
	  		
	  		#second, we stratify the data according to the biological covariates and compare across technical covariates
			#these differential analyses are required to do the batch corrections.
	  		covariate_list = aggregate(as.formula(paste("Sample ~", model$BiologicalCovariates)),annots,FUN=length)
		  	covariate_list = covariate_list[,names(covariate_list)!="Sample", drop = FALSE];
		  	sub_model = model;
		  	sub_model$BiologicalCovariates = "NULL";
		  	#go through each combination and for each, run a normalization
		  	#concatenate them all together
		  	res[["stratified"]] = lapply(1:nrow(covariate_list),function(i){
	  			cat(paste0("Stratified batch correction; running technical covariate model for ",paste(covariate_list[i,],collapse=" "),"\n"))
		  		sub = annots[,get(names(covariate_list))] == covariate_list[i,]
		  		sub_annots = annots[sub,]
				#tech_covariate_list = aggregate(as.formula(paste("Sample ~", model$TechnicalCovariates)),annots,FUN=length)
		  		#tech_covariate_list = covariate_list[,names(covariate_list)!="Sample", drop = FALSE];
		  		return(run_differential_analysis_internal(sub_model,sub_annots,only_validate_metadata));
		  	})
		  
	  	} else stop(paste("Unknown Batch correction type:",model$BatchCorrectionType));
	  	if (!only_validate_metadata)
	 		saveRDS(res, paste0(opt$output_dds_dir,"/",model$Name,opt$output_suffix,".rds"))
	})
}

#now we load all the data we need and run the differential analysis
cat("Loading Sequence Data")
data <- readRDS(opt$input)

if (add_specific_model){
	models = models[Name%in%specific_model_name,]
	specific_model_name = models$Name;
}

#first we validate, so we catch any errors in the models file quickly
run_differential_analysis(models,T)


run_differential_analysis(models,F)

cat("Done\n")