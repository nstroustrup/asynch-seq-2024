# batch correction
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(DESeq2)
  library(stringr)
  library(optparse)
})

use_defaults = F
if (!use_defaults)
	setwd("/users/nstroustrup/nstroustrup/projects/asynch-seq-2022")
	

source("./scripts/helpers/preprocess.R")
source("./scripts/helpers/batch_correction.R")

opt <- list(name = "single_and_pooled_worms",
            deseq_dir = "./data/differential_analysis/deseq/single_and_pooled_worms",
            models = "./data/differential_analysis/models/single_and_pooled_worms.csv",
            counts = "./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds",
            output = "./data/batch_corrected_counts/single_and_pooled_worms.rds",
            tissue_normalization = T
            )

# handle arguments
option_list = list(
  
  make_option(c("-n", "--name"), 
              type="character",
              default=NULL, 
              help="Name of experiment."),
  
  make_option(c("-i", "--deseq_dir"), 
              type="character",
              default=NULL, 
              help="Input directory containing model DESeq2 RDS objects."),
  
make_option(c("-S", "--input_suffix"),
              type="character",
              default=NULL,
              help="suffix for input filenames."),

  make_option(c("-c", "--counts"), 
              type="character",
              default=NULL, 
              help="Formatted counts RDS file."),
  
  make_option(c("-m", "--models"), 
              type="character",
              default=NULL, 
              help="Input models CSV file."),
  
  make_option(c("--tissue_normalization"), 
              default=FALSE,
              action="store_true",
              help="Should tissue-specific normalization be performed?"),
  
  make_option(c("-o", "--output"), 
              type="character",
              default=NULL,
              help="Output RDS file."))
if (!use_defaults){
	opt_parser <- OptionParser(option_list = option_list, prog = "Batch Correction")
	opt <- parse_args(opt_parser)
}
if (is.null(opt$name)) {
  stop("Missing name of experiment (--name)")
}

if (is.null(opt$models)) {
  stop("Missing models CSV file (--models)")
}
if (is.null(opt$counts)) {
  stop("Missing counts file (--counts)")
}

if (is.null(opt$output)) {
  stop("Missing output RDS file (--output)")
}

# load data
cat("Loading data\n")
models <- fread(opt$models)
if (1){
	
	raw_counts <- readRDS(opt$count)
	
	fn = sapply(models$Name,function(model_name){
		fn = paste0(opt$deseq_dir,"/",model_name,opt$input_suffix,".rds");
		if(!file.exists(fn)){
			print(paste("Cannot find deseq file",fn))
			return(F);
		}
		return(T)
	})
	if (any(fn)) warning("Some deseq files are missing")
	#names(dds_list) = models$Name
	#models <- models[names(dds_list), on = "Name"]
	#models = models[!(TissueSpecificNormalization == TRUE & BatchCorrect==FALSE),]

	#if (any(models$TissueSpecificNormalization==T)) {
	  tissues_df <- fread("./data/tissue_unique_genes/genes_unique_to_tissues.csv.gz")
	  tissues <- pull(tissues_df, "Tissue", "GeneName")
	#}
}
if (1){
# batch corrections
#Important explanation
#When we run correctCountsQuantile(), the correction is performed.
#The full regression design is taken from the dds object from the --input file, but
#this dds object in term takes the full design from the models file (usually the same as specified in --models commandline)
#The biological design is taken from the models file specified in the --models commandline
#The batch correction is performed on the *difference* between these two designs
#e.g if the full design is ~ Strain + BiologicalReplicate
#and the biological design is ~Strain
#Then the regression will normalize for systematic differences between BiologicalReplicate
#across samples.

cat("Started performing batch corrections\n")
batch_counts_list <- lapply(models$Name, function(name) {
  cat(paste0(name,"\n"));
  model <- models[name, on = "Name"]

  ##dds <- dds_list[[name]]

  
  fn = paste0(opt$deseq_dir,"/",name,opt$input_suffix,".rds")
  if (!file.exists(fn))
     return(NULL)
  dds = readRDS(fn)
  if (is.null(dds)){
	  cat(paste0("No differential analysis found for ",name,"\n"));
	  return(NULL)
  }
  if (model$TissueSpecificNormalization==T){
  	cat(paste0("No batch correction requested for TS analysis ",name,"\n"));
  	if (model$BatchCorrectionType != "standard")
  		stop("Exotic batch correction types are not implemented for tissue-specific normalization")
  	 uniq_tissues <- sort(names(dds))
  	 if ("counts" %in% uniq_tissues)
  	 	stop("Tissue-specific analysis has been specified but the deseq model file contained non-tissue-specific analysis for this model.")
	  ts <- lapply(uniq_tissues, function(tissue) {
	    #see important explanation above
	    cat(paste0(tissue,"..."))
	    out <- correctCountsQuantile(dds[[tissue]], biological_design = as.formula(paste("~", model$BiologicalCovariates))) 
	  })
	  cat("\n");
  	names(ts) = uniq_tissues
  	return(ts)
  }
 
  if (model$BatchCorrect == FALSE){
  	cat(paste0("No Batch Correction Requested for ",name,"\n"));
  	#see important explanation above
  	
  	#when there is only one level of a covariate, differential analysis just returns the raw data.
  	if (class(dds) %in% c("data.frame","data.table")) 
  		return (dds);
 
  	return(DESeq2::counts(dds))
  }
  if (is.null(model$TechnicalCovariates) || model$TechnicalCovariates %in% c("","NULL")){
 	 cat(paste0("Batch Correction requested with no technical covariates for ",name,"\n"));
  	return(NULL)
  }
  #now we do the batch correction
  if (model$BatchCorrectionType == "standard"){
	  cat("\n")
	  cat(name)
	  cat(paste("\nTechnical Covariate:", model$TechnicalCovariates, "\n\n"))
	  #see important explanation above
	  out <- correctCountsQuantile(dds, biological_design = as.formula(paste("~", model$BiologicalCovariates)))
	  return(out)
  }else if (model$BatchCorrectionType == "stratify_by_bio_covariate"){
  	#here we batch correct *within* each biological covariate
  	#but do not assume batch correction factors should be consistant across biological covariates
  	#this should be used only in very specific cases, e.g comparing across different sequencing methods
  	#collected across replicates as is done in the SmartSeq2 vs. MARS-seq comparison
  	
  	if (class(dds) != "list" || ! "stratified" %in% names(dds))
  		stop("Improperly formatted dds");
  		
	#go through each combination of covariates and for each, run a normalization
  	out_c = lapply(dds[["stratified"]],function(dds_x){
  		
  		#only to fix a problem that's already solved while the fix is running
  		#if (class(dds_x) %in% c("matrix")){
  		#browser()
		#	dds_x = list(counts = dds_x,nf = rep(1,ncol(dds_x)))
		#}
		
		#when there is only one level of a covariate, differential analysis just returns the raw data.
		if (class(dds_x) == "list"){
			#browser()
			return(data.frame(t(sweep(dds_x$counts, 2, dds_x$nf, "/"))))
		}
			
		out <- correctCountsQuantile(dds_x, biological_design = ~1)
	 	return(data.frame(t(out)))
	})
	
	out = t(as.matrix(rbindlist(out_c,fill=T)))
	out = out[apply(out,1,FUN=function(x)!any(is.na(x))),]
	
	names_c = lapply(out_c,function(x){rownames(x)})
	colnames(out) = do.call("c",names_c)
  	
  	#put columns in original order and return
  	return(out[,colnames(DESeq2::counts(dds[["unstratified"]])) ]);
  	
  }else{
  	stop(paste("Unknown batch correction method specified:",model$BatchCorrectionType))
  }
})
# names(batch_counts_list) <- paste0("counts_", models$Name)
names(batch_counts_list) <- models$Name
}

cat("Performing normalizations\n")
nf_list <- lapply(models$Name[models$TissueSpecificNormalization == F], function(name) {
	if (is.null(batch_counts_list[[name]]))
		return(NULL);
 	model <- models[name, on = "Name"]
	annots <- copy(raw_counts$annots[eval(parse(text = model$Subset))])
  	if (model$use_optional_strain_prefix){
  		print("Using extra Strain Suffixes")  #allows us to analyze the same data in different ways in different models
	  	annots$Strain = factor(str_trim(apply(annots,1,function(x)paste(as.character(x[["Strain"]]),x[["optional_strain_prefix"]]))))
  		annots$Group = factor(str_trim(apply(annots,1,function(x)paste0(as.character(x[["Group"]]),".",x[["optional_strain_prefix"]]))))
 	}
	# drop levels
	for (j in colnames(annots)) {
		value <- annots[[j]]
		if (is.factor(value)) {
			set(annots, j = j, value = droplevels(value))
		}
  	}
  	cat(paste0(name," "));
  
  	normalizationFactor(batch_counts_list[[name]], 
                      groups = annots$Group)
})  
if (any(models$TissueSpecificNormalization == F))
names(nf_list) <- models$Name[models$TissueSpecificNormalization == F]

cat("Performing tissue specific normalizations\n")
tissue_nf_list <- lapply(models$Name, function(name) {
    if (is.null(batch_counts_list[[name]]))
		return(NULL);
	cat(paste0(name," "));
   model <- models[name, on = "Name"]
   annots <- raw_counts$annots
   if (model$use_optional_strain_prefix){
	print("Using extra Strain Suffixes (2)")  #allows us to analyze the same data in different ways in different models
	annots$Strain = factor(str_trim(apply(annots,1,function(x)paste(as.character(x[["Strain"]]),x[["optional_strain_prefix"]]))))
	annots$Group = factor(str_trim(apply(annots,1,function(x)paste0(as.character(x[["Group"]]),".",x[["optional_strain_prefix"]]))))
   }
 	
   setkey(annots,Sample)
   #we need to get the group of each sample, so first we get the sample names
   #and then we look up the group using those sample names.
   #tissue specific regressions have multiple tables so we take the first tissue.
   if (class(batch_counts_list[[name]]) == "list"){
   #browser()	
   	groups = droplevels(annots[colnames(batch_counts_list[[name]][[1]]),]$Group)
   }else groups = droplevels(annots[colnames(batch_counts_list[[name]]),]$Group)
   
   if (model$TissueSpecificNormalization ==T){
    return(normalizationFactorTissueSpecific_list(batch_counts_list[[name]], 
    				       groups = groups))
   }else if (opt$tissue_normalization == T){
 
   	return(normalizationFactorTissueSpecific(batch_counts_list[[name]], tissues,
                                      groups = groups))
   }else{
   	return(NULL)
   }
  })
  names(tissue_nf_list) <- models$Name;


cat("Save results\n")
saveRDS(list(batch_counts_list = batch_counts_list, 
             nf_list = nf_list, 
             tissue_nf_list = tissue_nf_list), 
        opt$output)

cat("Done\n")
