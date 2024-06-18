
library(DESeq2)
library(dplyr)
library(scran)
library(data.table)
library(pbapply)
library(Hmisc)
library(optparse)
pbo <- pboptions(type="txt")
 
 
run_default_settings = F
if (run_default_settings){
	setwd(r"(Z:\nstroustrup\projects\asynch-seq-2022)")
}else setwd("/users/nstroustrup/nstroustrup/projects/asynch-seq-2022")

source("scripts/helpers/preprocess.R")
source("scripts/helpers/differential_analysis.R")
source("scripts/helpers/gene_variability.R")

opt <- list(model = "",
             nboots = 250,
             modelindex = 4,
             covariate_specific_tech_fits = F,
             reanalyze_boot = F
            )
            
            
option_list = list(
   make_option(c("-m", "--model"), 
               type="character",
               default="", 
              help="Name of model in differential analysis models"),
              
	 make_option(c("-n", "--nboots"), 
		     type="character",
		     default=NULL, 
      help="Number of bootstrap replicates."),
              
	 make_option(c("-i", "--modelindex"), 
		     type="character",
		     default=-1, 
      help="Index of model in differential analysis models."),
      
	make_option(c("-c", "--covariate_specific_tech_fits"), 
		     type="logical",
		     default=F, 
	      help="Should residual overdispersion be calculated relative to separate fits for each covariate?"),
	
	make_option(c("-r", "--reanalyze_boot"), 
		     type="logical",
		     default=F, 
	      help="Re-analyze stored boot reps")

      )

if (!run_default_settings){
	opt_parser <- OptionParser(option_list = option_list, prog = "Fit negative binomial model")
	opt <- parse_args(opt_parser)
}
opt$modelindex = as.integer(opt$modelindex);
if (opt$model == "" && opt$modelindex == -1)
	stop("Either a model or a model index must be specified.");
if (opt$model != "" && opt$modelindex != -1)
	stop("Either a model or a model index must be specified.");
calculate_separate_tech = opt$covariate_specific_tech_fits

nboots = opt$nboots;

if (1){
	
	models <- read.csv("./data/differential_analysis/models/single_and_pooled_worms.csv")
	if (opt$modelindex != -1){
			if (opt$modelindex > dim(models)[1])
				stop("Invalid model index")
				
		if (opt$reanalyze_boot){
			transl = c("n2_sw","bd_sw");
		}else{
			transl = c("n2_sw","bd_sw");
			#if (calculate_separate_tech){
			#	transl = 
			#}else transl =   
		}
		if (opt$modelindex > length(transl)){
			stop(paste0("Excess model index ",opt$modelindex))
		}
		#print(paste("Translating index",opt$modelindex,"to",transl[opt$modelindex]));
			
		cur_group = transl[opt$modelindex];
	}else{
		cur_group = opt$model;
	}
}

load_biocounts = function(){
	bc <- readRDS("./data/batch_corrected_counts/single_and_pooled_worms.rds")
	if (!(cur_group %in% names(bc$batch_counts_list)))
		stop(paste("The requested group:",cur_group,"is not present in the batch corrected counts list."))
	bio_counts <- bc$batch_counts_list[[cur_group]];
	
	rm(bc);
	gc();
	if (models$TissueSpecificNormalization[models$Name==cur_group]){
	    uniq_tissues <- sort(names(bio_counts))
	    ts <- lapply(uniq_tissues, function(tissue) {
		return(bio_counts[[tissue]])
		})
	   bio_counts = Reduce(rbind,ts);
	}
	return(bio_counts);
}

#groups_to_process = models$Name[models$EstimateGeneVariation == TRUE]

print(paste0("Loading ",cur_group))
data <- readRDS("./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds")
bio_counts = load_biocounts();

if (is.null(bio_counts))
	stop("No counts");

cur_bio_covariate = models$BiologicalCovariates[models$Name==cur_group]
batch_correct = models$BatchCorrect[models$Name==cur_group]

annots <- data$annots[Sample %in% colnames(bio_counts)]
setkey(annots, Sample)
annots[, Strain := droplevels(Strain)]
annots[, Genotype := droplevels(Genotype)]
annots[, Day := droplevels(Day)]
annots[, Group := droplevels(Group)]
annots[, Temperature := droplevels(Temperature)]   
annots[, Food := droplevels(Food)]  

tech_counts <- data$ercc[, annots$Sample]
tech_counts <- filterGenes(tech_counts, groups = NULL, min_count = 5, min_prop = .5)
tech_counts <- tech_counts[! rownames(tech_counts) %in% c("ERCC-00113", "ERCC-00136"), ]

rm(data)
gc()

sample_list <- lapply(levels(annots$Group), function(group) annots[annots$Group==group & annots$Exclude == F,]$Sample)
names(sample_list) <- levels(annots$Group)
groups <- annots$Group
grouped_annots <- annots[, .N, by = .(Group, Strain, Day, Temperature, Food,Pooled, BiologicalReplicate)]


bio_counts <- filterGenes(bio_counts[, annots$Sample], groups = annots$Group, min_count = 1, min_prop = .8)
#gene_list = c("WBGene00009306","WBGene00011158","WBGene00020088","WBGene00022042","WBGene00022501","WBGene00000149","WBGene00000254","WBGene00000293","WBGene00001209")
#gene_list = c("WBGene00022048","WBGene00009583")
  
  
suffix = ifelse(calculate_separate_tech,"","_shared_tech");
if (!opt$reanalyze_boot){
	if (is.null(bio_counts))
		stop("No Counts Loaded.");
	gene_N = dim(bio_counts)[1]
	split_groups = factor(floor(1:gene_N/5000))

	gene_groups = split(1:gene_N, split_groups)
	print(paste0("Processing ",cur_group, " with ", nboots, " bootstrap replicates; ", ifelse(calculate_separate_tech,"covariate-specific spike-ins","common spike-ins")))
	# iterate through gene sets 5000 genes at a time, to keep memory usage down
	resdisp_parts = pblapply(1:length(gene_groups),function(split_group){

		# bootstraping estimates of NB and residual over-dispersion
		resdisp <- lapply(seq_len(nboots), function(boot_index) {
		  print(paste("Boot",boot_index))
		  flush.console()
		  # select a random set of samples for each group (listed separately in sample_list)
		  boot_sample_list <- lapply(sample_list, function(samples) 
		    sample(samples, length(samples), replace = TRUE))

		  #boot_sample_list <- lapply(sample_list, function(samples) samples)

		  #merge together all samples from all groups together
		  boot_samples <- unlist(unlist(boot_sample_list), use.names = FALSE)
		  boot_annots <- copy(annots[boot_samples, on = "Sample"])
		  boot_bio_counts <- bio_counts[, boot_samples]
		  boot_tech_counts <- tech_counts[, boot_samples]
		  rm(boot_samples)
		  gc()

		  # rename samples so they are unique
		  boot_annots[, Sample := paste0("BootSample", stringr::str_pad(1:.N, width = 5, side = "left", pad = "0"))]
		  colnames(boot_bio_counts) <- boot_annots$Sample
		  colnames(boot_tech_counts) <- boot_annots$Sample
		  boot_groups <- boot_annots$Group
		  #separate out samples in each group.
		  boot_sample_list <- lapply(levels(boot_groups), function(group) 
		    boot_annots[group, on = "Group"]$Sample)
		  names(boot_sample_list) <- levels(boot_groups)
		  
		  # estimate negative binomial and residual over-dispersion
		  out <- estimateResDisp(bio_counts = boot_bio_counts, 
					 tech_counts = boot_tech_counts, 
					 annots = boot_annots, 
					 sample_list = boot_sample_list, 
					 groups = boot_groups,
					 gene_list=gene_groups[[split_group]],
					 calculate_separate_tech=calculate_separate_tech)

		  if (! is.null(out$params)) {
		    out$params[, BootIndex := boot_index]
		  }

		  if (! is.null(out$trend)) {
		    out$trend[, BootIndex := boot_index]
		  }
		  out
		})

		gc()
		return(resdisp)
	})
	print("Finished Bootstrapping.")
	flush.console()

	#Data is formatted like this:
	#list( list( list(trend,params),list(trend_params),...,nboots),list(list(trend,params),list(trend,params),...,nboots),...,length(grene_groups) )
	trend_fun = function(x){rbindlist(lapply(x,`[[`, "trend"))};
	trend = rbindlist(lapply(resdisp_parts, trend_fun))

	print("Writing Trend Info.")
	# format trend data (the ercc parametric fits, mostly) results into data.frame
	#trend <- rbindlist(lapply(resdisp, `[[`, "trend"))
	trend <- merge(trend, grouped_annots, by = "Group")
	setkey(trend, Group)
	fwrite(trend, paste0("./data/gene_variability/trend_",cur_group,suffix,".csv.gz"))
}else{
	print(paste("Re-analyzing existing bootstrap replicates with",ifelse(calculate_separate_tech,"covariate-specific spike-ins","common spike-ins")))
	trend = fread(paste0("./data/gene_variability/trend_",cur_group,suffix,".csv.gz"))
}


# params <- rbindlist(lapply(resdisp, `[[`, "params"))

if (!opt$reanalyze_boot){
	param_fun = function(x)rbindlist(lapply(x,`[[`, "params"))
	params = rbindlist(lapply(resdisp_parts, param_fun))

	print("Merging estimate table.")
	params <- merge(params, grouped_annots, by = "Group")
	setkey(params, GeneName, BootIndex, Group)
	fwrite(params, paste0("./data/gene_variability/residual_overdispersion_",cur_group,suffix,".csv.gz"))
	rm(resdisp_parts, envir = .GlobalEnv)
	resdisp_parts <<- NULL
	gc()
}else{
	params = fread(paste0("./data/gene_variability/residual_overdispersion_",cur_group,suffix,".csv.gz"))
	params[, Day := factor(Day, levels(annots$Day))]
	params[, Strain := factor(Strain, levels(annots$Strain))]
	params[, Food := factor(Food, levels(annots$Food))]
	params[, Temperature := factor(Temperature, levels(annots$Temperature))]
}

#set up factor levels for differential analysis II.
tmp = set_reference_groups_from_model(models[models$Name==cur_group,],
				      params[,.(Strain, Day, Food, Temperature,Pooled)]
				      );
params[,Strain:=tmp$Strain][,Day:=tmp$Day][,Food:=tmp$Food][,Temperature:=tmp$Temperature][,Pooled:=tmp$Pooled]
rm(tmp)

# summarize bootstrap results and compute HVG p-values

# summarize trend results
print("Summarizing Trend Info.")
trend_df <- summarizeTrend(trend)
fwrite(trend_df, paste0("./data/gene_variability/summary_trend_",cur_group,suffix,".csv.gz"))
rm(trend)
gc()

print("Summarizing estimates.")
params_df <- summarizeResDisp(params)
fwrite(params_df, paste0("./data/gene_variability/summary_residual_overdispersion_",cur_group,suffix,".csv.gz"))


print("Running Differential Analysis.")
form <- paste("~", cur_bio_covariate)
if (batch_correct)
form = paste(form,"+ BiologicalReplicate")
form = as.formula(form);
if (!calculate_separate_tech){
	#params include both separate and consensus spike-ins
	#only use the consensus ones
	params = params[SequenceType != "Technical",]
}
d8_dva_df <- runDifferentialVariabilityAnalysis(params, form)
fwrite(d8_dva_df, paste0("./data/gene_variability/regression_effects_",cur_group,suffix,".csv.gz"))

print("Done!")