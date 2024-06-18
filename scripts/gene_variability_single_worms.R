
library(DESeq2)
library(dplyr)
library(scran)
library(splines)
library(data.table)
library(pbapply)
library(Hmisc)
library(optparse)
library(stringr)
library(qs)
pbo <- pboptions(type="txt")
 
 
run_default_settings = F
if (0 & run_default_settings){
	setwd(r"(Z:\nstroustrup\projects\asynch-seq-2022)")
}else setwd("/users/nstroustrup/nstroustrup/projects/asynch-seq-2022")

source("scripts/helpers/preprocess.R")
source("scripts/helpers/differential_analysis.R")
source("scripts/helpers/gene_variability.R")

if (1) opt <- list(model = "ts_d8_UV",
             nboots = 4,
             modelindex = 2,
             covariate_specific_tech_fits = F,
             reanalyze_boot = F,
             split_bootstraps_across_nodes = 2,
             differential_analysis_model=NULL,
             differential_analysis_model_name="",
             differential_analysis_model_foldchange_column="",
             differential_analysis_comparing_da_model=F,
             differential_analysis_RNAi_GeneName_to_remove=NULL,
             differential_analysis_RNAi_GeneSymbol_to_remove=""
            )
if (0) opt <- list(model="n2_sw",
             nboots = 250,
             modelindex = -1,
             covariate_specific_tech_fits = T,
             reanalyze_boot = T,
             split_bootstraps_across_nodes = 50,
            )            
            
if (0) opt <- list(model="n2_sw",
     nboots = 250,
     modelindex = 1,
     covariate_specific_tech_fits = F,
     reanalyze_boot = F,
     split_bootstraps_across_nodes = 100,
     differential_analysis_model = "data/differential_analysis/tables/net_validation/N2_nlp28_Set3_Day8.csv.gz",
     differential_analysis_model_name = "da_nlp28",
     differential_analysis_model_foldchange_column = "RNAi_nlp28_vs_EV_log2FoldChange",
     differential_analysis_comparing_da_model = T,
     differential_analysis_RNAi_GeneName_to_remove = NULL,
     differential_analysis_RNAi_GeneSymbol_to_remove = "nlp-28"
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
	      help="Re-analyze stored boot reps"),
	      
	make_option(c("-s", "--split_bootstraps_across_nodes"), 
			     type="numeric",
			     default=1, 
	      help="How many nodes should each bootstrap be distributed across"),
	      
	make_option(c("-d", "--differential_analysis_model"), 
               type="character",
			     default=NULL, 
      		help="Specify a differential analysis file to model, and then estimate variabilty of the residuals"),
      		
	make_option(c("-e", "--differential_analysis_model_name"), 
               type="character",
				     default=NULL,
		help="Specify the nickname of this reference, to be used in filenames"),
		
	make_option(c("-f", "--differential_analysis_model_foldchange_column"), 
               type="character",
				     default=NULL, 
		help="Specify the name of the column in the differential analysis file to use for log2fold change data"),
		
	make_option(c("-D", "--differential_analysis_comparing_da_model"), 
		       type="logical",
					     default=F, 
	help="Instead of comparing model covariates, consider only the reference condition and compare counts to da model results"),
	
	make_option(c("-y", "--differential_analysis_RNAi_GeneSymbol_to_remove"), 
			       type="character",
						     default=NULL, 
	help="Instead of comparing model covariates, consider only the reference condition and compare counts to da model results"),
	
	make_option(c("-g", "--differential_analysis_RNAi_GeneName_to_remove"), 
			       type="character",
						     default=NULL, 
	help="Instead of comparing model covariates, consider only the reference condition and compare counts to da model results"),
	
	make_option(c("-S", "--specialized_subset"), 
				       type="logical",
							     default=F, 
	help="Only analyze a programatically-defined subset"),
	make_option(c("-p", "--print_count"), 
				       type="logical",
							     default=F, 
	help="Print the number of nodes required to complete the distributed task")
      )     
     
if (!run_default_settings){
	opt_parser <- OptionParser(option_list = option_list, prog = "Fit negative binomial model")
	opt <- parse_args(opt_parser)
}
opt$split_bootstraps_across_nodes = as.integer(opt$split_bootstraps_across_nodes)
if (is.na(opt$split_bootstraps_across_nodes))
	stop("split_bootstraps_across_nodes must be a number")
opt$modelindex = as.integer(opt$modelindex);
if (opt$model == "" && opt$modelindex == -1)
	stop("Either a model or a model index must be specified.");
calculate_separate_tech = opt$covariate_specific_tech_fits

nboots = as.integer(opt$nboots);

if (1){
	
	models <- read.csv("./data/differential_analysis/models/single_and_pooled_worms.csv")
	preset_model = (opt$model != "")
	if (preset_model)
		cur_group = opt$model;
	if (opt$specialized_subset != ""){
		sub = c("d8_UV","d1_UV","d8_glp1","d1_glp1","d8_starvation","d3_starvation","d1_sw","d8_sw","d8_20v25C","d1_20v25C")
		models = models[models$Name %in% sub,]
		print(paste("Only analyzing subset of models:",paste(sub,collapse=",")))
	}

	if (opt$split_bootstraps_across_nodes == 1 || opt$reanalyze_boot){
		if (opt$model != "" & opt$modelindex != -1)
			stop("With all bootstrap replicates on a single node, you can't specify both the model and the model index")

		split_bootstrap_node_index = 0
		boot_file_suffix = ""
		boot_file_prefix = "temp/";
		bootstrap_replicates_to_run = 1:nboots;
		model_index = opt$modelindex 
	}
	else{
		if (opt$modelindex == -1)
			stop("Bootstrap split has been requested across multiple nodes--a model index must be specified for each node.")
		print(paste("Bootstrap replicates will be split across",opt$split_bootstraps_across_nodes, "nodes"))
		bootstrap_replicates_per_node = ceil(nboots/opt$split_bootstraps_across_nodes)
		print(paste("Running ",opt$split_bootstraps_across_nodes, "bootstrap replicates per processing node"))
		if (opt$model != ""){	#analyzing a specific model
			split_bootstrap_node_index = opt$modelindex;
		}else{ #analyzing all models
			model_index = floor((opt$modelindex-1)/opt$split_bootstraps_across_nodes)+1;
			split_bootstrap_node_index = (opt$modelindex-1)%%opt$split_bootstraps_across_nodes+1;
		}
		boot_file_suffix = paste0("_boot_",split_bootstrap_node_index);

		bootstrap_replicates_to_run = ((split_bootstrap_node_index-1)*bootstrap_replicates_per_node):(split_bootstrap_node_index*(bootstrap_replicates_per_node)-1)
		print(paste("This is node number ",split_bootstrap_node_index,"running bootstrap replicates",paste(range(bootstrap_replicates_to_run),collapse="-"),"for this model"))

		boot_file_prefix = "temp/";
		dir.create("./data/gene_variability/temp", showWarnings = FALSE)
	}
	if (opt$print_count){
		print(paste("Max model index ID:",(dim(models)[1]-1)*opt$split_bootstraps_across_nodes+1))
	}
	if (!preset_model){
		if (opt$split_bootstraps_across_nodes > 1)
			if (model_index> dim(models)[1])
				stop("Invalid model index")
		cur_group = models$Name[model_index];
	}
}
da_model = NULL;
if (!is.null(opt$differential_analysis_model)){
	print(paste("Using differential analysis model",opt$differential_analysis_model));
	if(is.null(opt$differential_analysis_model_foldchange_column))
		stop("To use a differential analysis model, you must specify the column from which to get logfoldchanges.")
		
	if(is.null(opt$differential_analysis_model_name))	
		stop("To use a differential analysis model, you must specify the nickname for this DA model, to be used in filenames.")
	
	if(is.null(opt$differential_analysis_RNAi_GeneName_to_remove) && is.null(opt$differential_analysis_RNAi_GeneSymbol_to_remove)){
		str = "You have not specified a transcript to remove when doing DA normalization.  This is usually important when using DA derived from RNAi"
		cat(paste(str,"\n"))
		warning(str)
		RNAi_GeneName_to_exclude = "";
	}else{
		if (!is.null(opt$differential_analysis_RNAi_GeneName_to_remove)){
			RNAi_GeneName_to_exclude = opt$differential_analysis_RNAi_GeneName_to_remove
		}else{
			symb = opt$differential_analysis_RNAi_GeneSymbol_to_remove;
			
			geneid <- fread("./data/annotations/c_elegans.PRJNA13758.WS265.geneIDs.txt.gz")
			geneid[V3 == "", V3 := V4]
			geneid <- geneid[, .(GeneName = V2, GeneSymbol = V3, GeneType  = V6)]
			geneid <- geneid[! duplicated(GeneName)]
			setkey(geneid, GeneSymbol)
			
			if (!symb %in% geneid$GeneSymbol)
				stop("Could not disambiguate specified differential_analysis_RNAi_GeneName_to_remove");
			RNAi_GeneName_to_exclude = geneid[symb,GeneName];
		}
	}	
		
	da_model = fread(opt$differential_analysis_model);
	if (!("GeneName" %in% names(da_model))) stop("The differential analysis file did not have a GeneName column")
	if (!(opt$differential_analysis_model_foldchange_column %in% names(da_model))) stop("The differential analysis file did not have the specified logfoldchange column")
	da_model = data.table(GeneName = da_model$GeneName,
			      FC = da_model[,get(opt$differential_analysis_model_foldchange_column)]
			      );
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
if(1){
data <- readRDS("./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds")
bio_counts = load_biocounts();
}




if (is.null(bio_counts))
	stop("No counts");

cur_bio_covariate = models$BiologicalCovariates[models$Name==cur_group]
batch_correct = models$BatchCorrect[models$Name==cur_group]

annots <- data$annots[Sample %in% colnames(bio_counts)]


tech_counts <- data$ercc[, annots$Sample]
tech_counts <- filterGenes(tech_counts, groups = NULL, min_count = 5, min_prop = .5)
tech_counts <- tech_counts[! rownames(tech_counts) %in% c("ERCC-00113", "ERCC-00136"), ]
colnames(tech_counts) = annots$Sample

rm(data)
gc()

if (opt$differential_analysis_comparing_da_model){
	#if we are comparing the data +- the da model, load only the control group
	cols = get_reference_groups_and_columns_from_model(models[models$Name==cur_group,],annots);
	if(length(cols) == 1){	annots = annots[get(cols[[1]]["column"]) == cols[[1]]["group"],]
	}else{			annots = annots[get(cols[[1]]["column"]) == cols[[1]]["group"] & get(cols[[2]]["column"]) == cols[[2]]["group"],]}
	
	#grab just the samples from the control group
	bio_counts = bio_counts[,annots$Sample]
	
	tech_counts = tech_counts[,annots$Sample]

	#now duplicate all these samples, 1 for the counts and the other for the model residuals
	col_names = colnames(bio_counts);
	colnames(bio_counts) = NULL;
	bio_counts = cbind(bio_counts,bio_counts);
	colnames(bio_counts) = c(col_names,paste0(col_names,".da"))
	
	col_names = colnames(tech_counts);
	colnames(tech_counts) = NULL;
	tech_counts = cbind(tech_counts,tech_counts);
	colnames(tech_counts) = c(col_names,paste0(col_names,".da"))
	
	#set up the annotations so that the algorithm compares da to no da
	annots2 = annots;
	print("Discarding covariates and evaluating the effect of the differential analysis model.");
	annots$Group = apply(annots,1,function(x)paste0(x[["Group"]],".da_model"))
	annots$da_model_residuals = T
	annots$Sample = paste0(annots$Sample,".da")
	
	annots2$Group = apply(annots2,1,function(x)paste0(x[["Group"]],".no_da_model"))
	annots2$da_model_residuals = F
	ref_group = annots2$Group[1];
	
	annots = rbind(annots,annots2)
	annots$Group = factor(annots$Group);
	annots$Group = relevel(annots$Group,ref=annots2$Group[1])
}else annots$da_model_residuals = T	#If we are running a model and there isn't something fancy going on (eg comparing model to non-model,
					#then we run the model across all samples.


if (!is.null(da_model)  && !opt$reanalyze_boot){

	if (RNAi_GeneName_to_exclude != ""){
		if (!RNAi_GeneName_to_exclude %in% rownames(bio_counts)){
			stop("The DA transcript specified to be excluded from the da model was not present in the counts file.")
		}else{
			cat(paste("Excluding ",RNAi_GeneName_to_exclude,"from DA analysis.\n"))
		}
	}

	genes_to_use = intersect(rownames(bio_counts),da_model$GeneName);
	genes_to_use = genes_to_use[genes_to_use!= RNAi_GeneName_to_exclude]
	if (length(genes_to_use) < .75*nrow(bio_counts))
		warning(paste("There is a relatively small overlap between the genes detected in the differential analysis file and the counts file: ",length(genes_to_use)/nrow(bio_counts)))
	bio_counts = bio_counts[genes_to_use,]

	print("Running da model")
	#here, if requested, we use a differential analysis model to "remove" the effect
	#of one set of correlations from the count matrix, and run the subsequent analysis on the residuals.

	#note that we only do this for samples (columns0 where da_model_resulals is marked as true
	#this lets us compare model residuals to counts (eg data "normalized" by the model vs the actual counts)

	res = ns_fit_da_linear_model(bio_counts[,annots$da_model_residuals],da_model,differential_analysis_column_name="FC",RNAi_GeneName_to_exclude=RNAi_GeneName_to_exclude,keep_x=F)

	if (res$relative_variance > 1.05){
		warning(paste("The relative variance after the da model was ",res$relative_variance))
	}else{
		if (nrow(bio_counts) != nrow(res$residuals))
			bio_counts = bio_counts[rownames(res$residuals),]
		nms_1 = rownames(bio_counts)
		bio_counts[,annots$da_model_residuals] = 2^(res$residuals + res$log2_intercepts)
		rownames(bio_counts) = nms_1
	}
}
if (opt$reanalyze_boot){
	print("Aggregating bootstrap results and running differential variation analysis.")
}


if (models$use_optional_strain_prefix[models$Name==cur_group]==T){
  	print("Using extra strain suffixes")  #allows us to analyze the same data in different ways in different models
	annots$Strain = factor(str_trim(apply(annots,1,function(x)paste(as.character(x[["Strain"]]),x[["optional_strain_prefix"]]))))
	annots$Group = factor(str_trim(apply(annots,1,function(x)paste0(as.character(x[["Group"]]),".",x[["optional_strain_prefix"]]))))
  }
  
setkey(annots, Sample)
annots[, Strain := droplevels(Strain)]
annots[, Genotype := droplevels(Genotype)]
annots[, Day := droplevels(Day)]
annots[, Group := droplevels(Group)]
annots[, Temperature := droplevels(Temperature)]   
annots[, Food := droplevels(Food)]  

sample_list <- lapply(levels(annots$Group), function(group) annots[annots$Group==group & annots$Exclude == F,]$Sample)
names(sample_list) <- levels(annots$Group)
groups <- annots$Group
grouped_annots <- annots[, .N, by = .(Group, Strain, Day, Temperature, Food,Pooled, BiologicalReplicate)]


bio_counts <- filterGenes(bio_counts[, annots$Sample], groups = annots$Group, min_count = 1, min_prop = .8)
#gene_list = c("WBGene00009306","WBGene00011158","WBGene00020088","WBGene00022042","WBGene00022501","WBGene00000149","WBGene00000254","WBGene00000293","WBGene00001209")
#gene_list = c("WBGene00022048","WBGene00009583")
  
  
suffix = ifelse(calculate_separate_tech,"","_shared_tech");
suffix = ifelse(!is.null(da_model),
		paste0("_",opt$differential_analysis_model_name,suffix),
		suffix)
suffix = ifelse(opt$differential_analysis_comparing_da_model,paste0("_","mod_comp",suffix),suffix)
print(paste0("Writing with suffix \"",suffix,"\""))

if (!opt$reanalyze_boot){
	if (is.null(bio_counts))
		stop("No Counts Loaded.");
	#we can estimate variance for all the genes, so just split them all into rounds of about 5000 genes
	gene_N = nrow(bio_counts)
	split_groups = factor(floor(1:gene_N/5000))
	gene_groups = split(1:gene_N, split_groups)

	print(paste0("Processing ",cur_group, " with ", nboots, " bootstrap replicates; ", ifelse(calculate_separate_tech,"covariate-specific spike-ins","common spike-ins")))
	# iterate through gene sets 5000 genes at a time, to keep memory usage down
	resdisp_parts = pblapply(1:length(gene_groups),function(split_group){

		# bootstraping estimates of NB and residual over-dispersion
		resdisp <- lapply(bootstrap_replicates_to_run, function(boot_index) {
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
		  out <- estimateResVar(bio_counts = boot_bio_counts, 
					 tech_counts = boot_tech_counts, 
					# annots = boot_annots, 
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
	fwrite(trend, paste0("./data/gene_variability/",boot_file_prefix,"trend_",cur_group,suffix,boot_file_suffix,".csv.gz"))
}else{
	if (opt$split_bootstraps_across_nodes == 1){
		#all bootstrap replicates were done by a single node, so loadingh is easy.
		print(paste("Re-analyzing existing bootstrap replicates with",ifelse(calculate_separate_tech,"covariate-specific spike-ins","common spike-ins")))
		trend = fread(paste0("./data/gene_variability/trend_",cur_group,suffix,".csv.gz"))
	}else{
		#we need to merge files from across replcates
		path = paste0("./data/gene_variability/",boot_file_prefix);
		fls = list.files(path=path,
				 pattern=paste0("trend_",cur_group,suffix,"_boot*"))
		trend = rbindlist(lapply(fls,function(x){
			fread(paste0(path,x))
			}));
		if (length(unique(trend$BootIndex)) < opt$nboots)
			stop(paste("Only",nrow(trend),"bootstrap replicates could be loaded from disk. Filenames: ",paste(fls,collapse="\n")))
			
		num_per_boot = aggregate(Group~BootIndex,trend,FUN=length);
		if (length(unique(num_per_boot$Group)) != 1){
			print(num_per_boot);
			warning("An uneven number of rows were found across different bootstrap replicates!")
		}
		
		fwrite(trend, paste0("./data/gene_variability/trend_",cur_group,suffix,".csv.gz"))
	}
}

if (!opt$reanalyze_boot){
	param_fun = function(x)rbindlist(lapply(x,`[[`, "params"))
	params = rbindlist(lapply(resdisp_parts, param_fun))

	print("Merging estimate table.")
	params <- merge(params, grouped_annots, by = "Group")
	setkey(params, GeneName, BootIndex, Group)
	fwrite(params, paste0("./data/gene_variability/",boot_file_prefix,"residual_overdispersion_",cur_group,suffix,boot_file_suffix,".csv.gz"))
	rm(resdisp_parts, envir = .GlobalEnv)
	resdisp_parts <<- NULL
	gc()
}else{
	if (opt$split_bootstraps_across_nodes == 1){
		#all bootstrap replicates were done by a single node, so loadingh is easy.
		params = fread(paste0("./data/gene_variability/residual_overdispersion_",cur_group,suffix,".csv.gz"))
	}else{
		#we need to merge files from across replcates
		path = paste0("./data/gene_variability/",boot_file_prefix);
		fls = list.files(path=path,
				 pattern=paste0("residual_overdispersion_",cur_group,suffix,"_boot*"))
		params = rbindlist(lapply(fls,function(x){
			fread(paste0(path,x))
		}));
		if (length(unique(params$BootIndex)) < opt$nboots)
			stop(paste("Only",nrow(params),"bootstrap replicates could be loaded from disk. Filenames: ",paste(fls,collapse="\n")))
			
		num_per_boot = aggregate(Group~BootIndex,params,FUN=length);
		if (length(unique(num_per_boot$Group)) != 1){
			print(num_per_boot);
			warning("An uneven number of rows were found across different bootstrap replicates!")
		}
		fwrite(params, paste0("./data/gene_variability/residual_overdispersion_",cur_group,suffix,".csv.gz"))
	}
	params[, Day := factor(Day, levels(annots$Day))]
	params[, Strain := factor(Strain, levels(annots$Strain))]
	params[, Food := factor(Food, levels(annots$Food))]
	params[, Temperature := factor(Temperature, levels(annots$Temperature))]
}

if (opt$split_bootstraps_across_nodes > 1 && !opt$reanalyze_boot){
	print("Since multiple nodes per bootstrap were specified, no further analysis is run by this node.")
	print("Once all nodes are complete, re-run script with --reanalyze_boot set")
	
}else{
	#set up factor levels for differential analysis
	if (!opt$differential_analysis_comparing_da_model){
		tmp = set_reference_groups_from_model(models[models$Name==cur_group,],
							     params[,.(Strain, Day, Food, Temperature,Pooled)]
						      );

		params[,Strain:=tmp$Strain][,Day:=tmp$Day][,Food:=tmp$Food][,Temperature:=tmp$Temperature][,Pooled:=tmp$Pooled]
		rm(tmp)
	}else{
		params$da_model_residuals = NA;
		params$da_model_residuals[grepl("\\.no_da_model",params$Group)] = "no_da_model";
		params$da_model_residuals[grepl("\\.da",params$Group)] = "da_model";
		if(any(is.na(params$da_model_residuals)))
			stop("Could not parse da model application pattern");
		params$da_model_residuals = factor(params$da_model_residuals,levels=c("no_da_model","da_model"))
		if (length(unique(params$da_model_residuals))<2)
			stop("Could not find both da and non-da groups in the parameter files.")
		
	}

	# summarize bootstrap results and compute HVG p-values

	# summarize trend results
	print("Summarizing Trend Info.")
	trend_df <- summarizeTrend(trend)
	fwrite(trend_df, paste0("./data/gene_variability/summary_trend_",cur_group,suffix,".csv.gz"))
	rm(trend)
	gc()

	print("Summarizing estimates.")
	params_df <- summarizeResVar(params)
	fwrite(params_df, paste0("./data/gene_variability/summary_residual_overdispersion_",cur_group,suffix,".csv.gz"))


	print("Running Differential Analysis.")
	if (!opt$differential_analysis_comparing_da_model){
		form <- paste("~", cur_bio_covariate)
	}else form <- paste("~da_model_residuals")
	if (batch_correct) 
		form = paste(form,"+ BiologicalReplicate")
	form = as.formula(form);
	if (!calculate_separate_tech){
		#params include both separate and consensus spike-ins
		#only use the consensus ones
		params = params[SequenceType != "Technical",]
	}
	
	
	d8_dva_df <- runDifferentialVariabilityAnalysis(params[bad_fit==F,], form)
	fwrite(d8_dva_df, paste0("./data/gene_variability/regression_effects_",cur_group,suffix,".csv.gz"))

}
print("Done!")