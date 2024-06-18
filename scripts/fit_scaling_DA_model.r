library(dplyr)
library(data.table)
library(pbapply)
library(ggrepel)
library(optparse)
library(stringr)
library(qs)


run_default_settings = F
skip_data_loading = F
if (run_default_settings){
	setwd(r"(Z:\nstroustrup\projects\asynch-seq-2022)")
}else setwd("/users/nstroustrup/nstroustrup/projects/asynch-seq-2022")

source("scripts/helpers/gene_variability.R")
source("scripts/helpers/preprocess.R")

#load single worm count data
geneid <- fread("data/annotations/c_elegans.PRJNA13758.WS265.geneIDs.txt.gz")
geneid[V3 == "", V3 := V4]
geneid <- geneid[, .(GeneName = V2, GeneSymbol = V3, GeneType = V6)]
geneid <- geneid[! duplicated(GeneName)]
setkey(geneid, GeneName)


#max = 3818
#max = 870
#max = 1334
#max = 69
#max = 253
#max = 660
dummy_rep = 1
#for (dummy_rep in 1:172){

opt = list(da_model = "N2_vha11_Set3_Day8",
	   single_worm_model="d8_sw",
	   single_worm_model_subgroup="QZ0",
	   nboots=1,
	   recalc_center = F,
	   modelindex = -1,
	   restricted_subset = F,
	   gene_set="",
	   min_mean_gene_count=0,
	   min_min_gene_count = 30,
	   gene_set_name="test",
	   print_count=T
	   )	   
	   
	   
option_list = list(
   make_option(c("-d", "--da_model"), 
               type="character",
               default="", 
              help="Name of RNAi or intervention to use"),
   make_option(c("-s", "--single_worm_model"), 
               type="character",
               default="", 
              help="Name of the single worm model in the differential_analysis model specifications"), 
              
   make_option(c("-S", "--single_worm_model_subgroup"), 
               type="character",
               default="", 
              help="Name of the single worm model subgroup in the differential_analysis model specifications"),
              
       make_option(c("-G", "--gene_set"), 
		     type="character",
		     default="", 
      		help="CSV file with a column with the name GeneName containing the genes to use for fitting"),
      make_option(c("-C","--min_mean_gene_count"),
      			type="integer",
      			default=25,
      			help="Minimum average abundance threshold to select genes for fitting the model"),
      make_option(c("-D","--min_min_gene_count"),
      			type="integer",
      			default=0,
      			help="Minimum average abundance threshold to select genes for fitting the model"),
       make_option(c("-g", "--gene_set_name"), 
      	                     type="character",
      	                     default="", 
              help="label for the geneset being used for fitting"),
              
	 make_option(c("-n", "--nboots"), 
		     type="integer",
		     default=NULL, 
      help="Number of bootstrap replicates."),
      
       make_option(c("-R", "--restricted_subset"), 
	      		     type="logical",
	      		     default=F, 
      help="only look at top hits"),
      
       make_option(c("-P", "--print_count"), 
      	      		     type="logical",
      	      		     default=F, 
      help="only look at top hits"),
      
	 make_option(c("-i", "--modelindex"), 
		     type="integer",
		     default=-1, 
      help="Index of model in differential analysis models."),
      
       make_option(c("-r", "--recalc_center"), 
      		     type="logical",
      		     default=F, 
      	help="recalculate center")
     )  
     
if (!run_default_settings){
	opt_parser <- OptionParser(option_list = option_list, prog = "Fit negative binomial model")
	opt <- parse_args(opt_parser)
}
#print(opt)
## First, figure out which RNAis/interventions to use

da_models_df <- fread("data/differential_analysis/models/net_validation.csv")
da_models_df = da_models_df[!grepl("uls60",Name)&!grepl("ts_",Name),]

#zero out all RNAis!
#da_models_df = da_models_df[rep(F,nrow(da_models_df))]
#da_models_df = da_models_df[c(rep(F,nrow(da_models_df)))]


da_int_models_df <- fread("data/differential_analysis/models/single_and_pooled_worms.csv")
da_int_models_to_load = c("d1_sw","d8_sw","n2_sw","d1_glp1","d8_glp1","d1_UV","d8_UV","d8_20v25C","d1_20v25C",
"n2_series_d2","n2_series_d4","n2_series_d10","n2_series_d2_d6","n2_series_d6_d10")

da_int_models_metadata = data.frame(model = da_int_models_to_load,
		  RNAi = c("daf-2_vs_N2","daf-2_vs_N2","N2 d8_vs_d1","glp-1_vs_N2",
		  "glp-1_vs_N2","Live_vs_UV","Live_vs_UV","25_vs_20C","25_vs_20C",
		  "n2_series day 2 vs 4","n2_series day 4 vs 6","n2_series day 10 vs 12","n2_series day 2 vs 6","n2_series day 6 vs 10"),
		  Day = c(1,8,8,1,8,1,8,3,8,8,1,8,2,4,10,2,6),
		  Genotype = rep("QZ0",length(da_int_models_to_load)), #everything is relative to wildtype (like each RNAi is relative to the EV)
		  Set = rep(100,length(da_int_models_to_load)))
rownames(da_int_models_metadata) = da_int_models_metadata$model;

#da_int_models_metadata = da_int_models_metadata[grepl("series",da_int_models_metadata$model),]

setkey(da_int_models_df,Name)
#da_int_models_df = da_int_models_df[c(da_int_models_metadata$model,"n2_series","n2_series_d8_ref"),]
extra_int_models_added_in = 0
da_int_models_df = da_int_models_df[c(da_int_models_metadata$model),]

#zero out int models
#da_int_models_df = da_int_models_df[rep(F,nrow(da_int_models_df)),]
#da_int_models_df = da_int_models_df[grepl(Name,"series_d"),]

n2_series_model_to_load = data.table(
	  Name = c("n2_d1_vs_8","n2_d2_vs_8","n2_d4_vs_8","n2_d6_vs_8","n2_d10_vs_8","n2_d12_vs_8",
		     "n2_d8_vs_1","n2_d2_vs_1","n2_d4_vs_1","n2_d6_vs_1",'n2_d10_vs_1','n2_d12_vs_1'),
          Source = c(rep("n2_series_d8_ref",6),rep("n2_series",6)),
	  Day = c(1,2,4,6,10,12,8,2,4,6,10,12),
	  Genotype = rep("QZ0",12),	
	  Set = rep(100,12)
	  )
#zero out all n2 series
#n2_series_model_to_load = n2_series_model_to_load[c(rep(F,nrow(n2_series_model_to_load))),]

rnai_models_to_run = c();#da_models_df$Name;
int_models_to_run = c();#da_int_models_df$Name
n2_series_models_to_run = c();#n2_series_model_to_load$model
count_set = ""
min_mean_gene_count = opt$min_mean_gene_count
min_min_gene_count = opt$min_min_gene_count
if (min_min_gene_count > 0)
	print(paste("Restricting genes using a minimum minimum count threshold of ",min_min_gene_count))
if (min_mean_gene_count > 0)
	print(paste("Restricting genes using a minimum mean count threshold of ",min_mean_gene_count))

if (opt$restricted_subset==T){
	#only use subset!
	
	to_use = c("d1_glp1",  "d8_glp1",  "N2_acox1.1_Set1_Day8",  "N2_aexr1_Set3_Day8",  "N2_algn14_Set3_Day8",  "N2_atp5_Set5_Day8", 
		"N2_C08H9.15_Set3_Day8",  "N2_C50B6.7_Set1_Day8",  "N2_C50F4.4_Set3_Day8",  "N2_ccch1_Set3_Day8",  "N2_col122_Set3_Day8",  
		"N2_cyp13A10_Set1_Day8",  "N2_ech4_Set1_Day8",  "N2_eef1A.2_Set5_Day8",  "N2_egl3_Set4_Day8",  "N2_F16C3.2_Set3_Day8",  
		"N2_F41C3.5_Set4_Day8",  "N2_flp27_Set3_Day8",  "N2_flp5_Set3_Day8",  "N2_fmo5_Set3_Day8",  "N2_glc1_Set2_Day8",  
		"N2_idh1_Set4_Day8",  "N2_mak1_Set3_Day8",  "N2_mig6_Set3_Day8",  "N2_msrp7_Set5_Day8",  "N2_nlp15_Set3_Day8", 
		"N2_nlp28_Set3_Day8",  "N2_pap1_Set5_Day8",  "N2_pcn1_Set5_Day8",  "N2_pgal1_Set4_Day8",  "N2_pghm1_Set4_Day8", 
		"N2_rop1_Set1_Day8",  "N2_rpl6_Set5_Day8",  "N2_rpl7_Set5_Day8",  "N2_sbt1_Set4_Day8",  "N2_sdc2_Set2_Day8",
		"n2_series_d2",  "N2_srp2_Set3_Day8",  "N2_T08B2.11_Set3_Day8",  "N2_vha11_Set3_Day8")  
	#to_use = c("d1_glp1","d8_glp1","n2_series_d2","n2_series_d4","n2_series_d2_d6","n2_series_d6_d10","n2_series_d10","N2_nlp28_Set3_Day8","N2_eef1A.2_Set5_Day8","N2_pcn1_Set5_Day8","N2_egl3_Set4_Day8","N2_vha11_Set3_Day8")
	#to_use=c()
	da_models_df = da_models_df[ Name %in% 	to_use,]	
	#print(dim(da_int_models_df))
	#print(da_int_models_df$Name)	
	da_int_models_df = da_int_models_df[Name %in% 	to_use,]
	#print(dim(da_int_models_df))	
	#use all n2 series models
	n2_series_model_to_load = n2_series_model_to_load[Name %in% 	to_use,]
	#browser()
}
if (opt$da_model != ""){

	if (opt$da_model %in% da_models_df$Name){
		rnai_models_to_run = opt$da_model
		n2_series_da_models_to_run = c();
	}else{
		rnai_models_to_run = c();
		if (opt$da_model %in% da_int_models_to_load){
			rnai_models_to_run = c();
			int_models_to_run = opt$da_model
			n2_series_da_models_to_run = c();
		}else{
			
			int_models_to_run = c();
			rnai_models_to_run = c();
			if (opt$da_model %in% n2_series_model_to_load$model){
				n2_series_da_models_to_run = opt$da_model;
			}else stop(paste("Could not find specified model",opt$da_model))
		}
	}
	
}else{
	#we run all models
	rnai_models_to_run = da_models_df$Name;
	int_models_to_run = da_int_models_df$Name
	n2_series_da_models_to_run = n2_series_model_to_load$model
}
#models = fread("data/differential_analysis/models/single_and_pooled_worms.csv")

single_worm_data_possibilities =  
	list(
		d8_UV=list(Food=c("NEC937","NEC937_UV")),
		#d1_UV=list(Food=c("NEC937","NEC937_UV")),
		#n2_series_all = list(Day=c(1,2,4,6,8,10,12))#,
		#n2_series_all = list(Day=c(1,8))#, just 1 and 8
		d8_glp1=list(Strain=c("QZ0","CB4037")),
		#d1_glp1=list(Strain=c("QZ0","CB4037")),
		#d1_sw=list(Strain=c("QZ0","QZ120")),
		#d8_sw=list(Strain=c("QZ0","QZ120"))#,
		d8_20v25C=list(Temperature=c("20","25"))
		#d1_20v25C=list(Temperature=c("20","25"))
	)
if(0 & opt$restricted_subset==T){
	cat("Using restricted subset\n")
	to_use = c("n2_series","d8_sw","d8_sw","d8_glp1","d1_glp1")
	single_worm_data_possibilities = single_worm_data_possibilities[names(single_worm_data_possibilities)[ ! names(single_worm_data_possibilities)%in% to_use]]
}

if (opt$single_worm_model!= ""){
	if (!opt$single_worm_model %in% names(single_worm_data_possibilities))
		stop("Invalid single worm model")
	single_worm_data_to_use = single_worm_data_possibilities[opt$single_worm_model]
	
	if (opt$single_worm_model_subgroup != ""){
		if (!opt$single_worm_model_subgroup %in% single_worm_data_possibilities[[opt$single_worm_model]][[1]])
			stop("Invalid single worm model subgroup")
		#restrict the sbugroups to a single value
		 single_worm_data_to_use[[opt$single_worm_model]][[names(single_worm_data_to_use[[opt$single_worm_model]])]] = opt$single_worm_model_subgroup
	}
}else{
	#we test all single worm models
	single_worm_data_to_use = single_worm_data_possibilities
}
subgroup_counts = sapply(single_worm_data_to_use,function(x)length(x[[1]]))
set_lookup = cumsum(subgroup_counts)-subgroup_counts

total_single_worm_sets = sum(subgroup_counts)

RNAi_N = length(da_models_df$Name);
int_N = length(da_int_models_df$Name)

if (opt$modelindex != -1 || opt$print_count){
	#figure out which DA model to use--it can either be an RNAi, an intervention (int), or a n2_series comparrison
	da_model_num = floor((opt$modelindex-1)/total_single_worm_sets)+1;
	
	num_model_sets = RNAi_N+int_N-extra_int_models_added_in+nrow(n2_series_model_to_load)
	if (opt$print_count){
		print(paste0("Max model number: ",num_model_sets," models * ", total_single_worm_sets, " single worm cohorts = ",num_model_sets*total_single_worm_sets))
		print(paste0("(Model sets: ",RNAi_N," RNAi,  ", int_N, " interventions, ",nrow(n2_series_model_to_load), " n2_series)"))
		stop()
	}
	if (da_model_num <= RNAi_N){
		rnai_models_to_run = da_models_df$Name[da_model_num]
		n2_series_da_models_to_run = c();
		int_models_to_run = c();
	}else{
		rnai_models_to_run = c();
		if (da_model_num <= RNAi_N+int_N-extra_int_models_added_in){
			int_models_to_run = da_int_models_df$Name[da_model_num-RNAi_N]
			n2_series_da_models_to_run = c();
		}else {
			int_models_to_run = c();
			if (da_model_num <= num_model_sets){
				n2_series_da_models_to_run = n2_series_model_to_load$Name[da_model_num-RNAi_N-int_N+extra_int_models_added_in]
			} else stop("Invalid model number")
		}
	}
	#now figure out what single worm data set to use
	single_worm_model_num = (opt$modelindex-1)%%total_single_worm_sets+1
	tmp = which(set_lookup<single_worm_model_num)
	single_worm_model_name = names(tmp)[length(tmp)]
	single_worm_submodel_num = single_worm_model_num-set_lookup[single_worm_model_name]
	single_worm_submodel_name = single_worm_data_to_use[[single_worm_model_name]][[1]][single_worm_submodel_num]
	
	if (single_worm_model_name %in% unique(n2_series_model_to_load$source)){
		n2_series_single_worm_model = n2_series_model_to_load[Source == single_worm_model_name & Day==single_worm_submodel_name,Name]
	}else n2_series_single_worm_model = ""
}else{
#rnai_models_to_run = da_models_df$Name;
#	int_models_to_run = da_int_models_df$Name
#	n2_series_da_models_to_run = n2_series_model_to_load$model
	n2_series_single_worm_model = ""



	#rnai_models_to_run = da_models_df$Name
	#int_models_to_run = da_int_models_df$Name
	#n2_series_da_models_to_run = da_int_models_df$Name
	#n2_series_single_worm_model = ""
	single_worm_model_name = opt$single_worm_model
	single_worm_submodel_name = opt$single_worm_model_subgroup;
	single_worm_submodel_num = which((single_worm_data_to_use[[single_worm_model_name]][[1]]) == single_worm_submodel_name)
}
print(paste("Processing ",opt$modelindex,paste(rnai_models_to_run,collapse=","),paste(int_models_to_run,collapse=","),n2_series_da_models_to_run,single_worm_model_name,single_worm_submodel_name,n2_series_single_worm_model))

#}
#stop()

#OK! We now have rnai_models_to_run, int_models_to_run, single_worm_model_name, single_worm_submodel_name, n2_series_single_worm_model, and n2_series_da_models_to_run
#and we can process all of them!

genes_to_use = c();
gene_set_name = ""
if (opt$gene_set != ""){
	genes_to_use = fread(opt$gene_set)$GeneName;
	if (opt$gene_set_name != ""){
		gene_set_name = opt$gene_set_name
	}else 	gene_set_name = "custom_gene_set";
	cat(paste0("Fitting using custom gene set ",gene_set_name," containing ",length(genes_to_use)," genes.\n"))
}else{
	if (opt$gene_set_name != "")
		gene_set_name = opt$gene_set_name
}
#Load in differential analysis data

cat("Loading DA data...\n")
if (length(rnai_models_to_run) == 0){
	da_df = data.frame();
}else{	
	if (!skip_data_loading){
		da_df <- rbindlist(lapply(rnai_models_to_run, function(model) {
		  out <- fread(paste0("data/differential_analysis/tables/net_validation/", model,".csv.gz"))
		  colnames(out)[3:ncol(out)] <- c("log2FoldChange", "lfcSE", "pvalue", "padj")
		  out[, Name := model]
		  out
		}))
		da_df <- merge(da_models_df[, .(Name, RNAi, Genotype,Set,Day)], da_df, by = "Name")


		 RNAi_name_fix = matrix(c("C15C74", "C15C7.4"   ,
		  "K02F25",   "K02F2.5" ,
		  "R1022",   "R102.2" ,
		  "Y44A6D2",  "Y44A6D.2", 
		  "Y9C2VA1", "Y9C2UA.1" ,
		  "casy1", "casy-1"   ,
		  "egl21","egl-21"    ,
		  "flp12",   "flp-12" ,
		  "flp14",    "flp-14" , 
		  "flp1","flp-1",
		  "flp5","flp-5",
		  "flp9","flp-9",
		  "ida1","ida-1",
		  "ins24","ins-24",
		  "ins26","ins-26",
		  "ins30","ins-30",
		  "ins6","ins-6",
		  "mec12","mec-12",
		  "nlp15","nlp-15",
		  "nlp21","nlp-21",
		  "nlp50","nlp-50",
		  "pdf1","pdf-1",
		  "pgal1","pgal-1",
		  "pghm1","pghm-1",
		  "sbt1","sbt-1",
		  "snet1","snet-1",
		  "snt4","snt-4",
		  "ttr29","ttr-29",
		  "zig2","zig-2"),ncol=2,byrow=T)

		  for ( i in 1:nrow(RNAi_name_fix)){da_df[RNAi==RNAi_name_fix[i],]$RNAi = RNAi_name_fix[i,2]}
	}
}

#browser()
annots = fread("data/annotations/sample_annotations_single_and_pooled_worms.csv")

#Load in differential analysis for all other interventiosn (glp-1, aging, etc)
if (length(int_models_to_run) == 0){
	da_int_df = data.frame();
}
tmp = list()
if (length(n2_series_da_models_to_run)>0){
	setkey(n2_series_model_to_load,Name)
	#series_model_to_load$RNAi = series_model_to_load$model
	
	if ("n2_series_d8_ref" %in% n2_series_model_to_load[n2_series_da_models_to_run,Source]){
		#load timseries data
		 out <- fread(paste0("data/differential_analysis/tables/single_and_pooled_worms/n2_series_d8_ref.csv.gz"))
		 tmp[["n2_d1_vs_8"]] = out[,1:6] # 1 vs 8
		 tmp[["n2_d2_vs_8"]] = out[,c(1:2,7:10)] # 2 vs 8
		 tmp[["n2_d4_vs_8"]] = out[,c(1:2,11:14)] # 4 vs 8
		 tmp[["n2_d6_vs_8"]] = out[,c(1:2,15:18)] # 6 vs 8
		 tmp[['n2_d10_vs_8']] = out[,c(1:2,19:22)] # 10 vs 8
		 tmp[['n2_d12_vs_8']] = out[,c(1:2,23:26)] # 12 vs 8

		for (n in names(tmp)){
		  colnames(tmp[[n]])[3:6] <- c("log2FoldChange", "lfcSE", "pvalue", "padj")
		  tmp[[n]][,Name:=n]
		}
	}
	if ("n2_series" %in% n2_series_model_to_load[n2_series_da_models_to_run,Source]){
		 out <- fread(paste0("data/differential_analysis/tables/single_and_pooled_worms/n2_series.csv.gz"))
		 tmp[["n2_d8_vs_1"]] = out[,1:6] # 1 vs 8
		 tmp[["n2_d2_vs_1"]] = out[,c(1:2,7:10)] # 2 vs 8
		 tmp[["n2_d4_vs_1"]] = out[,c(1:2,11:14)] # 2 vs 8
		 tmp[["n2_d6_vs_1"]] = out[,c(1:2,15:18)] # 2 vs 8
		 tmp[['n2_d10_vs_1']] = out[,c(1:2,19:22)] # 2 vs 8
		 tmp[['n2_d12_vs_1']] = out[,c(1:2,23:26)] # 2 vs 8
		for (n in names(tmp)){
		  colnames(tmp[[n]])[3:6] <- c("log2FoldChange", "lfcSE", "pvalue", "padj")
		  tmp[[n]][,Name:=n]
		}
	}
	tmp = tmp[n2_series_da_models_to_run];
}

if (length(int_models_to_run) > 0){
	tmp2 = lapply(int_models_to_run, function(model) {
	  out <- fread(paste0("data/differential_analysis/tables/single_and_pooled_worms/", model, ".csv.gz"))
	  colnames(out)[3:6] <- c("log2FoldChange", "lfcSE", "pvalue", "padj")
	  out = out[,1:6]
	  out[, Name := model]
	  out;
	})
	names(tmp2) = int_models_to_run
	tmp = c(tmp,tmp2)
	rm(tmp2)
	
}
if (length(tmp) > 0){
	da_int_df <- rbindlist(tmp,fill=T)

	#add in metadata
	series_data = grepl("n2_d",da_int_df$Name)
	if (any(!series_data)){
		da_int_df$RNAi = da_int_models_metadata[da_int_df$Name,"RNAi"]
		da_int_df$Day = da_int_models_metadata[da_int_df$Name,"Day"]
		da_int_df$Genotype = da_int_models_metadata[da_int_df$Name,"Genotype"]
		da_int_df$Set = da_int_models_metadata[da_int_df$Name,"Set"]
	}

	if (any(series_data)){
		setkey(n2_series_model_to_load,Name)
		da_int_df$RNAi[series_data] = n2_series_model_to_load[da_int_df$Name[series_data],]$Name
		da_int_df$Day[series_data] = n2_series_model_to_load[da_int_df$Name[series_data],]$Day
		da_int_df$Genotype[series_data] = n2_series_model_to_load[da_int_df$Name[series_data],]$Genotype
		da_int_df$Set[series_data] = n2_series_model_to_load[da_int_df$Name[series_data],]$Set
	}
}

if (nrow(da_int_df) > 0 && nrow(da_df) > 0 ){
	da_df = rbind(da_int_df,da_df)
}else if (nrow(da_int_df) > 0){
	da_df = da_int_df;
}else if (nrow(da_df) > 0){
	da_df = da_df
}else stop("No DA data!")


#we add a few more gene names to exclude,
#because these have additional off-target RNAi artifacts
da_df$extra_RNAi_artifacts_to_exclude = NA;
if(any(grepl("col122",da_df$Name)))
da_df$extra_RNAi_artifacts_to_exclude[grepl("col122",da_df$Name)] = "col-122";
if(any(grepl("vha11",da_df$Name)))
	da_df$extra_RNAi_artifacts_to_exclude[grepl("vha11",da_df$Name)] = "vha-3";
#browser()
#we have all the DA data now.
cat("Loading counts...\n")
if (!skip_data_loading){
	bc <- readRDS("data/batch_corrected_counts/single_and_pooled_worms.rds") 
	if (!is.null(single_worm_model_name)){
		single_worm_data_to_use = single_worm_data_to_use[single_worm_model_name]
	}
	if (single_worm_submodel_name != ""){
		single_worm_data_to_use[[single_worm_model_name]][[1]] = single_worm_data_to_use[[single_worm_model_name]][[1]][single_worm_submodel_num]
	}

	#make sure no ERCCs have snuck in from somewhere
	counts = lapply(names(single_worm_data_to_use),function(nm){
		bc$batch_counts_list[[nm]][!grepl("ERCC",rownames(bc$batch_counts_list[[nm]])),]
	})
	names(counts) = names(single_worm_data_to_use)
	rm(bc);
	gc()

	setkey(annots,Sample)

	nf <- lapply(names(counts),function(nm){
		#if (min_min_gene_count > 0){
		#	genes_for_normalization <- apply(counts[[nm]],1,function(cnts2)!any(cnts2 < min_min_gene_count))
		#	cat(paste("Normalizing using",length(which(genes_for_normalization)),"of", nrow(counts[[nm]]),"genes...\n"))
		#}
		#normalizationFactor(counts[[nm]][genes_for_normalization,], groups = annots[colnames(counts[[nm]]),paste(Strain,Temperature,Day,Food)])
		normalizationFactor(counts[[nm]], groups = annots[colnames(counts[[nm]]),paste(Strain,Temperature,Day,Food)])
		})
	names(nf) = names(counts);
	#browser()


	norm <- lapply(names(counts),function(nm){
		normalizeCounts(counts[[nm]], nf[[nm]])
		})
	names(norm) = names(counts);
	if (!skip_data_loading)
		rm(counts)
	gc()
}

for (current_da in unique(da_df$Name)){
	da = da_df[Name == current_da,]
	if (length(genes_to_use)>0){
		common_genes = intersect(genes_to_use,da$GeneName)
		da = da[GeneName %in% common_genes,]
	}
	setkey(geneid,GeneSymbol)
	RNAi_GeneName_to_exclude = c(geneid[unique(da$RNAi),GeneName],geneid[unique(da$extra_RNAi_artifacts_to_exclude),GeneName])
	RNAi_GeneName_to_exclude = RNAi_GeneName_to_exclude[!is.na(RNAi_GeneName_to_exclude)]
	#browser()
	for(current_single_worm_model in names(single_worm_data_to_use)){
		setkey(annots,Sample)
		cur_annots = annots[colnames(norm[[current_single_worm_model]]),]
		if (current_single_worm_model %in% c("n2_series","n2_series_d8_ref","n2_series_all","n2_series_d8_ref_all")){
			cur_annots = cur_annots[grepl("TS_QZ0",Sample),]  #use only the single TS replicate shared among all timepoints--d8 and d1 by default contain all animals from all reps
		}
		current_single_worm_submodel_column = names(single_worm_data_to_use[[current_single_worm_model]])
	
		lapply(single_worm_data_to_use[[current_single_worm_model]][[1]],function(current_single_worm_submodel){
			samples = cur_annots[get(current_single_worm_submodel_column)==current_single_worm_submodel,Sample]
			cur_counts = norm[[current_single_worm_model]][,samples]
			print(paste("Fitting a population of ",ncol(cur_counts),"individuals"))
			if (length(genes_to_use)>0){
				common_genes = intersect(common_genes,rownames(cur_counts))
				cur_counts = cur_counts[common_genes,]
				cat(paste0("DA has ",nrow(da)," genes; counts have ", nrow(cur_counts), " genes\n"));
			}
			#cur_counts <- cur_counts[apply(cur_counts,1,function(x)!any(x < 30)), ]	#performs poorly on counts with low expression
			#cat(paste0(nrow(cur_counts)," genes remain after abundance check\n"))
			#print(unique(cur_annots[get(current_single_worm_submodel_column)==current_single_worm_submodel,get(current_single_worm_submodel_column)]))
			
			filename = paste(current_da,current_single_worm_model,current_single_worm_submodel,sep="_")
			if (length(genes_to_use)==0 && gene_set_name == ""){
				path = "data/scaling_models/all_grid/"
			}else{
				path = paste0("data/scaling_models/",gene_set_name,"/")
			}
			dir.create(path, showWarnings = F)
			
			all_results_fname = paste0(path,filename,"=all",".qs")
			stats_fname = paste0(path,filename,"=stats",".qs")
			boot_fname = paste0(path,filename,"=boots",".qs")
			
			if (!opt$recalc_center){
				cat(paste("Bootstrapping against",current_single_worm_model,current_single_worm_submodel_column,current_single_worm_submodel,"...\n"))
				
				 res = ns_fit_da_linear_model_bootstrapped(opt$nboots,
							   cur_counts,
							   da,
							   min_mean_count_for_fitting=min_mean_gene_count,
							   min_min_count_for_fitting=min_min_gene_count,
							   differential_analysis_column_name = "log2FoldChange",
							   RNAi_GeneName_to_exclude=RNAi_GeneName_to_exclude,
							   allow_recursion_for_outliers=2,
						   keep_x=F,
						   keep_residuals=F
						)
				res$da = current_da	
				res$min_min_gene_count = min_min_gene_count
				res$single_worm_model = current_single_worm_model
				res$single_worm_submodel = current_single_worm_submodel
				
				#save different data sources separately to speed up loading
							
				cat(paste("Saving to:",all_results_fname,"\n",stats_fname,"\n",boot_fname,"\n"))
				
				qsave(res[c("dat","stats","da","single_worm_model","single_worm_submodel","outliers_removed","min_min_gene_count")],all_results_fname)
			
				qsave(res[["boot_data"]],boot_fname)
			}else{
				cat(paste("Updating residuals for ",current_single_worm_model,current_single_worm_submodel_column,current_single_worm_submodel,"...\n"))
				res = list()
				res[["stats"]] = qread(stats_fname)
				cat(paste("Saving to:",stats_fname,"\n"))
			}
			p = as.data.frame(res[["stats"]]);
			p$da = res$da
			p$single_worm_model = res$single_worm_model
			p$single_worm_submodel = res$single_worm_submodel
			qsave(p,stats_fname)
			print(paste("Change in relative variance:",res$stats["relative_per_gene_variance","mean_est"],"\n"))
			
		})
	}
}
print("Done")
 