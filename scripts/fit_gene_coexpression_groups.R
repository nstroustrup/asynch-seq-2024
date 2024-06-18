library(dplyr)
library(data.table)
library(optparse)
library(pbapply)
library(dendextend)
library(qs)
#library(ggplot2)
#library(cowplot)
#library(ComplexHeatmap)
#library(BASiCS)
#library(dendextend)
#library(ggrepel)
#library(reshape2)
#library(ggrastr)
#require(scales)
#library(RColorBrewer)
#library(gridExtra)
#library(circlize)
#Run community detection
library(igraph)
#theme_set(theme_cowplot())


#source("scripts/helpers/clustering.R")
#source("scripts/notebooks_helpers/gene_network.R")


run_default_settings= F
if (run_default_settings){
  setwd("Z:/nstroustrup/projects/asynch-seq-2022")
}else 
  setwd("/users/nstroustrup/nstroustrup/projects/asynch-seq-2022")
source("scripts/helpers/gene_variability.R")
source("scripts/helpers/preprocess.R")
## Load data

### Gene names

geneid <- fread("./data/annotations/c_elegans.PRJNA13758.WS265.geneIDs.txt.gz")
geneid[V3 == "", V3 := V4]
geneid <- geneid[, .(GeneName = V2, GeneSymbol = V3, GeneType  = V6)]
geneid <- geneid[! duplicated(GeneName)]


### Gene annotations

gene_annotations <- readRDS("./data/annotations/gene_annotations.rds")
gene_annotations <- gene_annotations[c("Wormcat", "Wormexp", "KEGG")]


opt <- list(model = "n2_sw",
            model_index = -1,
            tissue_normalization = F,
            compile_and_label = F,
            da_model = "da_d8_glp1",
            da_model_comparison=F
)
if (0)
  opt <- list(model = "",
              model_index = 1,
              tissue_normalization = F,
              compile_and_label = F,
              da_model = "",
              da_model_comparison=F
  )
if (1)
opt <- list(model = "",
            model_index = 21,
            tissue_normalization = F,
            compile_and_label = T,
            da_model = "",
            da_model_comparison=F
)


# handle arguments
option_list = list(
  
  make_option(c("-m", "--model"), 
              type="character",
              default="", 
              help="Specific model to process"),
  
  make_option(c("-i", "--model_index"), 
              type="numeric",
              default=-1, 
              help="Specific model to process, by index"),
  
  make_option(c("-t", "--tissue_normalization"), 
              type="logical",
              default=F, 
              help="run tissue_specific analysis"),
  
  make_option(c("-r", "--compile_and_label"), 
              type="logical",
              default=F, 
              help="compile and label previous reults"),
  
  make_option(c("-d", "--da_model"), 
              type="character",
              default="", 
              help="compile and label previous reults"),
  make_option(c("-D", "--da_model_comparison"), 
              type="logical",
              default=F, 
              help="Compare covariance before and after da model")
)

if (!run_default_settings){
  opt_parser <- OptionParser(option_list = option_list, prog = "Fit co-expression groups")
  opt <- parse_args(opt_parser)
}
build_communities = !opt$compile_and_label
if (opt$model == "" && opt$model_index == -1 && !opt$compile_and_label){
  cat("No model or model index were specified, so we will run and analyze everything.\n");
  opt$compile_and_label = T;
}
if ((opt$model != "" || opt$da_model != "") && opt$model_index != -1)
  stop("Either model/da_model or model_index should be specified, not both")
  
if ((opt$model != "" || opt$da_model != "") &&  opt$compile_and_label)
  stop("Model specified but compile and label was also specified: nothing to be done!")

### Set up model comparison structures to guide analysis

models = read.csv("./data/differential_analysis/models/single_and_pooled_worms.csv");

all_model_names = c("d8_glp1","d1_glp1","n2_sw","daf2_sw","d1_sw","d8_sw")

setkey(geneid, GeneSymbol)

da_models = data.table(name = c("da_nlp28","da_atp5","da_d8_glp1"),
                       filename = c("./data/differential_analysis/tables/net_validation/N2_nlp28_Set3_Day8.csv.gz",
                                    "./data/differential_analysis/tables/net_validation/N2_atp5_Set5_Day8.csv.gz",
                                    "./data/differential_analysis/tables/single_and_pooled_worms/d8_glp1.csv.gz"),
                       colName = c("RNAi_nlp28_vs_EV_log2FoldChange","RNAi_atp5_vs_EV_log2FoldChange","Strain_CB4037_vs_QZ0_log2FoldChange"),
                       GeneName = c(geneid["nlp-28",GeneName],geneid["atp-5",GeneName],""))
da_model_name = ""
if (opt$model != ""){
  if (! opt$model %in% all_model_names)
    stop("Model name not found in list")
  model_names = opt$model
  da_model_name = opt$da_model
  output_suffix = opt$model;
}else if (opt$model_index != -1){
  if (opt$model_index > (nrow(da_models)+1)*length(all_model_names) || opt$model_index < 1)
    stop("Invalid model index")
  model_names = all_model_names[[((opt$model_index-1) %% length(all_model_names))+1]];
  da_models_index = floor((opt$model_index-1)/length(all_model_names))
  output_suffix = model_names;
  if(da_models_index == 0){
    da_model_name = "";
  }else{ 
    da_model_name = da_models[da_models_index,name];
  }
}else{ 
  model_names = all_model_names;
}
if (da_model_name != ""){
  output_suffix = paste0(output_suffix,"_",da_model_name);
  if (opt$da_model_comparison)
    output_suffix = paste0(output_suffix,"_DC");
}

use_tissue_specific = opt$tissue_normalization
file_suffix = "";
if (use_tissue_specific){
  file_suffix = "_tissue_specific";
  cat("Using tissue normalized data.\n")
}else cat ("Using all Genes (not just tissue specific)")



coexp_ID_from_number = function(i){
  sapply(i,function(x){ifelse(x<10,paste0("G0",x),paste0("G",x))})
}

makeGraph <- function(rho, theta, beta = 2) {
  adj <- abs(rho)^beta
  adj[upper.tri(adj, diag = TRUE)] <- 0
  adj[adj < theta] <- 0
  gr <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE, weighted = TRUE)
}

tissues_df <- fread("./data/tissue_unique_genes/genes_unique_to_tissues.csv.gz")
tissues <- pull(tissues_df, "Tissue", "GeneName")

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

if (build_communities){
  cat(paste("Building communities for",paste0(model_names,collapse=",")))
  
  if (da_model_name != ""){
    cat(paste(" using da model",da_model_name))
  }
  cat("\n")
  #make a list of which samples match up to each subgroup (biological covariate group)
  #load counts
    cat("Loading counts...\n")
    sw <- readRDS("./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds")
    annotations = sw$annots;
    samples_list = lapply(model_names,function(group){
      cur_annots = sw$annots[eval(parse(text=models$Subset[models$Name == group])),]
      subgroup_vals = as.character(as.data.frame(cur_annots)[,models$BiologicalCovariates[models$Name == group]])
     # browser()
      subgroups = unique(subgroup_vals)
      res = lapply(subgroups,function(subgroup){
          cur_annots$Sample[subgroup_vals == subgroup]
      }) 
      names(res) = subgroups
      res
    })
    names(samples_list) = model_names;
  if (0){#load raw counts
    # counts_list <- sapply(samples_list, function(samples) sw$counts[, samples])
    counts_list <- lapply(names(samples_list),function(x){
        res = lapply(names(samples_list[[x]]),function(y){
            sw$counts[, colnames(sw$counts)%in% samples_list[[x]][[y]]]
        })
        names(res) = names(samples_list[[x]])
        res;
    })
    names(counts_list) = names(samples_list)
    rm(sw)
  }else{#load batch-corrected counts
    rm(sw)
    bc <- readRDS("./data/batch_corrected_counts/single_and_pooled_worms.rds") 
    counts_list <- lapply(model_names,function(x){
        res = lapply(names(samples_list[[x]]),function(y){
   
           bc$batch_counts_list[[x]][, colnames(bc$batch_counts_list[[x]])%in% samples_list[[x]][[y]]]
      })
      names(res) = names(samples_list[[x]])
      res;
    })
    names(counts_list) = model_names
    rm(bc)
  }
  cat("Performing normalization...\n")
  nf <- lapply(model_names,function(group){
     ret = lapply(counts_list[[group]],function(x){
       normalizationFactor(x, groups = annotations[Sample %in% colnames(x)]$Group)
     })
     names(ret) = names(counts_list[[group]])
     ret;
  })
     
  names(nf) = model_names;
  
  norm_list <- lapply(model_names,function(group){
      res = lapply(names(counts_list[[group]]),function(subgroup){
      
        normalizeCounts(counts_list[[group]][[subgroup]], nf[[group]][[subgroup]])
      })
      names(res) = names(counts_list[[group]])
      res
  })
  names(norm_list) = model_names;
  
  #normalize each tissue indepenently
  tissue_norm_list <- lapply(model_names, function(group) {
        ret = lapply(counts_list[[group]],function(x){
          k <- normalizationFactorTissueSpecific(x, tissues,annotations[Sample %in% colnames(x)]$Group)
          ret = normalizeCountsTissueSpecific(x, k, tissues)
          ret;
        })
        names(ret) = names(counts_list[[group]])
        ret;
  })
  names(tissue_norm_list) = model_names;
          
  
  if (use_tissue_specific){
    cv_to_use = tissue_norm_list
  }else{
    cv_to_use = norm_list
  }
  rm(tissue_norm_list)
  rm(norm_list)
  

  
  #n2 series needs to be split up into pairwise analyses
  if ("n2_series" %in% model_names){
    #split the series into a set of pairwise comparrisons
    
    cv_to_use[["n2_series_1_6"]] = list(`1` = cv_to_use[["n2_series"]][["1"]],`6` = cv_to_use[["n2_series"]][["6"]])
    cv_to_use[["n2_series_1_10"]] = list(`1` = cv_to_use[["n2_series"]][["1"]],`10` = cv_to_use[["n2_series"]][["10"]])
    cv_to_use[["n2_series_1_12"]] = list(`1` = cv_to_use[["n2_series"]][["1"]],`12` = cv_to_use[["n2_series"]][["12"]])
    cv_to_use[["n2_series_2_10"]] = list(`2` = cv_to_use[["n2_series"]][["2"]],`10` = cv_to_use[["n2_series"]][["10"]])
    cv_to_use[["n2_series_2_12"]] = list(`2` = cv_to_use[["n2_series"]][["2"]],`12` = cv_to_use[["n2_series"]][["12"]])
    cv_to_use[["n2_series"]] = NULL
    model_names = names(cv_to_use)
  }
  
  
  #add extra data to study the effect of removing nlp-28 effect from day 8 and atp-5 effect from day 1
  
  if (da_model_name != ""){
    setkey(geneid,GeneSymbol)
    da_info = da_models[name == da_model_name,]
    GeneName = da_info$GeneName
    da_data = fread(da_info$filename)
    
    for (gn in model_names){
      tmp = cv_to_use[[gn]]
      if (opt$da_model_comparison){
        new_name = paste0(gn,"_",da_model_name,"_DC")
        cat(paste0("Running DA model on ",gn," submodel ",names(tmp)[[1]],"...\n"))
        da_model_results = ns_fit_da_linear_model(tmp[[1]],da_data,RNAi_GeneName_to_exclude=da_info$GeneName,keep_x=F,keep_residuals=T,differential_analysis_column_name=da_info$colName)
        cv_to_use[[new_name]] = list(tmp[[1]],2^da_model_results$residuals+da_model_results$log2_intercepts)
        names( cv_to_use[[new_name]] ) = c(names(tmp)[1],paste0(names(tmp)[1],"_",da_model_name))
      }else{
        new_name = paste0(gn,"_",da_model_name)
        cat(paste0("Running DA model on ",gn," submodel ",names(tmp)[[1]],"...\n"))
        da_model_results_1 = ns_fit_da_linear_model(tmp[[1]],da_data,RNAi_GeneName_to_exclude=da_info$GeneName,keep_x=F,keep_residuals=T,differential_analysis_column_name=da_info$colName)
        cat(paste0("Running DA model on ",gn," submodel ",names(tmp)[[2]],"...\n"))
        da_model_results_2 = ns_fit_da_linear_model(tmp[[2]],da_data,RNAi_GeneName_to_exclude=da_info$GeneName,keep_x=F,keep_residuals=T,differential_analysis_column_name=da_info$colName)
        if (da_model_results_1$relative_variance < .95){
          d1 = 2^da_model_results_1$residuals+da_model_results_1$log2_intercepts;
        }else d1 = tmp[[1]]
        if (da_model_results_2$relative_variance < .95){
          d2 = 2^da_model_results_2$residuals+da_model_results_2$log2_intercepts;
        }else d2 = tmp[[2]]
        cv_to_use[[new_name]] = list(d1,d2)
        names( cv_to_use[[new_name]] ) = names(tmp)
      }
      
        cv_to_use[[gn]] = NULL;
    }
    model_names = names(cv_to_use)
  }
  
  #filter counts to genes expressed at a minimum level in both covariate groups
  
  min_mu = 50
  cat(paste0("Filtering genes lower than ",min_mu," counts...\n"))
  
  # select genes with high mean
  intersect_genes <-lapply(cv_to_use, function(x) Reduce(intersect,lapply(x,rownames)))
  names(intersect_genes) = names(cv_to_use)
  
  mu_list <- lapply(names(cv_to_use), function(group) lapply(names(cv_to_use[[group]]),function(subgroup){ 
    apply(cv_to_use[[group]][[subgroup]][intersect_genes[[group]], ], 1, mean)}))
  names(mu_list) = names(cv_to_use)
  
  select_genes <- lapply(mu_list,function(x){names(which(rowMeans(simplify2array(x)) > min_mu))})
  names(select_genes)= names(cv_to_use)
  cat("Gene Counts:\n")
  print(lapply(select_genes,length))
 
  filter_tissue_norm_list <- lapply(names(cv_to_use), function(group){
      ret = lapply(names(cv_to_use[[group]]),function(subgroup)cv_to_use[[group]][[subgroup]][select_genes[[group]], ])
      names(ret) = names(cv_to_use[[group]])
      ret;
      })
  names(filter_tissue_norm_list)= names(cv_to_use)

  #calculate correlation matricies
  cat("Calculating correlation matricies...\n")
  filter_tissue_rho_list <- lapply(filter_tissue_norm_list, function(x) lapply(x,function(y)cor(t(y), method = "s")))
  names(filter_tissue_rho_list)= names(filter_tissue_norm_list)
  
  filter_tissue_combined_rho_list <- lapply(filter_tissue_norm_list, function(x){ 
                                            y = do.call("cbind",x);
                                            cor(t(y), method = "s")})
  names(filter_tissue_combined_rho_list)= names(filter_tissue_norm_list)
  
  
  
  to_fit = names(filter_tissue_rho_list);
  community_info = lapply(to_fit,function(group){
    cat(paste("Processing",group,"\n")) 
    subgroup_names = names(filter_tissue_rho_list[[group]])
   
    if (length(subgroup_names)>2)
      stop(paste("Model",group, "has too many subgroups:",paste(subgroup_names,collapse=" ")))
    communities= do.call("c",lapply(subgroup_names,function(subgroup){
      cat(paste("Processing subgroup ",subgroup,"\n")) 
      cur_counts = filter_tissue_norm_list[[group]][[subgroup]]
  
      covar = filter_tissue_rho_list[[group]][[subgroup]]
      
      cat(paste("Fitting",subgroup,"\n"))
      flush.console()
      g = makeGraph(covar, theta = 0.3, beta = 2)
      comun <- cluster_leiden (g, objective_function ="modularity",weights = E(g)$weight,n_iterations=300,resolution_parameter=1)
      s = sizes(comun);
      s = s[s>5]
      if (length(s)==0) stop("All communities were too small")
      s = s[order(s,decreasing=T)]
      
     lapply(names(s),function(x){
         community_genes = comun$names[comun$membership == x] 
        
         sub_counts = cur_counts[community_genes,]
         pca_res = prcomp(t(sub_counts),scale=T,center=T)
         pca_info = list(scale=pca_res$scale,center=pca_res$center,rotation=pca_res$rotation,frac_var = pca_res$sdev^2/sum(pca_res$sdev^2));
         list(all_genes = community_genes, projection_info = pca_info,fit_subgroup = subgroup,group=group)
      })
    }))
    names(communities) = coexp_ID_from_number(1:length(communities))
    for (i in 1:length(communities)){
      communities[[i]]$within_group_community_id = i;
    }
    communities
  })
  names(community_info) = to_fit
  if (! (opt$compile_and_label && build_communities)){ #if we do everything in one run, we don't need to save to disk
   filename = paste0("./data/gene_variability/communities/community_effect_",output_suffix,file_suffix,".Rds")
   saveRDS(community_info,filename);
   filename = paste0("./data/gene_variability/communities/count_data_",output_suffix,file_suffix,".Rds")
   saveRDS(filter_tissue_norm_list,filename);
  }
}
  if (!opt$compile_and_label){
    cat("Finished.\n");
  }else{
    if (!(opt$compile_and_label && build_communities)){
      path = paste0("./data/gene_variability/communities/");
      fls = list.files(path=path,
                       pattern=paste0("community_effect_.*",file_suffix,".Rds"))
      if (length(fls) == 0)
        stop("Could not find any files to aggregate.")
      cat(paste0("Loading files:\n",paste(fls,collapse="\n")))
      cat ("\n")
      community_info = Reduce(append,lapply(fls,function(x){
        readRDS(paste0(path,x))
      }))
      
      fls = list.files(path=path,
                       pattern=paste0("count_data_.*",file_suffix,".Rds"))
      

      filter_tissue_norm_list = Reduce(append,lapply(fls,function(x){
        readRDS(paste0(path,x))
      }))
      if (any(!names(filter_tissue_norm_list)  %in% names(community_info)))
        stop("Incomplete count data found")
      if (any(! names(community_info) %in% names(filter_tissue_norm_list) ))
        stop("Extraneous count data found"
             )
      cat(paste("Loaded models",paste(names(community_info),collapse=", ")))
      cat ("\n")
    }
  
  #build data structure with lots of data on communities that we'll then use to merge strongly overlapping communities
  cat("Comparing communities...")
  all_communities = do.call("c",community_info)
  community_pca_correlation = matrix(NA,nrow=length(all_communities),ncol=length(all_communities))
  community_jaccard_index = matrix(NA,nrow=length(all_communities),ncol=length(all_communities));
  
  last_per = 0
  for (community_1 in 1:length(all_communities)){
    cur_per = floor(100*(community_1-1)/length(all_communities))
    if (cur_per > last_per+4){
      cat(paste(cur_per,"%.."))
      last_per = cur_per;
    }
    cat(".")
    for (community_2 in community_1:length(all_communities)){
     # if (community_1 == community_2) next;
      coms = list(all_communities[[community_1]],all_communities[[community_2]])
      com_names = list(paste(coms[[1]]$group,coms[[1]]$fit_subgroup),paste(coms[[2]]$group,coms[[2]]$fit_subgroup))
      
      community_source_counts = lapply(coms,function(x){filter_tissue_norm_list[[x$group]][[x$fit_subgroup]]})
      cur_genes = lapply(coms,function(x){x$all_genes})
      
      
      community_jaccard_index[community_1,community_2] = community_jaccard_index[community_2,community_1] = jaccard(cur_genes[[1]],cur_genes[[2]])
     
     
      cur_pca_info = lapply(coms,function(x){x$projection_info})
      
      
      #collect the genes found in common between both data sets
      cur_common_genes = list(cur_genes[[1]][cur_genes[[1]] %in% rownames(community_source_counts[[2]])],
                                cur_genes[[2]][cur_genes[[2]] %in% rownames(community_source_counts[[1]])])
      if (length(cur_common_genes[[1]]) < .95*length(cur_genes[[1]]) ||
          length(cur_common_genes[[2]]) < .95*length(cur_genes[[2]]) ){
        #print(paste("Not enough gene overlap between ",paste(com_names,collapse=" ")))
        
        next;
      }
      #get the PCA project matrix that works for the genes in common between both data sets.
      cur_pca_info_common_genes = lapply(1:2,function(cur_community){
        gene_matchup = match(cur_common_genes[[cur_community]],cur_genes[[cur_community]])
        pca_info = list(scale=cur_pca_info[[cur_community]]$scale[gene_matchup],center=cur_pca_info[[cur_community]]$center[gene_matchup],rotation=cur_pca_info[[cur_community]]$rotation[gene_matchup,]);
        pca_info
      })
    
      # stop()
      #calculate the PCA scores of each individual in each of the two data sets from which the communities were detected
      score_correlations= unlist(lapply(community_source_counts,function(sub_counts){ 
        scores = rbindlist(lapply(1:2,function(cur_community){
          sub_pca_info = cur_pca_info_common_genes[[cur_community]];
          sub_counts = t(scale(t(sub_counts[ cur_common_genes[[cur_community]],]),center=T,scale=T))
          scores = t(sub_counts) %*% sub_pca_info$rotation[,1]
          scores = data.frame(t(unlist(scores)))
          }))
        cor(t(scores[1,]),t(scores[2,]), method = "s")
        
      }))
    
     #score correlations is a vector with two numbers.
     #the first is the correlation of the two community scores on the source data for community 1
     #the second is the correlation of the two community scores on the source data for community 2
     community_pca_correlation[community_1,community_2] = community_pca_correlation[community_2,community_1] = max(score_correlations)
     
  
  
    }
  }
  
  cat("\nEstimating merge via hierarchical clustering...\n")
  nms = sapply(all_communities,function(x)paste(x$group,x$fit_subgroup,x$within_group_community_id))
  
  rownames(community_pca_correlation) = colnames(community_pca_correlation) = rownames(community_jaccard_index) = colnames(community_jaccard_index) =nms
  
  jaccard_index_cutoff = .30
  community_pca_correlation_cutoff= .9
  
 
  #MERGE TOGETHER USING two rounds of hierarchical clustering
  #The first uses both the jaccard index and the community PCA correlation.  
  #the second round just uses jaccard index (so we don't have to re-run the PCA analysis for the merged sets)
  
  combined_distance = community_jaccard_index;
  colnames(combined_distance) = names(all_communities)
  rownames(combined_distance) = names(all_communities)
  for (i in 1:length(all_communities)){
    for (j in i:length(all_communities)){
      if (i!=j && !is.na(community_pca_correlation[i,j])) combined_distance[i,j] = max(abs(community_pca_correlation[i,j]),community_jaccard_index[i,j]) 
       }   
  }
  merge_groups_community_correlation = rep(-1,length(all_communities));
  g_clust = as.dendrogram(hclust(as.dist(1-combined_distance),method="ward.D2"))
  
  similarity_cutoff = 1-seq(0,.95,by=.05)
  num_groups = sapply(similarity_cutoff,function(y)quantile(table(cutree(g_clust,h=1-y)),1))
  plot(similarity_cutoff,num_groups,xlab="similarity threshold",ylab="size of largest group")
  
  grps = cutree(g_clust,h=.75);
  grp_sizes = table(grps)[order(table(grps),decreasing=T)]
  
  groups = unique(grps);
  updated_jaccard = matrix(NA,nrow=length(groups),ncol=length(groups))
  for (i in 1:length(groups)){
    for (j in i:length(groups)){
      #if (i != j){
        coms_1 = names(grps)[which(grps == groups[i])]
        genes_1 = unique(do.call('c',lapply(coms_1,function(x){all_communities[[x]]$all_genes})))
        coms_2 = names(grps)[which(grps == groups[j])]
        genes_2 = unique(do.call('c',lapply(coms_2,function(x){all_communities[[x]]$all_genes})))
        updated_jaccard[i,j] = updated_jaccard[j,i] = jaccard(genes_1,genes_2)
     # }
    }
  }
  g_clust_2 = as.dendrogram(hclust(as.dist(1-updated_jaccard),method="ward.D2"))
  grps_2 = cutree(g_clust_2,h=.77);
  grp_sizes_2 = table(grps_2)[order(table(grps_2),decreasing=T)]
  
  merge_groups_community_correlation = rep(-1,length(all_communities));
  grp_name = 1;
  new_groups = unique(grps_2)
  for (i in 1:length(new_groups)){
    comsets = as.integer(names(grps_2)[which(grps_2 == new_groups[i])])
    coms = names(grps)[which(grps %in% comsets)]
    if (length(coms) <= 1) next;
    merge_groups_community_correlation[match(coms,names(all_communities))] = grp_name
    grp_name = grp_name + 1;
  }
  grp_sizes_3 = table(merge_groups_community_correlation)[order(table(merge_groups_community_correlation),decreasing=T)]
  cat(paste("Preparing to Merge",length(all_communities),"communities into",length(unique(merge_groups_community_correlation))-1+length(which(merge_groups_community_correlation==-1)),"communities\n"))
  
  #each value of merge groups tells the new group for the corresponding community in all_communties
  #eg. merge_group[i] is the group for all_communities[i]
  merge_groups = merge_groups_community_correlation
  merge_g= unique(merge_groups[merge_groups!=-1])
  
  #figure out the new names for the global naming scheme
  label_matchup_reference = data.table(new_coexp_group_ID = rep(NA,length(all_communities)),old_coexp_group_ID = rep(NA,length(all_communities)),group = rep(NA,length(all_communities)),fit_subgroup =  rep(NA,length(all_communities)), N =  rep(NA,length(all_communities)),merged_max_N =  rep(NA,length(all_communities)),group_subgroup_size = rep(NA,length(all_communities)),all_communities_index = rep(NA,length(all_communities)))
  
  #assign merged groups their IDS
  for (i in 1:length(merge_g)){
     merged_group_members = which(merge_groups==merge_g[i]);
     label_matchup_reference$new_coexp_group_ID[merged_group_members] = i;
     label_matchup_reference$old_coexp_group_ID[merged_group_members] = sapply(merged_group_members,function(x)all_communities[[x]]$within_group_community_id)
     label_matchup_reference$group[merged_group_members] =              sapply(merged_group_members,function(x)all_communities[[x]]$group)
     label_matchup_reference$fit_subgroup[merged_group_members] =       sapply(merged_group_members,function(x)all_communities[[x]]$fit_subgroup)
     label_matchup_reference$N[merged_group_members] =                  sapply(merged_group_members,function(x)length(all_communities[[x]]$all_genes))
     label_matchup_reference$merged_max_N[merged_group_members] =       max(label_matchup_reference$N[merged_group_members]);
     label_matchup_reference$group_subgroup_size[merged_group_members] = length(merged_group_members)
     label_matchup_reference$all_communities_index[merged_group_members] = merged_group_members
     
  }
  #longer_start_label is the ID for the first co-expression group that doesn't need to be merged with others
  
  loner_start_label = length(merge_g)+1
  loners = which(merge_groups == -1)
  if (length(loners)>0){
    label_matchup_reference$new_coexp_group_ID[loners] = loner_start_label-1 + 1:length(loners);
    label_matchup_reference$old_coexp_group_ID[loners] = sapply(loners,function(x)all_communities[[x]]$within_group_community_id)
    label_matchup_reference$group[loners] =             sapply(loners,function(x)all_communities[[x]]$group)
    label_matchup_reference$fit_subgroup[loners] =      sapply(loners,function(x)all_communities[[x]]$fit_subgroup)
    label_matchup_reference$N[loners] =                 sapply(loners,function(x)length(all_communities[[x]]$all_genes))
    label_matchup_reference$merged_max_N[loners] =      sapply(loners,function(x)length(all_communities[[x]]$all_genes))
    label_matchup_reference$group_subgroup_size[loners]  = 1;
    label_matchup_reference$all_communities_index[loners] = loners
  }
  
  #find the number of genes referred to by each new label
  all_new_labels = unique((label_matchup_reference[, .(new_coexp_group_ID,merged_max_N)]))
  all_new_labels = all_new_labels[order(all_new_labels$merged_max_N,decreasing=T),]
     
  new_ids_in_order = rep(NA,nrow(all_new_labels))
  #Reassign the new ids so they appear in descending order
  for (i in 1:nrow(all_new_labels)){
    new_ids_in_order[label_matchup_reference$new_coexp_group_ID ==all_new_labels$new_coexp_group_ID[i]] = i;
  }
  label_matchup_reference$new_coexp_group_ID = new_ids_in_order;
  label_matchup_reference = label_matchup_reference[order(label_matchup_reference$new_coexp_group_ID),]
  
  #provide a merged genelist
  all_new_labels = unique(label_matchup_reference[, .(new_coexp_group_ID,merged_max_N)])
  merged_genelists = list();
  for (i in all_new_labels$new_coexp_group_ID){
    group_members = label_matchup_reference[new_coexp_group_ID == i,all_communities_index]
    merged_genelists[[i]] = unique(do.call("c",lapply(group_members,function(x)all_communities[[x]]$all_genes)))
  }
 
  #Great!  Here we now have all the new co-expression group IDs and we can use them to merge the groups.
  
  merged_community_info = list();
  
  
  for (cur_group in unique(label_matchup_reference$group)){
    #cur_group = "n2_sw"
    print(paste("Merging coexpression groups in",cur_group));
    flush.console();
    cur_ref = label_matchup_reference[label_matchup_reference$group==cur_group,]
    new_coexpr_groups = unique(cur_ref$new_coexp_group_ID)
    
    
    
    merged_community_info[[cur_group]] = list();
    for (cur_coexp in new_coexpr_groups){
      #cat (paste0(cur_coexp,"..."))
      new_id = coexp_ID_from_number(cur_coexp);
      source_groups_to_merge = unlist(cur_ref[new_coexp_group_ID == cur_coexp,old_coexp_group_ID])
      if (length(source_groups_to_merge) == 1){
        #coexpr groups that don't need to be merged are easy!  We just copy them over using their new name
        merged_community_info[[cur_group]][[new_id]] = community_info[[cur_group]][[source_groups_to_merge]]
        source_gene_info = data.frame(GeneName = community_info[[cur_group]][[source_groups_to_merge]]$all_genes);
        source_gene_info$fit_subgroup = community_info[[cur_group]][[source_groups_to_merge]]$fit_subgroup
        source_gene_info$source_coexp_group_id = source_groups_to_merge
        merged_community_info[[cur_group]][[new_id]]$core_genes = merged_community_info[[cur_group]][[new_id]]$all_genes
       
        merged_community_info[[cur_group]][[new_id]]$source_gene_info = source_gene_info
        merged_community_info[[cur_group]][[new_id]]$global_group_community_id = new_id;
      }else{
        source_gene_info = rbindlist((lapply(source_groups_to_merge,function(x){
          res = data.frame(GeneName = community_info[[cur_group]][[x]]$all_genes);
          res$fit_subgroup = community_info[[cur_group]][[x]]$fit_subgroup
          res$source_coexp_group_id = x
          res
        })))
        # print("BWA")
     #   browser()
        source_gene_count = unlist(lapply(source_groups_to_merge,function(x)length(community_info[[cur_group]][[x]]$all_genes)))
        source_fit_subgroups =  unlist(lapply(source_groups_to_merge,function(x)community_info[[cur_group]][[x]]$fit_subgroup))
        #choose the fit_subgroup that has the most genes
        #we use this to recalculate the PCA later.
        best_fit_subgroup = unique(source_fit_subgroups[which(source_gene_count==max(source_gene_count))[1]])
        all_genes = unique(source_gene_info$GeneName)
        gene_frequency = table(source_gene_info$GeneName)
        core_genes = names(gene_frequency)[gene_frequency>median(gene_frequency)]
        
        cur_counts = filter_tissue_norm_list[[cur_group]][[best_fit_subgroup]]
        
        core_genes = core_genes[core_genes %in% rownames(cur_counts)]

        if (length(core_genes) == 0){ 
          core_genes = names(gene_frequency)[gene_frequency>=median(gene_frequency)]
          core_genes = core_genes[core_genes %in% rownames(cur_counts)]
        }
        if (length(core_genes) == 0) stop(paste0("Could not find any core genes for ",cur_group,new_id))
        sub_counts = cur_counts[core_genes,]
        pca_res = prcomp(t(sub_counts),scale=T,center=T)
        pca_info = list(scale=pca_res$scale,center=pca_res$center,rotation=pca_res$rotation,frac_var = pca_res$sdev^2/sum(pca_res$sdev^2));
        
  
         
        merged_community_info[[cur_group]][[new_id]] = 
          list(all_genes = all_genes, 
               core_genes = core_genes,
               projection_info = pca_info,
               fit_subgroup = best_fit_subgroup,
               group = cur_group,
               source_gene_info = source_gene_info,
               global_group_community_id = new_id);
       # browser()
      }
    }
  }
  
  generate_community_flatfile = function(com_info,gene_type){
    rbindlist(lapply(names(com_info),function(cur_group){
    #  print(cur_group)
      rbindlist(lapply(names(com_info[[cur_group]]),function(cur_com){
        cur_com_data = com_info[[cur_group]][[cur_com]]
        if (!gene_type %in% names(cur_com_data))
          stop(paste("In",cur_group,cur_com,"gene_type must be either all_genes or core_genes.  Found ",paste(names(cur_com_data),collapse=" ")))
        r = data.frame(GeneName = cur_com_data[[gene_type]])
        r$fit_subgroup = cur_com_data$fit_subgroup;
        r$fit_subgroup_2 = ifelse (length(cur_com_data$source_fit_subgroups)>1,setdiff(cur_com_data$source_fit_subgroups,cur_com_data$fit_subgroup),"")
        com_int = as.integer(substring(cur_com,2,4))
        r$commmunity = com_int
        r$group = cur_group
        r$coexp_group_ID = cur_com
        r$coexp_group_name = paste(cur_com_data$fit_subgroup,com_int)
        r = merge(r,tissues_df[,.(GeneName,Tissue)],by="GeneName",all.x=T)
        r
      }))
    }))
  }
  #add some extra (redundant) metadata to the communty_info structure so it 
  #has the exact same structure as the merged_community_info
  for (cur_group in names(community_info)){
   for (cur_com in names(community_info[[cur_group]])){
        source_gene_info = data.frame(GeneName = community_info[[cur_group]][[cur_com]]$all_genes);
        source_gene_info$fit_subgroup = community_info[[cur_group]][[cur_com]]$fit_subgroup
        source_gene_info$source_coexp_group_id = -1
        community_info[[cur_group]][[cur_com]]$source_gene_info = source_gene_info
        community_info[[cur_group]][[cur_com]]$core_genes = community_info[[cur_group]][[cur_com]]$all_genes
        community_info[[cur_group]][[cur_com]]$global_group_community_id = -1; 
   }
  }
  
  cat("Writing files...\n")
  saveRDS(community_info,paste0("./data/gene_network/model_specific_community_info",file_suffix,".Rds"))
  saveRDS(merged_community_info,paste0("./data/gene_network/global_community_info",file_suffix,".Rds"))
  saveRDS(filter_tissue_norm_list,paste0("./data/gene_network/community_info_count_data",file_suffix,".Rds"))
  
  
  community_membership = generate_community_flatfile(community_info,"all_genes");
  write.csv(community_membership,paste0("./data/gene_network/model_specific_community_membership_all_genes",file_suffix,".csv"),row.names=F)
  community_membership = generate_community_flatfile(community_info,"core_genes");
  write.csv(community_membership,paste0("./data/gene_network/model_specific_community_membership_core_genes",file_suffix,".csv"),row.names=F)

  merged_community_membership = generate_community_flatfile(merged_community_info,"all_genes");
  write.csv(merged_community_membership,paste0("./data/gene_network/global_community_membership_all_genes",file_suffix,".csv"),row.names=F)
  merged_community_membership = generate_community_flatfile(merged_community_info,"core_genes");
  write.csv(merged_community_membership,paste0("./data/gene_network/global_community_membership_core_genes",file_suffix,".csv"),row.names=F)
  
  cat("Re-calculating correlation matricies...\n")
  filter_tissue_rho_list <- lapply(filter_tissue_norm_list, function(x) lapply(x,function(y)cor(t(y), method = "s")))
  names(filter_tissue_rho_list)= names(filter_tissue_norm_list)
  
  
  gene_network = list(filter_tissue_norm_list = filter_tissue_norm_list,
       filter_tissue_rho_list = filter_tissue_rho_list,
       gene_annotations = gene_annotations,
       geneid = geneid,
       community_info=community_info,
       merged_community_info=merged_community_info,
       merged_community_membership=merged_community_membership,
       community_membership=community_membership
       )
  
  qsave(gene_network,paste0("./data/gene_network/community_networks",file_suffix,".qs"))
  
}
