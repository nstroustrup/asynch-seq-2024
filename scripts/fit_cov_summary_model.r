
library(data.table)

library(reshape2)


source("../scripts/notebooks_helpers/gene_variability.R")
source("../scripts/helpers/enrichment_analysis.R")
source("../scripts/helpers/basics.R")
source("../scripts/helpers/clustering.R")
source("../scripts/helpers/preprocess.R")
source("../scripts/notebooks_helpers/gene_network.R")

alpha <- 1e-4 # significance level
genotype_colors <- c(N2 = "#990000ff", `daf-2(e1368)` = "#0b5394ff")
cluster_colors <- c("firebrick", "sienna1", "cornflowerblue", "darkorchid", "black")

## Load data

### Gene names

geneid <- fread("../data/annotations/c_elegans.PRJNA13758.WS265.geneIDs.txt.gz")
geneid[V3 == "", V3 := V4]
geneid <- geneid[, .(GeneName = V2, GeneSymbol = V3, GeneType  = V6)]
geneid <- geneid[! duplicated(GeneName)]
setkey(geneid, GeneName)

### Gene annotations

gene_annotations <- readRDS("../data/annotations/gene_annotations.rds")
gene_annotations <- gene_annotations[c("Wormcat", "Wormexp", "KEGG")]

### Set up model comparison structures to guide analysis


models = read.csv("../data/differential_analysis/models/single_and_pooled_worms.csv");

model_names = c("d8_glp1","d1_glp1","d8_UV","d1_UV","d_UV","d8_20v25C","d1_20v25C","n2_sw","daf2_sw","d1_sw","d8_sw","d_25C","n2_series_all")

#make a list of which samples match up to each subgroup (biological covariate group)

sw <- readRDS("../data/formated_counts/counts_and_annots_single_and_pooled_worms.rds")
  annotations = sw$annots;
  samples_list = lapply(model_names,function(group){
    print(group)
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
  bc <- readRDS("../data/batch_corrected_counts/single_and_pooled_worms.rds") 
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
tissues_df <- fread("../data/tissue_unique_genes/genes_unique_to_tissues.csv.gz")
tissues <- pull(tissues_df, "Tissue", "GeneName")
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

library(rms)
to_plot = names(counts_list)
to_plot = to_plot[!to_plot %in% c("n2_series_all","d8_HT115")]
cov_stats = lapply(names(counts_list),function(model){
   
    gp1_counts = counts_list[[model]] 
     print(paste(model,paste(dim(gp1_counts),collapse=",")))
    gp1_norm_counts = lapply(gp1_counts,function(x){
          nf = normalizationFactor(x)
          normalizeCounts(x,nf)
    })
    
    gp1_shuffled_norm_counts = lapply((gp1_counts),function(x){
          for (i in 1:nrow(x))
            x[i,] = x[i,sample(1:ncol(x),ncol(x))]  
          nf = normalizationFactor(x)
          normalizeCounts(x,nf)
    })
    names(gp1_shuffled_norm_counts) = paste0(names(gp1_shuffled_norm_counts),"_shuffled")
    gp1_norm_counts = c(gp1_norm_counts,gp1_shuffled_norm_counts)
    glp1_norm_counts = gp1_norm_counts
    
    if (0){
    names(gp1_norm_counts) = names(gp1_counts)
    glp1_norm_counts[["sim"]] = as.matrix(rbindlist(lapply(genes_to_use,function(g){
      dat = params_single_df[GeneName==g & Day==8,]
      data.frame(t(rnbinom(n=500,mu=2^dat$LogMean,size=1/2^(dat$LogDisp))))
    })))
    rownames(glp1_norm_counts[["sim"]]) = genes_to_use
    }
    
    min_mu = 20
    # select genes with high mean
    gn = function(x,y){print(paste(x,y));intersect(rownames(x),rownames(y))}
    intersect_genes = rownames(glp1_norm_counts[[1]])
    for (i in 1:length(glp1_norm_counts))
      intersect_genes = intersect(intersect_genes,rownames(glp1_norm_counts[[i]]))
    
    mu_list <- sapply(names(glp1_norm_counts),function(subgroup){ 
      apply(glp1_norm_counts[[subgroup]][intersect_genes, ], 1, mean)
    })
    colnames(mu_list) = names(glp1_norm_counts)
    
    select_genes <- rownames(mu_list)[apply(mu_list,1,function(x)mean(x) > min_mu)]
    
    glp1_filter_norm_list  = lapply(names(glp1_norm_counts),function(subgroup)glp1_norm_counts[[subgroup]][select_genes, ])
    names(glp1_filter_norm_list) = names(glp1_norm_counts)
    
    aging_da = gene_variability_da[model==model_to_use]
    genes_to_use = intersect(aging_da$GeneName,params_single_df$GeneName)
    genes_to_use_2 = intersect(genes_to_use,select_genes)
    
    all_dat = rbindlist(lapply(glp1_filter_norm_list,function(x){
      m = cor(t(x[genes_to_use_2,]), method = "s")
      #take a random 1% of all correlations, otherwise the calculation gets too slow.
      x=sample.int(n=nrow(m)*ncol(m),size=nrow(m)*ncol(m)/1000,replace=F)-1
      x_i = x%%nrow(m)+1
      x_j = floor(x/nrow(m))+1
      not_diago = x_i != x_j#remove diagonal
      x_i= x_i[not_diago] 
      x_j= x_j[not_diago] 
      data.frame(value=m[x_i,x_j],Var1=rownames(m)[x_i],Var2=colnames(m)[x_j])
    }))
      
    all_dat$change_1 = aging_da[as.character(cv_focused$Var1),var_effect,]
    all_dat$change_2 = aging_da[as.character(cv_focused$Var2),var_effect,]
    all_dat$grp = grp
    
    all_dat[,change_type := ifelse(change_1=="." & change_2==".","neither",ifelse(change_1!="." & change_2!=".","one","both"))]
    
    tp = all_dat[change_1=="." & change_2==".",][sample(1:nrow(all_dat),1000000),]
    grps = unique(all_dat$grp)
    colors = c("#F8766D","#00BA38","black","dark gray")
    names(colors) = grps;
    lty = c("solid","solid","solid","dash")
    lty = c(1,1,2,2)
    names(lty) = grps
    ggplot(tp,aes(x=value,color=grp,fill=grp))+geom_density(alpha=.1,adjust =1,lwd=1) + xlim(-1,1)+geom_vline(xintercept=0,lwd=1,col="dark gray",lty=2)+xlab("Corrleation between gene-pairs")+scale_color_manual(values=colors)+scale_linetype_manual(values=lty)+scale_fill_manual(values=colors)
    ggsave(paste0("../figures/04/cov_plots/",model,"correlation_effect.pdf"),width=5,height=3)
    ggplot(tp,aes(x=abs(value),color=grp,fill=grp,linetype=grp))+stat_ecdf(lwd=1) + xlim(0,1)+xlab("Correlation Magnitude\nbetween gene-pairs")+scale_color_manual(values=colors)+scale_fill_manual(values=colors)+scale_linetype_manual(values=lty)
    ggsave(paste0("../figures/04/cov_plots/",model,"_correlation_effect_inset.pdf"),width=5,height=3)
    res = aggregate(abs(value)~grp,tp,mean)
    names(res)[2] = "mean"
    res2 = aggregate(abs(value)~grp,tp,median)
    names(res2)[2] = "median"
    res3 = aggregate(abs(value)~grp,tp,sd)
    names(res3)[2] = "sd"
    res4 = aggregate(abs(value)~grp,tp,length)
    names(res4)[2] = "N"
    r = merge(res,res2,by="grp")
    r = merge(r,res3,by="grp")
    r = merge(r,res4,by="grp")
    r$model = model
   
    tp$grp_f = factor(tp$grp,levels = names(gp1_norm_counts))
    tp$abs_val = abs(tp$val)
    nl_res = bj(Surv(abs_val)~grp_f,tp[!is.na(tp$abs_val) & tp$abs_val!=0,])
    list(nl_res = nl_res, means = r)
})
saveRDS(cov_stats, "../data/gene_variability/intervention_cov_stats.rds")
