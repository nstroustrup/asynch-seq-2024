suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(scran)
  library(igraph)
  library(optparse)
})
run_default = F
setwd("/users/nstroustrup/nstroustrup/projects/asynch-seq-2022")
source("./scripts/helpers/preprocess.R")
source("./scripts/helpers/graphs.R")

# test args
 opt <- list(
   counts_and_annots = "./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds",
   batch_corrected_counts = NULL,
   tissue_corrected_counts =  "./data/tissue_specific_regression/rds/single_and_pooled_worms.rds",
   counts_name = "n2_d8_sw",
   condition1_name = "n2",
   condition2_name = NULL,
   condition1_subset = 'Day == 8 & NumberOfWorms == 1 & Strain == "QZ0" & Temperature==20 & Food=="NEC937"',
   condition2_subset = NULL,
   output_rds =  "./data/gene_network/all_genes_gene_network_d8.rds",
   output_csv =  "./data/gene_network/all_genes_communities_d8.csv.gz",
   output_covariance = "./data/gene_network/all_genes_communities_d8_compact_file.rds",
min_community_size = 20,
theta = .6,
beta = 2,
min_mu = 30,
trials_infomap = 200
 )

# handle args
option_list = list(
  
  make_option(c("-i", "--counts_and_annots"),
              type="character",
              default=NULL,
              help="Input DESeq2 object list RDS file."),
  
  make_option(c("--batch_corrected_counts"),
              type="character",
              default=NULL,
              help="Batch corrected counts RDS file."),
  
  make_option(c( "--tissue_corrected_counts"),
              type="character",
              default=NULL,
              help="Tissue and batch corrected counts RDS file."),
  
  make_option(c("--counts_name"),
              type="character",
              default=NULL,
              help="Name of counts matrix to use (in case of batch or tissue-corrected)."),
  
  make_option(c("--condition1_name"),
              type="character",
              default=NULL,
              help="Name of first condition."),
  
  make_option(c("--condition2_name"),
              type="character",
              default=NULL,
              help="Name of second condition."),
  
  make_option(c("--condition1_subset"),
              type="character",
              default=NULL,
              help="Name of first condition."),
  
  make_option(c("--condition2_subset"),
              type="character",
              default=NULL,
              help="Name of second condition."),
  
  make_option(c("--min_community_size"),
              type="integer",
              default=10L,
              help="Minimum community size."),
  
  make_option(c("--theta"),
              type="numeric",
              default=0.5,
              help="Correlation^beta threshold."),
  
  make_option(c("--beta"),
              type="numeric",
              default=2,
              help="Correlation power."),
  
  make_option(c("--min_mu"),
              type="numeric",
              default=30,
              help="Minimum count means."),
  
  make_option(c("--trials_infomap"),
              type="integer",
              default=200L,
              help="Clustering trials."),
  
  make_option(c("--output_csv"),
              type="character",
              default=NULL,
              help="Output CSV file."),
  
  make_option(c("--output_rds"),
              type="character",
              default=NULL,
              help="Output RDS file."),

  make_option(c("--output_covariance"),
	      type="character",
	      default=NULL,
	      help="output covariance matrix.")
)
if (!run_default){
   opt_parser <- OptionParser(option_list = option_list, prog = "Network Analysis")
   opt <- parse_args(opt_parser)
}
if (is.null(opt$counts_and_annots)) {
  stop("Formated counts RDS file (--counts_and_annots)")
}

# handle parameters
final_theta <- opt$theta
thetas <- c(final_theta)
names(thetas) <- paste0("Theta", thetas)
final_theta_name <- paste0("Theta", final_theta)

# load data
cat("Loading data\n")

data <- readRDS(opt$counts_and_annots)
for (var in names(data)) {
  assign(var, data[[var]])
}
rm(var, data)


tissues_df <- fread("./data/tissue_unique_genes/genes_unique_to_tissues.csv.gz")
tissues <- pull(tissues_df, "Tissue", "GeneName")



if(is.null(opt$batch_corrected_counts) & is.null(opt$tissue_corrected_counts)) {
  cat("Using raw counts.\n")
  x <- counts
  
} else if (! is.null(opt$batch_corrected_counts) & ! is.null(opt$tissue_corrected_counts)) {
  stop("You specified both batch_corrected_counts and tissue_corrected_counts. I don't know what to do. Aborting.")
} else {
  
  if (! is.null(opt$batch_corrected_counts)) {
    cat("Using batch-corrected counts.\n")
    
    data <- readRDS(opt$batch_corrected_counts)
    sel_counts <- data$batch_counts_list[[opt$counts_name]]
    rm(data)
    
  } else {
    
    cat("Using tissue-specific counts from regression.\n")
    data <- readRDS(opt$tissue_corrected_counts)
    if (!opt$counts_name %in% names(data))
    stop(paste("Counts name", opt$counts_name, " not found in tissue corrected counts"))    
sel_counts <- data[[opt$counts_name]]$counts
    sel_norm <- data[[opt$counts_name]]$norm
    rm(data)
    
    
  }
  
}

# counts list
cat("Formatting data\n")

counts_list <- list(
  condition1 = sel_counts[, annots[eval(parse(text = opt$condition1_subset))]$Sample]);
names(counts_list) = opt$condition1_name;
if(!is.null(opt$condition2_subset)){
  counts_list$condition2 = sel_counts[, annots[eval(parse(text = opt$condition2_subset))]$Sample]
names(counts_list) <- c(opt$condition1_name, opt$condition2_name)
}
rm(sel_counts)

if (is.null(opt$tissue_corrected_counts)) {
  tissue_norm_list <- lapply(counts_list, function(x) {
    k <- normalizationFactorTissueSpecific(x, tissues)
    normalizeCountsTissueSpecific(x, k, tissues)
  })
} else {
  tissue_norm_list <- list(
    condition1 = sel_norm[, annots[eval(parse(text = opt$condition1_subset))]$Sample])
    names(tissue_norm_list) = opt$condition1_name;
    if (!is.null(opt$condition2_subset)){
       tissue_norm_list$condition2 = sel_norm[, annots[eval(parse(text = opt$condition2_subset))]$Sample]
  
	names(tissue_norm_list) <- c(opt$condition1_name, opt$condition2_name)
  }
#  rm(sel_norm)
}

# select genes with high mean
if (!is.null(opt$condition2_subset)){
   intersect_genes <- Reduce(intersect, lapply(tissue_norm_list, rownames))
}else intersect_genes = rownames(tissue_norm_list[[1]])

mu_list <- lapply(tissue_norm_list, function(x) 
  apply(x[intersect_genes, ], 1, mean))

select_genes <- names(which(rowMeans(simplify2array(mu_list)) > opt$min_mu))
filter_tissue_norm_list <- lapply(tissue_norm_list, function(x) x[select_genes, ])

cat(paste0("We will look a total of ", length(select_genes), " genes\n"))

# compute correlations
cat("Computing correlations\n")
filter_tissue_rho_list <- lapply(filter_tissue_norm_list, function(x) cor(t(x), method = "s"))

# make graphs and union graphs
cat("Making individual graphs\n")
tissue_graph_list <- lapply(thetas, function(theta) {
  lapply(filter_tissue_rho_list, function(x) {
    makeGraph(x, theta = theta, beta = opt$beta)
  })
})
names(tissue_graph_list) = names(thetas)

cat("Making union graphs\n")
if (!is.null(opt$condition2_subset)){
tissue_union_graph_list <- lapply(tissue_graph_list, function(x) 
  makeUnionGraph(x[[1]], x[[2]]))
}else {
      tissue_union_graph_list = tissue_graph_list[[1]]
      names(tissue_union_graph_list) = names(tissue_graph_list)[1]
}
      
# detect communities
set.seed(63)
cat("Finding communities\n")
tissue_union_community_list <- lapply(tissue_union_graph_list, function(g) {
  # cluster_walktrap(g)
print(Sys.time()) 
 g = cluster_infomap(g, e.weights = E(g)$weight, nb.trials = opt$trials_infomap)
print(Sys.time())
g
})

names(tissue_union_community_list) = names(tissue_union_graph_list)
print("done with communities")
tissue_union_membership_list <- lapply(tissue_union_community_list, function(g) {
  m <- membership(g)
  m <- tapply(names(m), m, function(x) x)
  m <- m[sapply(m, length) >= opt$min_community_size]
  m <- m[order(- sapply(m, length))]
  names(m) <- paste0("Community", seq_along(m))
  # names(m) <- paste0("Community", 
  #                    stringr::str_pad(seq_along(m),
  #                                     width = ceiling(log10(length(m)+1)), 
  #                                     side = "left", 
  #                                     pad = "0"))
  m
})

# plot communities
# for (i in names(tissue_union_community_list)) {
#   for (j in names(tissue_union_community_list[[i]])) {
#     plot(tissue_union_community_list[[i]][[j]],
#          tissue_union_graph_list[[i]][[j]],
#          vertex.label = NA,
#          vertex.size = 0.5,
#          main = paste(i, j))
#   }
# }

# consensus clusters
tissue_graph <- tissue_graph_list[[final_theta_name]]
tissue_union_graph <- tissue_union_graph_list[[final_theta_name]]
tissue_communities <- tissue_union_membership_list[[final_theta_name]]
tissue_communities_vector <- unlist(lapply(seq_along(tissue_communities), function(i) {
  x <- tissue_communities[[i]]
  setNames(rep_len(i, length(x)), x)
}))

# tests
# source("./scripts/helpers/heatmap.R")
# library(ComplexHeatmap)
# plotSplitCorrelationHeatmap(filter_tissue_rho_list$n2, tissue_communities)
# net_val_annots <- fread("./data/annotations/sample_annotations_net_validation.csv")
# rnai <- unique(net_val_annots$RNAi)
# mean(rnai %in% geneid[unlist(tissue_communities)]$GeneSymbol)
# rnai_df <- net_val_annots[RNAi != "EV" & RNAi != "21U-RNA 6959"][! duplicated(RNAi), 
#                                                                  .(RNAiSymbol = RNAi, Set)]
# rnai_df <- rnai_df[order(Set, RNAiSymbol)]
# rnai_df[, RNAiName := geneid[RNAiSymbol, on = "GeneSymbol"]$GeneName]
# rnai_df[, TissueUnique := RNAiName %in% names(tissues)]
# rnai_df[, PassMeanFilter := RNAiName %in% rownames(filter_tissue_norm_list$n2)]
# rnai_df[, InCommunity := RNAiName %in% unlist(tissue_communities)]
# rnai_df[, .(mean(InCommunity)), by = .(Set)]
# fwrite(rnai_df, "~/Desktop/rnai.csv")

# print genes
for (community in names(tissue_communities)) {
  cat(community)
  cat("\n")
  cat(sort(geneid[tissue_communities[[community]]]$GeneSymbol))
  cat("\n\n")
}

# save results
cat("Saving results\n")

if (! is.null(opt$output_csv)) {
  comm_df <- rbindlist(lapply(names(tissue_communities), function(i)
    data.table(Community = i, GeneName = tissue_communities[[i]])))
  comm_df[, GeneSymbol := geneid[GeneName]$GeneSymbol]
  comm_df <- comm_df[, Tissue := tissues_df[GeneName, on = "GeneName"]$Tissue]
  comm_df <- comm_df[, .(Community, Tissue, GeneSymbol, GeneName)]
  
  fwrite(comm_df, opt$output_csv)
}

if (! is.null(opt$output_rds)) {
  out <- list(
    counts_list = counts_list,
    tissue_norm_list = tissue_norm_list,
    mu_list = mu_list,
    
    filter_tissue_norm_list = filter_tissue_norm_list,
    filter_tissue_rho_list = filter_tissue_rho_list,
    
    tissue_graph = tissue_graph,
    tissue_union_graph = tissue_union_graph,
    tissue_communities = tissue_communities)
  
  saveRDS(out, opt$output_rds)
}
if (!is.null(opt$output_covariance)){
   out = list(
    filter_tissue_rho_list = filter_tissue_rho_list,
    tissue_graph = tissue_graph,
    tissue_union_graph = tissue_union_graph,
    tissue_communities = tissue_communities
)
   saveRDS(out,opt$output_covariance)
}

cat("Done\n")
