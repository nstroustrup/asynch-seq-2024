suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
})

setwd("Z:\\nstroustrup\\projects\\asynch-seq-2022")

# Data comes from Coleen Murphy.
# https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007559
# https://worm.princeton.edu/media/download/all_tissue_prediction_scores.txt

# They didn't measure gene expression in germline. 
# We remove genes that are also highly expressed in the germline 
# and add a germline set using glp-1, glp-4 and tissue-specific ama-1 and rbp-2

# gene symbols
geneid <- fread("./data/annotations/c_elegans.PRJNA13758.WS265.geneIDs.txt.gz")
geneid[V3 == "", V3 := V4]
geneid <- geneid[, .(GeneName = V2, GeneSymbol = V3, GeneType  = V6)]
geneid <- geneid[! duplicated(GeneName)]
setkey(geneid, GeneName)


# load unique genes
soma_uniqgenes_df <- rbindlist(lapply(
  list.files("./data/tissue_unique_genes", pattern = "^unique"), 
  function(file) {
    d <- fread(paste0("./data/tissue_unique_genes//", file))[, .(GeneSymbol = Gene)]
    d[, Tissue := gsub("unique_to_|[.]csv", "", file)]
    d
  }
))
soma_uniqgenes_df <- merge(geneid[, .(GeneName, GeneSymbol)], 
                           soma_uniqgenes_df, by = "GeneSymbol")
setkey(soma_uniqgenes_df, Tissue)
soma_uniqgenes_df[, .N, by = .(Tissue)]

# check if duplicates
soma_uniqgenes_df[duplicated(GeneName)]

#calculate germline
{
	# load support vector machine weights
	svm_weights <- fread("./data/tissue_unique_genes/all_tissue_prediction_scores.txt")
	svm_weights <- svm_weights[, .(GeneName = geneid[symbol, on = "GeneSymbol"]$GeneName, 
				       germ_line, hypodermis, muscle = muscular_system, intestine, neurons = neuron)]
	svm_weights <- svm_weights[! is.na(GeneName)]

	# load glp-1, glp-4 and germline-specific RNAi
	glp1_glp4_da <- fread("./data/differential_analysis/tables/glp1_glp4/glp1_glp4.csv.gz")
	ama1_da <- fread("./data/differential_analysis/tables/tissue_specific_rnai/d8_Germline_ama1.csv.gz")
	rpb2_da <- fread("./data/differential_analysis/tables/tissue_specific_rnai/d8_Germline_rpb2.csv.gz")

	germline_df <- Reduce(function(x, y) merge(x, y, by = c("GeneName", "GeneSymbol"), all = TRUE),
			      list(glp1_glp4_da, ama1_da, rpb2_da))
	for (j in grep("log2FoldChange", colnames(germline_df), value = TRUE)) {
	  set(germline_df, i = which(is.na(germline_df[[j]])), j = j, value = 0)
	}
	for (j in grep("pvalue", colnames(germline_df), value = TRUE)) {
	  set(germline_df, i = which(is.na(germline_df[[j]])), j = j, value = 1)
	}
	for (j in grep("padj", colnames(germline_df), value = TRUE)) {
	  set(germline_df, i = which(is.na(germline_df[[j]])), j = j, value = 1)
	}

	# no evidence of being in germline
	log2fc <- 0
	alpha <- 0.01

	not_germline <- germline_df[(! (Genotype_glp1_vs_EV_log2FoldChange < log2fc & Genotype_glp1_vs_EV_padj < alpha)) &
				      (!  (Genotype_glp4_vs_EV_log2FoldChange < log2fc & Genotype_glp4_vs_EV_padj < alpha)) &
				      (! (RNAi_ama1_vs_EV_log2FoldChange < log2fc & RNAi_ama1_vs_EV_padj < alpha)) &
				      (! (RNAi_rpb2_vs_EV_log2FoldChange < log2fc & RNAi_rpb2_vs_EV_padj < alpha))]$GeneName
	length(not_germline)

	# strong expression in germline
	log2fc <- -1
	alpha <- 0.01
	theta <- 0

	germline <- intersect(germline_df[Genotype_glp1_vs_EV_log2FoldChange < log2fc & Genotype_glp1_vs_EV_padj < alpha & 
					    Genotype_glp4_vs_EV_log2FoldChange < log2fc & Genotype_glp4_vs_EV_padj < alpha &
					    RNAi_ama1_vs_EV_log2FoldChange < log2fc & RNAi_ama1_vs_EV_padj < alpha &
					    RNAi_rpb2_vs_EV_log2FoldChange < log2fc & RNAi_rpb2_vs_EV_padj < alpha]$GeneName,

			      svm_weights[hypodermis < theta & muscle < theta & intestine < theta & neurons < theta]$GeneName)
}



compute_ts_genes = function (soma_uniqgenes_df){
  # merge germline and soma
  uniqgenes_df <- rbind(soma_uniqgenes_df[GeneName %in% not_germline,], 
                        data.table(GeneSymbol = geneid[germline]$GeneSymbol,
                                   GeneName = germline,
                                   Tissue = "germ_line"))
  uniqgenes_df[, .N, by = .(Tissue)]
  
  additional_expert_annotation =  c("WBGene00006434", #also expressed in gut
              "WBGene00010131",   #also expressed in pharynx
              "WBGene00012463", #also expressed in neurons
              "WBGene00021114",#expressed in multiple tissues
              "WBGene00016977",#gonad and neurons
              "WBGene00011128",#intestine and neurons
              "WBGene00010130",#excretory cell and intestine
              "WBGene00006918",#excretory cell and intestine
              "WBGene00013025",#excretory cell and intestine
              "WBGene00010911",#also in neurons
              "WBGene00010756",#also in body wall muscle
              "WBGene00013127",#also in germline
              "WBGene00004930",#also coelomocytes
              "WBGene00011172",#and neuron
              "WBGene00013926",#and neuron
              "WBGene00002130",#also in germline
              "WBGene00001130",#also intestine
              "WBGene00001480",#also in neurons
              "WBGene00002106", #also in germline
              "WBGene00006919", #also in neurons,
              "WBGene00007924", #also distal tip cell
              "WBGene00009349", #also in neurons and embryos
              "WBGene00009724", #also in neurons
              "WBGene00010681", #also in muscle
              "WBGene00010759", #also pharynx
              "WBGene00012107", #also germline precursor cell
              "WBGene00012592", #also germline
              "WBGene00012592", #also germline
              "WBGene00013008", #also in neurons
              "WBGene00000176",#also in germline
              "WBGene00012673", #in neurons not germline
              "WBGene00007966" #in neurons not germline
  );
  uniqgenes_df = uniqgenes_df[! GeneName %in% additional_expert_annotation,]
  uniqgenes_df
}

uniqgenes_df = compute_ts_genes(soma_uniqgenes_df);

fwrite(uniqgenes_df, "./data/tissue_unique_genes/genes_unique_to_tissues.csv.gz", 
       row.names = FALSE, quote = FALSE)

