suppressPackageStartupMessages({
  library(dplyr)
  library(reshape2)
  library(data.table)
  library(scran)
  library(optparse)
  library(pbapply)
})
source("./scripts/helpers/preprocess.R")
pbo <- pboptions(type="txt")

# handle arguments
option_list = list(
  make_option(c("-n", "--name"), 
              type="character",
              default=NULL, 
              help="Name of experiment."),
  
  make_option(c("-c", "--counts"), 
              type="character",
              default=NULL, 
              help="Merged counts CSV file."),
  
  make_option(c("-a", "--annotations"), 
              type="character",
              default=NULL, 
              help="Sample annotations CSV file."),
  
  make_option(c("-o", "--output"), 
              type="character",
              default=NULL,
              help="Output RDS file."),
  
  make_option(c("--tissue_normalization"), 
              default=FALSE,
              action="store_true",
              help="Should tissue-specific normalization be performed?"),
  
  make_option(c("--bootstrap"), 
              default=FALSE,
              action="store_true",
              help="Should sample bootstrap be performed?"))

if (1){
opt_parser <- OptionParser(option_list = option_list, prog = "Format Counts")
opt <- parse_args(opt_parser)
}else{
if (0){
opt <- list(name = "single_and_pooled_worms",
            counts = "./data/annotated_counts/counts_single_and_pooled_worms.csv.gz",
            annotations = "./data/annotations/sample_annotations_single_and_pooled_worms.csv",
            output = "./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds",
            tissue_normalization = T,
            bootstrap=F
            )
}else{
opt <- list(name = "net_validation",
            counts = "./data/annotated_counts/counts_net_validation.csv.gz",
            annotations = "./data/annotations/sample_annotations_net_validation.csv",
            output = "./data/formated_counts/counts_and_annots_net_validation.rds",
            tissue_normalization = T,
            bootstrap=F
            )
}
 
}


if (is.null(opt$name)) {
  stop("Missing name of experiment (--name)")
}

if (is.null(opt$counts)) {
  stop("Missing merged counts CSV file argument (--counts)")
}

if (is.null(opt$annotations)) {
  stop("Missing sample annotations CSV file argument (--annotations)")
}

if (is.null(opt$output)) {
  stop("Missing output RDS file argument (--output)")
}

# load data
cat("Reading data\n")

geneid <- fread("./data/annotations/c_elegans.PRJNA13758.WS265.geneIDs.txt.gz")
geneid[V3 == "", V3 := V4]
geneid <- geneid[, .(GeneName = V2, GeneSymbol = V3, GeneType  = V6)]
geneid <- geneid[! duplicated(GeneName)]
setkey(geneid, GeneName)

df <- fread(opt$counts)
annots <- fread(opt$annotations)
tb = table(annots$Sample)
if (any(tb>1)){
	print(tb[which(tb>1)])
	stop("Sample names are repeated in the annotation file.")
}

if ("Exclude" %in% colnames(annots)) {
  cat(paste0("Excluding ",nrow(annots[Exclude == TRUE]), " samples for which Exclude is marked as True\n"))
  annots <- annots[Exclude == FALSE,]
}
setkey(annots, Sample)

# format counts
cat("Formating counts and filtering low-quality samples\n")
allcounts <- acast(df, GeneName ~ Sample, value.var = "Count")
rm(df)

allcounts <- allcounts[, annots$Sample]
counts <- allcounts[grepl("^WB", rownames(allcounts)), ]
ercc <- allcounts[grepl("^ERCC", rownames(allcounts)), ]
rm(allcounts)

setkey(annots,Sample)
lenient_count_criteria_requested_for_sample = annots[colnames(counts),LenientCountCriteria]


# remove samples with low library size
sample_keep <- 
	(!lenient_count_criteria_requested_for_sample & colSums(counts) > 5e5 & 
	  colSums(ercc) > 100 & 
	  colSums(ercc) < 2e5) | 
  	(lenient_count_criteria_requested_for_sample & colSums(counts) > 5e4 & 
	  colSums(ercc) > 100 & 
	  colSums(ercc) < 1e4)
cat(paste("Removing",length(which(!sample_keep)),"samples",":"))
print(annots[!sample_keep,])

if (opt$name == "net_validation") {
  sample_keep <- sample_keep &
    colnames(counts) %in% annots[EvidenceForKnockdown == TRUE,]$Sample
}

counts <- counts[, sample_keep]
ercc <- ercc[, sample_keep]
rm(sample_keep)

annots <- annots[Sample %in% colnames(counts)]
setkey(annots, Sample)

counts <- counts[, annots$Sample]
ercc <- ercc[, annots$Sample]

# boostrap
if (opt$bootstrap) {
  cat("Performing bootstrap\n")
  sel_samples <- sample(ncol(counts), size = ncol(counts), replace = TRUE)
  counts <- counts[, sel_samples]
  ercc <- ercc[, sel_samples]
  annots <- annots[sel_samples, ]
  
  annots[, Sample := paste0(Sample, "_", seq_len(.N))]
  colnames(counts) <- annots$Sample
  colnames(ercc) <- annots$Sample
}

# format annotations
annots[, Directory := as.factor(Directory)]

if (opt$name == "single_and_pooled_worms" | opt$name == "single_worms_bootstrap") {
  annots[, Strain := relevel(as.factor(Strain), "QZ0")]
  annots[, Genotype := relevel(as.factor(Genotype), "N2")]
  annots[, Day := relevel(as.factor(Day), "8")]
  annots[, BiologicalReplicate := as.factor(BiologicalReplicate)]
  annots[, Temperature := as.factor(Temperature)]
  annots[, NumberOfWorms := as.factor(NumberOfWorms)]
  annots[, Food := as.factor(Food)]
  annots[, Group := factor(paste0(Strain, ".t",Temperature,".d", Day, ".f",Food,".br", BiologicalReplicate,
                                  if_else(NumberOfWorms == 1, ".sw", ".pw")))]
  
} else if (opt$name == "beads") {
  annots[, Group := Beads]
  
} else if (opt$name == "tso_only") {
  annots[, Group := TSOOnly]
  
} else if (opt$name == "glp1_glp4") {
  annots[, Genotype := as.factor(Genotype)]
  annots[, Group := Genotype]
  
} else if (opt$name == "tissue_specific_rnai") {
  annots[, RNAi := factor(RNAi, levels = c("EV", "ama-1", "rpb-2", "daf-2"))]
  annots[, Tissue := relevel(factor(Tissue), "All")]
  annots[, Group := factor(paste0(RNAi, ".", Tissue))]
  
} else if (opt$name == "net_validation") {
  
  annots[, Genotype := relevel(as.factor(Genotype), "N2")]
  annots[, RNAi := relevel(as.factor(RNAi), "EV")]
  annots[, Set := as.factor(Set)]
  annots[, Group := factor(paste0(Genotype, ".", RNAi, ".", Set))]
  
 } else {
  stop("Are you sure the first argument is correct?")
}

# filter genes based on counts
cat("Filtering genes with low counts\n")
counts <- filterGenes(counts[, annots$Sample], group=annots$Group)
ercc <- filterGenes(ercc[, annots$Sample],groups=NULL, min_count = 5, min_prop = .8)

# remove non protein coding genes
if (! opt$name %in% c("beads", "tso_only")) {
  counts <- counts[rownames(counts) %in% geneid[GeneType %in% "protein_coding_gene"]$GeneName, ]
}
geneid[, GeneType := NULL]

# normalize
cat("Performing standard normalization\n")
nf <- normalizationFactor(counts, groups = annots$Group)
nf_ercc <- normalizationFactor(ercc, groups = annots$Group)

# tissue specific normalization
if (opt$tissue_normalization) {
  cat("Performing tissue-specific normalization\n")
  tissues_df <- fread("./data/tissue_unique_genes/genes_unique_to_tissues.csv.gz")
  tissues <- pull(tissues_df, "Tissue", "GeneName")
  nf_tissues <- normalizationFactorTissueSpecific(counts = counts, tissues = tissues, groups = annots$Group)
}

## save
cat("Saving results\n")
out <- list(
  geneid = geneid, 
  annots = annots,
  counts = counts,
  ercc = ercc,
  nf = nf,
  nf_ercc = nf_ercc)

if (opt$tissue_normalization) {
  out <- c(out, list(nf_tissues = nf_tissues))
}

saveRDS(out, opt$output)

cat("Done\n")
