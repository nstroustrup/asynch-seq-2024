library(DESeq2)
library(dplyr)
library(scran)
library(data.table)
library(pbapply)
library(Hmisc)
pbo <- pboptions(type="txt")

source("scripts/helpers/preprocess.R")
source("scripts/helpers/differential_analysis.R")
source("scripts/helpers/gene_variability.R")

nboots <- 50

data <- readRDS("./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds")
bc <- readRDS("./data/batch_corrected_counts/single_and_pooled_worms.rds")

bio_counts <- bc$batch_counts_list$d8_aw

annots <- data$annots[Sample %in% colnames(bio_counts)]
setkey(annots, Sample)
annots[, Strain := droplevels(Strain)]
annots[, Genotype := droplevels(Genotype)]
annots[, Day := droplevels(Day)]
annots[, Group := droplevels(Group)]

bio_counts <- bio_counts[, annots$Sample]
bio_counts <- filterGenes(bio_counts, min_count = 5, min_prop = .5)

tech_counts <- data$ercc[, annots$Sample]
tech_counts <- filterGenes(tech_counts, min_count = 5, min_prop = .5)
tech_counts <- tech_counts[! rownames(tech_counts) %in% c("ERCC-00113", "ERCC-00136"), ]

rm(bc, data)
gc()

sample_list <- lapply(levels(annots$Group), function(group) 
  annots[group, on = "Group"]$Sample)
names(sample_list) <- levels(annots$Group)
groups <- annots$Group
grouped_annots <- annots[, .N, by = .(Group, Strain, Day, Pooled, BiologicalReplicate)]

# bootstraping estimates of NB and residual over-dispersion
resdisp <- pblapply(seq_len(nboots), function(boot_index) {
  
  # select samples
  boot_sample_list <- lapply(sample_list, function(samples) 
    sample(samples, length(samples), replace = TRUE))
  boot_samples <- unlist(unlist(boot_sample_list), use.names = FALSE)
  boot_annots <- copy(annots[boot_samples, on = "Sample"])
  boot_bio_counts <- bio_counts[, boot_samples]
  boot_tech_counts <- tech_counts[, boot_samples]
  rm(boot_samples)
  
  # rename samples so they are unique
  boot_annots[, Sample := paste0("BootSample", stringr::str_pad(1:.N, width = 5, side = "left", pad = "0"))]
  colnames(boot_bio_counts) <- boot_annots$Sample
  colnames(boot_tech_counts) <- boot_annots$Sample
  boot_groups <- boot_annots$Group
  boot_sample_list <- lapply(levels(boot_groups), function(group) 
    boot_annots[group, on = "Group"]$Sample)
  names(boot_sample_list) <- levels(boot_groups)
  
  # estimate negative binomial and residual over-dispersion
  out <- estimateResDisp(bio_counts = boot_bio_counts, 
                         tech_counts = boot_tech_counts, 
                         annots = boot_annots, 
                         sample_list = boot_sample_list, 
                         groups = boot_groups)
  
  if (! is.null(out$params)) {
    out$params[, BootIndex := boot_index]
  }
  
  if (! is.null(out$trend)) {
    out$trend[, BootIndex := boot_index]
  }
  
  out
})

# format results into data.frame
trend <- rbindlist(lapply(resdisp, `[[`, "trend"))
trend <- merge(trend, grouped_annots, by = "Group")
setkey(trend, Group)
fwrite(trend, "./data/gene_variability/trend_single_vs_pooled.csv.gz")

params <- rbindlist(lapply(resdisp, `[[`, "params"))
params <- merge(params, grouped_annots, by = "Group")
setkey(params, GeneName, BootIndex, Group)
fwrite(params, "./data/gene_variability/residual_overdispersion_single_vs_pooled.csv.gz")

# summarize trend results
trend_df <- summarizeTrend(trend)
fwrite(trend_df, "./data/gene_variability/summary_trend_single_vs_pooled.csv.gz")

# summarize bootstrap results and compute HVG p-values
params_df <- summarizeResDisp(params)
fwrite(params_df, "./data/gene_variability/summary_residual_overdispersion_single_vs_pooled.csv.gz")
