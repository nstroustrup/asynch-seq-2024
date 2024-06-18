library(data.table)
library(stringr)

data <- readRDS("./data/formated_counts/counts_and_annots_tissue_specific_rnai.rds")
annots <- unique(data$annots[RNAi != "EV", .(Day, RNAi, Tissue)])

df <- rbindlist(lapply(seq_len(nrow(annots)), function(i) {
  a <- annots[i]
  
  name <- paste0("d", a$Day, "_",
                 str_replace_all(a$Tissue, " ", "_"), "_", 
                 str_remove(a$RNAi, "-"))
  subset <- paste0("Day == \'", a$Day, "' & Tissue == \'", a$Tissue, "\' & RNAi %in% c(\'EV\', \'", a$RNAi, "\')")
  
  data.table(Name = name, Subset = subset,
             BiologicalCovariates = "RNAi", TechnicalCovariates = "NULL",
             TissueSpecificNormalization = c(FALSE),
             Run = TRUE, 
             Day = a$Day, RNAi = a$RNAi, Tissue = a$Tissue)
}))
df

fwrite(df, "./data/differential_analysis/models/tissue_specific_rnai.csv", quote = TRUE)