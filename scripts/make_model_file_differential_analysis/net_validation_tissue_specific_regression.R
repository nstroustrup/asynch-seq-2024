library(data.table)
library(stringr)

data <- readRDS("./data/formated_counts/counts_and_annots_net_validation.rds")
annots <- unique(data$annots[RNAi != "EV" & EvidenceForKnockdown == TRUE,
                             .(RNAi, Genotype, Set,Day)])

df <- rbindlist(lapply(seq_len(nrow(annots)), function(i) {
  a <- annots[i]
  
  name <- paste0(str_remove_all(a$Genotype, "\\(|\\)|-"), "_", str_remove(a$RNAi, "-"), "_Set", a$Set,"_Day",a$Day)
  subset <- paste0("Set == \'", a$Set, "' & Genotype == \'", a$Genotype, "\' & RNAi %in% c(\'EV\', \'", a$RNAi, "\') & Day == ",a$Day)
  
  res = data.table(Name = paste0(c("", "ts_"), name), 
             Subset = subset,
             BiologicalCovariates = "RNAi", 
             TechnicalCovariates = "NULL",
             BatchCorrectionType="standard",
             LenientCountCriteria=F,
             ControlGroup1="EV",
             use_optional_strain_prefix=F,
             TissueSpecificNormalization = c(FALSE, TRUE),
             Run = TRUE, 
             remove_RNAi_target=T,
             Set = a$Set, Day = a$Day, RNAi = a$RNAi, Genotype = a$Genotype, Tissue = a$Tissue,
	     WaldTest=T)
}))
df <- df[order(Set, Genotype, RNAi)]
df


fwrite(df, "./data/tissue_specific_regression/models/net_validation.csv", quote = TRUE)
