library(data.table)

data <- readRDS("../../../scripts/WormAnnotations/data/rds/gene_annotations.rds")
for (var in names(data)) {
  assign(var, data[[var]])
}
rm(data)

raw_annotations_list <- list(
  `Biological Process` = go_df[Ontology == "P", .(GeneName, Annotation = Name)],
  # `Cellular Compartment` = go_df[Ontology == "C", .(GeneName, Annotation = Name)],
  # `Molecular Function` = go_df[Ontology == "F", .(GeneName, Annotation = Name)],
  # `Tissue` = tissue_df[, .(GeneName, Annotation = Tissue)],
  # `Phenotype` = pheno_df[, .(GeneName, Annotation = Phenotype)],
  `Wormcat` = wormcat_df[, .(GeneName, Annotation = Category3)],
  `Wormexp` = wormexp_df[, .(GeneName, Annotation)],
  `KEGG` = kegg_df[, .(GeneName, Annotation)],
  # `WormexpTarget` = wormexp_df[Category == "Targets", .(GeneName, Annotation)],
  TF = tf_df[, .(GeneName, Annotation = TFSymbol)]
)

annotations_list <- lapply(raw_annotations_list, function(a) {
  keep <- a[, .N, by = .(Annotation)][N >= 5]$Annotation
  a <- a[Annotation %in% keep]
  setkey(a, GeneName)
  a
})

saveRDS(annotations_list, "./data/annotations/gene_annotations.rds")
