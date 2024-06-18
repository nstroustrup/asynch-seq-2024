suppressPackageStartupMessages({
  library(data.table)
  library(pbapply)
})
pbo <- pboptions(type="txt")

# handle arguments
args <- commandArgs(trailingOnly=TRUE)
#args = "single_and_pooled_worms"
if (length(args) == 0) {
  stop("Please supply one argument")
}
name <- args[1]

# load annotations
annots <- fread(paste0("./data/annotations/sample_annotations_", name, ".csv"))
tb = table(annots$Sample)
if (any(tb>1)){
	print(tb[which(tb>1)])
	stop("Sample names are repeated in the annotation file.")
}

#first, if the user has requested that we remove the targets of RNAi from the counts file, we check that all the annotations are matched up.
if ("remove_RNAi_target" %in% colnames(annots)){

	cat("RNAi target removal has been requested in the sample annotations file. Confirming that all annotations can be disambiguated...\n")
	geneid <- fread("./data/annotations/c_elegans.PRJNA13758.WS265.geneIDs.txt.gz")
	geneid[V3 == "", V3 := V4]
	geneid <- geneid[, .(GeneName = V2, GeneSymbol = V3, GeneType = V6)]
	geneid <- geneid[! duplicated(GeneName)]
	setkey(geneid, GeneSymbol)

	bad_symbols = c();
	for (i in which(annots$remove_RNAi_target)){
		RNAi_GeneName = geneid[annots$RNAi[i],GeneName]
		if (is.na(RNAi_GeneName) | is.null(RNAi_GeneName)){
			if (annots$RNAi[i] %in% bad_symbols)
				next;
			print(paste("Could not find the GeneSymbol corresponding to RNAi ",annots$RNAi[i], " in the geneIDs table."))
			bad_symbols = c(annots$RNAi[i],bad_symbols)
		}
	}
	if (length(bad_symbols) > 0)
		stop("Un-recognizable RNAi names found in the sample annotations file. Please fix the names or set the \"remove_RNAi_target\" column value to false for these samples");
			
	print("Done.")
}

# read counts
cat("Reading counts\n")
df <- pblapply(1:nrow(annots), function(i) {
  a <- annots[i]
  fn <- paste0("../../microscopy/sequence_data/", 
               a$Directory, 
               "/data/counts/", 
               a$Basename,
               ".count.txt.gz")
  
  x <- tryCatch(fread(fn, 
                      header = TRUE, drop = 2:6,
                      colClasses = c("character", "character", "character", "character", "character", "integer", "integer")), 
                error = function(e) {
                  print(fn)
                  return(data.table())
                })
  
  if (nrow(x) == 0) {
    print(paste0(fn, " had zero rows. Is it missing?"))
    return(NULL)
  }
  
  colnames(x)[1] <- "GeneName"
  colnames(x)[2] <- "Count"
  x[, Basename := a$Basename]
  x[, Directory := a$Directory]
  
  
  if("remove_RNAi_target" %in% names(a))
  	if (a$remove_RNAi_target){
	RNAi_GeneName = geneid[a$RNAi,GeneName]
	wh = which(x$GeneName == RNAi_GeneName)
  	cat(paste("Setting",length(wh), "transcript(s)", RNAi_GeneName, "count to zero in",a$RNAi,"RNAi data\n"))
	x[wh,Count:= 0];	
   }
  x
})
df <- rbindlist(df)

# merge with sample names
cat("Giving sample names\n")
df <- merge(annots[, .(Sample, Basename, Directory)], df, by = c("Basename", "Directory"))

# remove zero counts
cat("Removing zero counts\n")
allzeros <- df[, .(mean(Count == 0)), by = .(GeneName)][V1 == 1]$GeneName
filtered_df <- df[! GeneName %in% allzeros, .(Sample, GeneName, Count)]
setkey(filtered_df, Sample, GeneName)

# write results
cat("Writing results\n")
fwrite(filtered_df, paste0("./data/annotated_counts/counts_", name, ".csv.gz"), quote = FALSE, row.names = FALSE)

cat("Done\n")