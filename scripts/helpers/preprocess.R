#Removes genes that are not present at the specified counts and proporiton in any group
filterGenes <- function(counts, groups, 
				min_count = 5, 
				min_prop_of_samples= 0.5) {
  
  if (is.factor(groups))
    groups <- droplevels(groups)
  
  if (is.null(groups)) {
     filter <- rowMeans(counts >= min_count) >= min_prop_of_samples
    out <- counts[filter, ]
  } else {
    filter <- tapply(colnames(counts), groups, function(i) {
      rowSums(counts[, i, drop = FALSE] >= min_count) >= min_prop_of_samples
    })
    filter <- simplify2array(filter)
    out <- counts[apply(filter, 1, any), ]
  }
  
  out
}

normalizationFactor <- function(counts, groups = NULL) {
   # XXX Major bug here!  Previously the drop levels occured after the group size check
   # which would then discard all group info.
   if (is.factor(groups))
    groups <- droplevels(groups)
    
  if (any(table(groups) <= 1)) {
    warnings("At least one group only has one sample. Not grouping to normalize.\n")
    groups <- NULL
  }
  
  
  
  if (is.null(groups)) {
    # nf <- DESeq2::estimateSizeFactorsForMatrix(counts)
    # nf <- scuttle::medianSizeFactors(counts)
    nf <- scran::calculateSumFactors(counts, min.mean = 3) 
    #s = colSums(counts);
    #nf = s/sum(s);
    warning("No groups provided to normalization factor!")
  } else {
    # nf <- DESeq2::estimateSizeFactorsForMatrix(counts)
    # nf <- scuttle::medianSizeFactors(counts, clusters = groups)
    nf <- scran::calculateSumFactors(counts, clusters = groups, min.mean = 3)
    #s = colSums(counts);
    #nf = s/sum(s);
  }
  names(nf) <- colnames(counts)
  nf
}

normalizationFactorTissueSpecific <- function(counts, tissues, groups = NULL) {
  uniq_tissues <- sort(unique(tissues))
  nf <- lapply(uniq_tissues, function(tissue) {
  
    genes <- names(tissues)[tissues == tissue]
    tissue_counts <- counts[rownames(counts) %in% genes, ]
    normalizationFactor(tissue_counts, groups = groups)
  })
  names(nf) <- uniq_tissues
  nf
}
normalizationFactorTissueSpecific_list <- function(counts,  groups = NULL) {
  nf <- lapply(names(counts), function(tissue) {
    tissue_counts <- counts[[tissue]];
    normalizationFactor(tissue_counts, groups = groups)
  })
  names(nf) <- names(counts)
  nf
}

normalizeCounts <- function(counts, nf) {
  sweep(counts, 2, nf, "/")
}

normalizeCountsTissueSpecific <- function(counts, nf_list, tissues) {
  uniq_tissues <- names(nf_list)
  out <- lapply(uniq_tissues, function(tissue) {
    genes <- names(tissues)[tissues == tissue]
    tissue_counts <- counts[rownames(counts) %in% genes, ]
    normalizeCounts(counts = tissue_counts, nf = nf_list[[tissue]])
  })
  out <- Reduce(rbind, out)
  out
}

find_covariate_levels = function(annots,control_group){
  
  #rename resp_disp values to "group1" or "group2"
  strain_levels = as.character(unique(annots$Strain));
  day_levels = as.character(unique(annots$Day));
  temperature_levels = as.character(unique(annots$Temperature));
  food_levels = as.character(unique(annots$Food));
  
  comp_levels = NULL;
  if (length(strain_levels) == 2){
    comp_levels = strain_levels;
    col_prefix = "Strain"; 
  }
  if (length(day_levels) == 2){
    comp_levels = day_levels;
    col_prefix = "Day"
  }
  if (length(temperature_levels) == 2){
    comp_levels = temperature_levels;
    col_prefix = "Temperature"
  }
  if (length(food_levels) == 2){
    comp_levels = food_levels;
    col_prefix = "Food"
  }
  if (is.null(comp_levels))
    return(NULL);
  
  control = which(comp_levels == control_group)
  if (control != 1)
    comp_levels = rev(comp_levels);
  return(list(covariate_column_name = col_prefix,levels=comp_levels))
}
