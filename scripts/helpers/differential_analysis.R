runDESeq <- function(counts, annots, nf, formula, nb_test = TRUE,
                     fitType = "parametric") {
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                        colData   = annots,
                                        design    = formula)
  
  # run differential analysis
  DESeq2::sizeFactors(dds) <- nf
  dds <- DESeq2::estimateDispersions(dds, fitType = fitType) # estimate overdispersion parameter
  if (nb_test) dds <- DESeq2::nbinomWaldTest(dds) # run NB GLM and do Wald Test on coefficients
  dds
}

formatDESeq <- function(dds, geneid = NULL,...) {
  
  if (is.null(dds)) return(NULL)
  
  da <- lapply(DESeq2::resultsNames(dds)[-1], function(i) {
    r <- DESeq2::results(dds, name = i,...)
    r <- as.data.table(r)
    r <- r[, .(log2FoldChange, lfcSE, pvalue, padj)]
    r[is.na(padj), padj := 1]
    colnames(r) <- paste0(i, "_", colnames(r))
    r
  })
  da <- Reduce(cbind, da)
  da[, GeneName := rownames(dds)]
  if (! is.null(geneid)) {
    da <- merge(geneid, da, by = "GeneName", all.x = FALSE, all.y = TRUE)
  }
  da
}
#sets the level order of the apropriate columns in annots such that subsequent analysis will 
#use the ControlGroup columns specified in model as the reference group
 set_reference_groups_from_model = function(model,annots){
      reference_groups = c(model$ControlGroup1,model$ControlGroup2)
      reference_groups = reference_groups[!is.na(reference_groups)]
      found_col_for_ref_group = rep(F,length(reference_groups))
      
       
    
      if (any(nchar(reference_groups)>0)){
  	   
  	    #find all the columns that appear as covariates in the model
  	    form = as.formula(paste("~", model$BiologicalCovariates))
	    covariate_columns <- labels(terms(form))
	    #remove cross terms
	    covariate_columns = covariate_columns[!grepl("\\*",covariate_columns) & !grepl("\\:",covariate_columns)]
  
  	    #go through each covariate column and check to see if it contains one of the reference groups
  	    for (col in covariate_columns){
  		m = which(reference_groups%in%unlist(unique(as.character(annots[,get(col)]))))
  		if (length(m) == 0) next;
  		if (length(m) > 1) stop(paste("Column",col,"has values that match both control groups",reference_groups[1],reference_groups[2]))
  		if (found_col_for_ref_group[m] == T)
  			stop("Control group",reference_groups[m]," matches multiple columns")
  		found_col_for_ref_group[m] =T;
  		#set the reference group!
  		print(paste("Setting",reference_groups[m], "as the control group for covariate",col))
  		annots[[col]] = factor(annots[,get(col)]);
  		annots[[col]]= relevel(annots[,get(col)],ref=reference_groups[m])
	    }
	    missing = which(found_col_for_ref_group==F)
	    if (length(missing)>0){
	    	browser()
	    	print(paste("Could not find columns for reference group(s)",paste(reference_groups[missing],collapse=" ")))
	    	stop()
	    }
      }else print("No reference groups specified")
      annots
  }

get_reference_groups_and_columns_from_model = function(model,annots){
      reference_groups = c(model$ControlGroup1,model$ControlGroup2)
      reference_groups = reference_groups[!is.na(reference_groups)]
      found_col_for_ref_group = rep(F,length(reference_groups))
      
       
      cols = list();
      i = 1;
      if (any(nchar(reference_groups)>0)){
  	   
  	    #find all the columns that appear as covariates in the model
  	    form = as.formula(paste("~", model$BiologicalCovariates))
	    covariate_columns <- labels(terms(form))
	    #remove cross terms
	    covariate_columns = covariate_columns[!grepl("\\*",covariate_columns) & !grepl("\\:",covariate_columns)]
  
  	    #go through each covariate column and check to see if it contains one of the reference groups
  	    for (col in covariate_columns){
  		m = which(reference_groups%in%unlist(unique(as.character(annots[,get(col)]))))
  		if (length(m) == 0) next;
  		if (length(m) > 1) stop(paste("Column",col,"has values that match both control groups",reference_groups[1],reference_groups[2]))
  		if (found_col_for_ref_group[m] == T)
  			stop("Control group",reference_groups[m]," matches multiple columns")
  		found_col_for_ref_group[m] =T;
  		cols[[i]] = c(col,reference_groups[m]);
  		names(cols[[i]]) = c("column","group")
	    }
	    missing = which(found_col_for_ref_group==F)
	    if (length(missing)>0){
	    	browser()
	    	print(paste("Could not find columns for reference group(s)",paste(reference_groups[missing],collapse=" ")))
	    	stop()
	    }
      }else print("No reference groups specified")
      cols
  }
