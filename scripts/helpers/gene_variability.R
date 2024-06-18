#find the overdispersion and mean for the mRNAs (bio_counts) and spike-ins (tech_counts)
#and then generate a best estimate of the "true" dispersion for each mRNA, unconfounded
#by the mean-variance relationship of the underlying sequencing technology,
#by estimating the difference between each mRNA's dispersion and the 
#dispersion measured for the spike-ins at that abundance. 
#This difference is calculated using a SSasymp model fit to the spike-ins
#gene_list is an optional parameter that limits analysis to a specific set of genes
#library("gamlss")
library("optimization")

#Given a data frame with LogMean and LogDisp parameter estimates and a SequenceType column specifying spike-ins or mRNA
#calculate the local minimum LogDisp measured around each count specified in new_x below
#This is done by binning the biological counts and spike-in counts separately and finding minimim of all the bins
#The local minimum can then be used to calculate residual overdispersion.
fit_bio_and_tech_running_minimum = function(good_params){
	new_x = c(1,1.5,2,2.5,3,3.5,4,5,6,7,10,25,50,75,100,250,500,750,1000,2500,5000,7500,10000,25000,50000,75000,100000,250000,500000)

	x_tech = good_params[SequenceType == "Technical",]
	tech_min = min(x_tech$LogMean)
	x_bio = good_params[good_params$LogMean <= tech_min,];

	x_tech_range_all = good_params[good_params$LogMean >= tech_min]

	run_min = function(step_range,params_for_binning_calc,params_for_min_calc){
	   xs = params_for_binning_calc$LogMean
	   xs = xs[order(xs,decreasing=T)]

	  x_steps = step_range; #find step size that leaves least left out after last bin.
	  x_left = length(xs)%%x_steps
	  x_num = floor(length(xs)/x_steps)
	  x_left = x_left[x_num>=1] #we need one left over for the index +1 
	  x_step = min(x_steps[which(x_left == min(x_left))])

	    x_i = seq(1,length(xs),by=x_step);
	    if (length(x_i) <3){
	    	warning(paste("When fitting model for ResVar model, the specified range for datapoints per step,",paste(range(step_range),collapse=":"),
	    			",yielded only", length(x_i)  , " subdivisions(s)"))
	    	min_cutoff = quantile(params_for_min_calc[,LogVar],.05)
	        return(params_for_min_calc[LogMean <= xb[1] & LogMean > xb[2] & (LogVar <= min_cutoff),.(LogMean,LogVar,SequenceType)])
	 
	    }
	    x_mins =rbindlist(lapply(2:(length(x_i)-1),FUN=function(i){
	   	      xb = xs[c(x_i[i-1],x_i[i+1])];
	   	      min_cutoff = min(params_for_min_calc[LogMean<= xb[1] & LogMean > xb[2],LogVar])
	   	      
	   	      return(params_for_min_calc[LogMean <= xb[1] & LogMean > xb[2] & (LogVar <= min_cutoff),.(LogMean,LogVar,SequenceType)])
	    }))
	    x_mins
	}
	
	datapoints_per_step = ifelse(nrow(x_bio) > 500,150:250,floor(nrow(x_bio)/10):floor(nrow(x_bio)/8))
	x_mins = rbind(run_min(datapoints_per_step,x_bio,x_bio),run_min(1:2,x_tech,x_tech_range_all))
	spl_basis = bs(x_mins$LogMean,df=4)
	   res = lm(LogVar ~ bs(LogMean,df=4), data = as.data.table(x_mins))

	   rng = 2^range(x_mins$LogMean)
	   x_within_range = new_x>=rng[1] & new_x <= rng[2]
	   pred_res = rep(NA,length(new_x));
	   pred_res[x_within_range] = predict(res,newdata=data.frame(LogMean = log2(new_x[x_within_range])))
	   overall_max = max(pred_res[x_within_range])
	   overall_max_i = min(which(pred_res[x_within_range] == overall_max)) #smallest mean that reaches the max val
	   force_monotonic = which(x_within_range);
	   force_monotonic = force_monotonic[force_monotonic>=overall_max_i]
	   pred_res[force_monotonic]= overall_max
	   
	   names(pred_res) = paste0("x",as.character(new_x));
	   #browser()
	   pred_res;
	}

#given a set of the local minimum Respdisp, generate a smooth predictor that can be used to estimate the "baseline" LogVar at any point.
generate_predictive_model_from_min_fit = function(fit){
	new_x = c(1,1.5,2,2.5,3,3.5,4,5,6,7,10,25,50,75,100,250,500,750,1000,2500,5000,7500,10000,25000,50000,75000,100000,250000,500000)
	if (length(fit) != length(new_x)) stop("fit must be the running minimum at preset points")
	dat = data.frame(LogMean = log2(new_x),LogVar= fit)
	dat = dat[!is.na(dat$LogVar),] 
	smooth.spline(dat$LogMean,dat$LogVar)
}


#calculate_separate-tech: should we calculate the residual overdispersion in respect to spike-ins unique to each co-variate group?
#Set to F to fit to the combined spike-in set across all groups
estimateResVar <- function(bio_counts, tech_counts, 
				sample_list, 
				groups,
				gene_list=NULL,
				calculate_separate_tech = F,
				trend_fit_type = "running_min")
				{
 
  
  .extractDESeqParams <- function(dds) {
    log2_mu <- drop(coef(dds))
    log2_delta <- log2(dispersions(dds))
    row_data <- rowData(dds)
    data.table(GeneName = rownames(dds), 
               AllZero = row_data$allZero, 
               BaseMean = row_data$baseMean, 
               BaseVar = row_data$baseVar, 
               LogMean = log2_mu, 
               LogDisp = log2_delta)
  }
  if (trend_fit_type == "SSasymp"){
	  .computeResVar <- function(fit, newdata) {
	    log_delta_hat <- as.vector(predict(fit, newdata = newdata))
	    epsilon = 2^(newdata$LogVar - log_delta_hat)-1 #we subtract one so that epsilon == 0 indicates no variation (in addition to spike-ins)
	   epsilon2=ifelse(is.na(epsilon),0,ifelse(epsilon>0,log2(epsilon),log2(.01)))  #minimum floor of stdev = .01
	
	    cbind(newdata, ResVar = epsilon2)
	  } 
  }else if (trend_fit_type == "running_min"){
  	.computeResVar <- function(fit, newdata) {
  	    log_delta_hat <- as.vector(predict(fit, newdata$LogMean)$y)
  	   
  	    epsilon <- 2^(newdata$LogVar - log_delta_hat)-1 #we subtract one so that epsilon == 0 indicates no variation (in addition to spike-ins)
  
	    epsilon2=ifelse(is.na(epsilon),0,ifelse(epsilon>0,log2(epsilon),log2(.01)))  #minimum floor of stdev = .01    
	 
  	    cbind(newdata, ResVar = epsilon2)
	  } 
  }else stop(paste("Unknown trend fitting type specified: ",trend_fit_type))
  
   #note that we need to use the whole transcriptome to generate nf factors
   #we can't use just the subset specified in gene_list
   bio_nf <- normalizationFactor(bio_counts, groups = groups)
   tech_nf <- normalizationFactor(tech_counts, groups = groups)
   
   if (!calculate_separate_tech)
   	consensus_tech_nf <- normalizationFactor(tech_counts,groups=rep(1,ncol(tech_counts)))

   if (!is.null(gene_list)){
   	
   	if (length(gene_list) == 0) {
   		warning("No genes selected in this round");
   		return(NULL)
   	}
   	bio_counts = bio_counts[gene_list,]
   }
cat(paste("Fitting",nrow(bio_counts),"genes\n"))
flush.console()


      #use MLE estimation to fit
       too_few_nonzero = function(x){length(which(x==0))>.90*length(x)}

       null_result = c(0,0,0); 
       names(null_result) = c("size","mu","");

       
       .fit_mu_delta = function(counts_to_use,nf,samples_list,sequence_type){
	       lapply( sample_list, function(samples) {
		   if (length(samples) == 1){
			stop("Encountered a sample set with only one sample.  Perhaps one of the groups only has a single sample?")
		   }
		   counts = normalizeCounts(counts_to_use[,samples], nf[samples])

		  # print("Fitting")
		   bio_list = t(apply(counts,1,function(x){
			   if (too_few_nonzero(x)) return(null_result);
			   tryCatch({
				f = MASS::fitdistr(as.integer(x),"negative binomial", lower=c(.1,10^-10),upper=c(max(x),Inf))
				c(f$estimate,f$loglik)
				},error=function(err){print(err);return(null_result)})
			}))
		   bad_fit = bio_list[,"size"] == 0;
		   LogVar = ifelse(bad_fit,0,
				   log2(1/bio_list[,"size"]*(bio_list[,"mu"]^2)+bio_list[,"mu"]) #from definition in dnbinom used by fitdistr()
			    )
		   data.table(GeneName = rownames(bio_list), 
			 AllZero = apply(counts,1,too_few_nonzero),
			 BaseMean = apply(counts,1,mean),
			 BaseVar = apply(counts,1,var),
			 bad_fit = bad_fit,
			 LogMean = ifelse(bad_fit,0,log2(bio_list[,"mu"])), 
			 LogDisp = ifelse(bad_fit,0,log2(1/bio_list[,"size"])),
			 LogVar = LogVar,
			 LogLik = ifelse(bad_fit,2,bio_list[,3]/log(2)),
			 SequenceType = sequence_type)
	
		})
	}
	bio_mu_delta <- .fit_mu_delta(bio_counts,bio_nf,samples_list,"Biological");
	tech_mu_delta <- .fit_mu_delta(tech_counts,tech_nf,samples_list,"Technical");
	if (!calculate_separate_tech)
	   consensus_tech_mu_delta = .fit_mu_delta(tech_counts,consensus_tech_nf,list(colnames(tech_counts)),"Technical")[[1]];
	
  
 	
	  if (!calculate_separate_tech){
		  print(paste("Fitting",trend_fit_type,"model using consensus spike-ins"))
		  tech_params_to_use_for_fitting = list(consensus_tech_mu_delta)
	  }else{
		tech_params_to_use_for_fitting = tech_mu_delta
		print(paste("Fitting",trend_fit_type,"model using covariate-group-specific spike-ins"))
	  }
	  flush.console()
	  if (trend_fit_type == "SSasymp"){
		  # fit overdispersion - mean trend using SSasymp model
		  fit_list <- lapply(tech_params_to_use_for_fitting, function(params) {
		    good_params = params[!(is.infinite(LogDisp) & LogLik == 0),]
		    iter <- 0
		    max_iter <- 100000
		    fit <- NULL
		    while(is.null(fit) & iter <= max_iter) {
		      iter <- iter + 1
		      init <- list(Asym  = - rlnorm(1, log(9), 0.4),
				   resp0 = rlnorm(1, log(0.6), 0.4),
				   lrc   = - rlnorm(1, log(2.4), 0.4))
		      fit <- tryCatch(nls(LogVar ~ Asym + (resp0-Asym) * exp(-exp(lrc)*LogMean), 
					  data = as.data.table(good_params), 
					  start = init),
				      error = function(e) {
					return(NULL)
				      })
		    }
		    return(list(fit=fit,coef=coef(fit)))
		})
	}else{
	    fit_list <- lapply(1:length(tech_params_to_use_for_fitting), function(i) {
	    	#this approach uses both mRNA and spike-ins for fitting
	    	good_params = rbind(bio_mu_delta[[i]],tech_params_to_use_for_fitting[[i]])
		good_params = good_params[!(is.infinite(LogVar) & LogLik == 0),]
		fit = fit_bio_and_tech_running_minimum(good_params);
		return(list(fit=generate_predictive_model_from_min_fit(fit),coef=fit));
	        })
	}
	

	if (!calculate_separate_tech){
		#we repeat the consensus tech each time.
		fit_list_tmp = fit_list
		fit_list <- lapply(tech_mu_delta,function(params){fit_list[[1]]})
	}
	names(fit_list) = names(sample_list)

	  if (any(sapply(fit_list, is.null))) {
	    return(list(params = NULL, trend = NULL))
	  }
	  
	  fit_params <- t(sapply(fit_list, function(x)x[["coef"]]))
	  fit_params_cnames = colnames(fit_params)
	  fit_params <- cbind(Group = names(fit_list), as.data.table(fit_params))
	  

	  bio_params <- lapply(seq_along(fit_list), function(j) {
	    .computeResVar(fit = fit_list[[j]]$fit, newdata = bio_mu_delta[[j]])
	  })
	  names(bio_params) <- names(fit_list)

	  tech_params <- lapply(seq_along(fit_list), function(j) {
	    .computeResVar(fit = fit_list[[j]]$fit, newdata = tech_mu_delta[[j]])
	  })
	  names(tech_params) <- names(fit_list)
	  if (!calculate_separate_tech){
		  consensus_tech_params <- lapply(seq_along(fit_list), function(j) {
		      .computeResVar(fit = fit_list[[j]]$fit, newdata = consensus_tech_mu_delta)
		    })
		  names(consensus_tech_params) <- names(fit_list)
	  }
	  all_params <- rbindlist(lapply(names(fit_list), function(group) {

	    col_order <- c("Group", "GeneName", "SequenceType", 
			   "BaseMean", "BaseVar", "AllZero",   
			   "LogMean", "LogDisp", "LogVar","LogLik","ResVar")

	    bio_p <- bio_params[[group]]
	    bio_p[, Group := group]
	    setcolorder(bio_p, col_order)

	    tech_p <- tech_params[[group]]
	    tech_p[, Group := group]
	    tech_p[, SequenceType := "Technical"]
	    setcolorder(tech_p, col_order)

	    if (!calculate_separate_tech){
		    consensus_tech_p <- consensus_tech_params[[group]]
		    consensus_tech_p[, Group := group]
		    consensus_tech_p[, SequenceType := "Consensus Technical"]
		    setcolorder(consensus_tech_p, col_order)
		    return(rbind(bio_p, tech_p,consensus_tech_p))
	    }else   return(rbind(bio_p, tech_p))

	  }))
  
  all_params[, Group := factor(Group)]
  all_params[, GeneName := factor(GeneName)]
  all_params[, SequenceType := factor(SequenceType)]
  
  list(params = all_params, trend = fit_params)
}

summarizeTrend <- function(trend,trend_fit_type = "running_min") {
   if (trend_fit_type == "SSasymp"){
	  summary_trend <- trend[, .(
	    Asym   = wtd.mean(Asym, N, na.rm = TRUE),
	    AsymSE = sqrt(wtd.var(Asym, N, na.rm = TRUE)),
	    resp0   = wtd.mean(resp0, N, na.rm = TRUE),
	    resp0SE = sqrt(wtd.var(resp0, N, na.rm = TRUE)),
	    lrc   = wtd.mean(lrc, N, na.rm = TRUE),
	    lrcSE = sqrt(wtd.var(lrc, N, na.rm = TRUE))),
	    by = .(Strain, Food, Day, Temperature,Pooled)]
  }else if (trend_fit_type == "running_min"){
	  summary_trend_mean = as.data.frame(trend[,-'Group'][, lapply(.SD, mean, na.rm = TRUE),by=.(Strain, Food, Day, Temperature,Pooled)])
	  summary_trend_sd = as.data.frame(trend[,-'Group'][, lapply(.SD, function(x)sqrt(var(x, na.rm = TRUE))),by=.(Strain, Food, Day, Temperature,Pooled)])
	  frac_na = trend[,-'Group'][, lapply(.SD, function(x)length(which(is.na(x)))/length(x)),by=.(Strain, Food, Day, Temperature,Pooled)]
	  too_high_na_rate = which(frac_na>=.5,arr.ind=T)
	  #don't remove metadata rows
	  too_high_na_rate = too_high_na_rate[too_high_na_rate[,2] >= 6 & too_high_na_rate[,2] <= 34,]
	  summary_trend_mean[too_high_na_rate] = NA;
	  summary_trend_sd[too_high_na_rate] = NA;
	  summary_trend = merge(summary_trend_mean,summary_trend_sd,
				by=c("Strain","Food","Day","Temperature","Pooled"),
				suffixes=c("","SE"))  
  }
  summary_trend
}

summarizeResVar <- function(params) {
 # params[,bad_fit:= (is.infinite(LogDisp) & LogLik == 0)]#All zeros indicate a model-fitting error; separate
  
  summary_params <- params[, .(
    BaseMean          = mean(BaseMean, na.rm = TRUE),
    BaseVar           = mean(BaseVar, na.rm = TRUE),
    AllZeroProp       = mean(AllZero, na.rm = TRUE),
    MissingValuesProp = mean(is.na(ResVar)),
    LogMean   = wtd.mean(LogMean, N, na.rm = TRUE),
    LogMeanSE = sqrt(wtd.var(LogMean, N, na.rm = TRUE)),
    LogDisp   = wtd.mean(LogDisp, N, na.rm = TRUE),
    LogDispSE = sqrt(wtd.var(LogDisp, N, na.rm = TRUE)),
    LogVar   = wtd.mean(LogVar, N, na.rm = TRUE),
    LogVarSE = sqrt(wtd.var(LogVar, N, na.rm = TRUE)),
    ResVar  = wtd.mean(ResVar, N, na.rm = TRUE),
    ResVarSE = sqrt(wtd.var(ResVar, N, na.rm = TRUE)),
    LogLik   = wtd.mean(LogLik, N, na.rm = TRUE),
    LogLikSE = sqrt(wtd.var(LogLik, N, na.rm = TRUE))),
    by = .(GeneName, SequenceType, Strain, Day, Food, Temperature,Pooled,bad_fit)]
  summary_params[, ResVarPValue := 2 * pnorm(q = abs(ResVar), mean = 0, sd = ResVarSE, lower.tail = FALSE)]
  summary_params[, ResVarAdjPValue := p.adjust(ResVarPValue, "fdr")]
  
  summary_params
}


# run regression to perform differential variability analysis
runDifferentialVariabilityAnalysis <- function(all_params, formula) {
  
  setkey(all_params, GeneName, BootIndex, Group)
  
  all_genenames <- as.character(sort(unique(all_params$GeneName)))
  gene_N = length(all_genenames)
  split_groups = factor(floor(1:gene_N/1000))
  
  gene_groups = split(1:gene_N, split_groups)
  print(paste0("Differential analysis for ",cur_group, " with ", max(all_params$BootIndex), " bootstrap replicates; ",ifelse(calculate_separate_tech,"covariate-specific spike-ins","common spike-ins")))
  # iterate through gene sets 1000 genes at a time, to keep memory usage down
  return(rbindlist(pblapply(1:length(gene_groups),function(split_group){

	  genenames <- all_genenames[gene_groups[[split_group]]]
	  params = all_params[GeneName %in% genenames,]
	  nboots <- max(params$BootIndex)

	  formula <- stats::update(formula, ResVar ~ .)
	 
	  err_debug = 0;
	  dva_df <- rbindlist(lapply(genenames, function(genename) {

	    results <- rbindlist(lapply(seq_len(nboots), function(boot_index) {
	    
	      dat <- params[.(genename, boot_index),]
	      X <- model.matrix(formula, dat)
	      colnames(X) <- gsub("\\(|\\)", "", colnames(X))
	      #We can also check dat$AllZero
	      #as all the NAs are generated by genes for which not enough counts are available.
	      not_na <- ! is.na(dat$ResVar)
	  
	      fit = NULL;
	      if (length(which(not_na))>1){
	      	      #note that model_matrix automagically removes rows that are NA, so we don't
	      	      #need to subset it.
		      fit <- tryCatch(lm.wfit(X, dat$ResVar[not_na], w = dat$N[not_na]),
				      error = function(e) NULL)
	      }
	      if (! is.null(fit)) {
		coeff <- coef(fit)
		err_debug=0
	      } else {
	      	if (length(which(not_na))>1){
	      		print("unexpected data fitting err");
	      	}
	      	err_debug=err_debug+1;
	      	if(err_debug == 250){
	      		print("Many data fitting errors in a row!");
	      	}
	   
		coeff <- setNames(rep_len(NA_real_, ncol(X)), colnames(X))
	      }

	      coeff_df <- as.data.table(t(coeff))
	      coeff_df[, BootIndex := boot_index]
	      coeff_df[, MissingValueProp := 1 - sum(dat$N[not_na]) / sum(dat$N)]

	      coeff_df
	    }))
	    results[, GeneName := genename]
	    results
	  }))
	 
	  X <- model.matrix(formula, params[1,])
	  coeffnames = gsub("\\(|\\)", "", colnames(X))
	  summary_dva_df <- lapply(coeffnames, function(coeffname) {
	    out <- dva_df[, .(logFoldChange = mean(get(coeffname), na.rm = TRUE),
			      lfcSE = sd(get(coeffname), na.rm = TRUE)), 
			  by = .(GeneName)]
	    out[, pvalue := 2 * pnorm(q = abs(logFoldChange), mean = 0, sd = lfcSE, lower.tail = FALSE)]
	    out[, padj := p.adjust(pvalue, "fdr")]
	    colnames(out)[-1] <- paste0(coeffname, "_", colnames(out)[-1])
	    out
	  })
	  summary_dva_df <- Reduce(function(x, y) merge(x, y, by = "GeneName"),
				   summary_dva_df)

	  summary_dva_df <- merge(dva_df[, .(MeanMissingValueProp = mean(MissingValueProp)), by = .(GeneName)], 
				  summary_dva_df, 
				  by = "GeneName", 
				  all = TRUE)
	  summary_dva_df
	})))
	
}


ns_fit_set_da_linear_model = function (counts,differential_analysis,differential_analysis_column_name = "log2FoldChange",RNAi_GeneName_to_exclude, keep_x=T,
                                   min_mean_count_for_fitting = 25,allow_recursion_for_outliers=0,keep_residuals=T,reference_group=NA){
  if (class(counts) != "list")
    stop("counts must be a list");
  
  offset = -1
  
  if (is.na(reference_group)){
    sample_sizes = unlist(lapply(counts,ncol))
    sums = do.call("cbind",lapply(counts,function(x)apply(x,1,function(x)sum(log2(x-offset)))))
    grand_mean = apply(sums,1,function(x)sum(x))/sum(sample_sizes)
  }else{
    grand_mean = rowMeans(log2(counts[[reference_group]]-offset))
  }
  res = lapply(counts,function(x){
    ns_fit_da_linear_model(x,
                            differential_analysis=differential_analysis,
                            differential_analysis_column_name=differential_analysis_column_name,
                            RNAi_GeneName_to_exclude=RNAi_GeneName_to_exclude,
                            keep_x=keep_x,
                            min_mean_count_for_fitting=min_mean_count_for_fitting,
                            allow_recursion_for_outliers=allow_recursion_for_outliers,
                            keep_residuals=keep_residuals,
                            preset_log_means=grand_mean)
    })
  names(res)=names(counts);
  res;
}

#take a counts matrix with i samples, and a differential analysis, V, which is a vector of log2FoldChanges, and
#fit the model log(y_i) = beta_i*V + e_i with y_i, V, and e_i being vectors of length G where g is the number of genes measured
ns_fit_da_linear_model = function (counts,
differential_analysis,
differential_analysis_column_name = "log2FoldChange",
RNAi_GeneName_to_exclude, 
genes_to_fit = c(),keep_x=T,
	min_mean_count_for_fitting = 25,
	min_min_count_for_fitting = 0,
	allow_recursion_for_outliers=0,
	keep_residuals=T,
	preset_log_means=NA,
	only_return_stats = F,
	log_transform=T,
	cor_type = "s",
	debug=c()){

  if (log_transform){
    offset = -1
    log_counts = log2(counts-offset)
  }else{
    offset = 0
    log_counts = counts
  }
  if (length(preset_log_means)==1 && is.na(preset_log_means)){
    log2_gene_means = apply(log_counts,1,mean)
  }else{
    log2_gene_means = preset_log_means
  }
  log_scaled = log_counts-log2_gene_means
  
  common_genes = intersect(rownames(log_scaled),differential_analysis$GeneName);
  if (!any(is.na(RNAi_GeneName_to_exclude))){
    common_genes = common_genes[! common_genes %in% RNAi_GeneName_to_exclude]
  }else warning("RNAi_GeneName_to_exclude is NA")
  
  
  setkey(differential_analysis,GeneName)
  cur_scaled = log_scaled[common_genes,]
  cur_da = differential_analysis[common_genes,]
  if (any(table(cur_da$GeneName)>1)){
    stop("Some genes appear more than once in the differential analysis data") 
  }
  if ("data.table" %in% class(cur_da)){
    cur_da_c = cur_da[,get(differential_analysis_column_name)]
  }else 
    cur_da_c = cur_da[,differential_analysis_column_name]
   
   if (cor_type == "p"){
   	ns_min_corr = function(x,dat,vec,ordered_vec)abs(cor(dat-x*vec,vec,method="p"))
   }else if (cor_type == "s"){
   	ns_min_corr = function(x,dat,vec,rank_vec)abs(cor(frank(dat-x*vec),rank_vec,method="p"))
   	#ns_min_corr = function(x,dat,vec,rank_vec)abs(cor(dat-x*vec,vec,method="s"))
   }
  
  #don't fit the model on low-expressed genes
  if (min_mean_count_for_fitting>0){
    abundant_genes = names(which(rowMeans(counts[common_genes,])>=min_mean_count_for_fitting))
   }else abundant_genes = common_genes
   if (min_min_count_for_fitting>0){
    abundant_genes <- intersect(abundant_genes,rownames(counts)[apply(counts,1,function(cnts2)!any(cnts2 < min_min_count_for_fitting))])
   }
   if (length(genes_to_fit)>0){
   	in_common = intersect(genes_to_fit,abundant_genes);
   	if (length(in_common) < 500){
   		stop(paste0("genes_to_fit and the supplied counts and DA data have few genes in common: ",length(in_common)))
   	}
   	genes_to_fit = 	in_common
  }else genes_to_fit = abundant_genes
   
   genes_to_fit =  common_genes %in% genes_to_fit 
   
  qt = quantile(abs(cur_da_c),.98)
  strong_da = abs(cur_da_c) >= qt
  if (length(which(qt&strong_da))/length(which(strong_da)) < .9)
  	stop("Less than 90% of the most differentially expressed transcipts have counts greater than the specified min_min_count_for_fitting and min_mean_count_for_fitting thresholds")
  if (sd(cur_da_c[genes_to_fit]) == 0){warning("Standard deviation of DA is zero"); 
    ft = rep(0,ncol(cur_scaled))
  }else{
      #print(paste("ignoring",length(which((!genes_to_fit))),"genes"))
        #if we are using spearman correlation, we can precomupte the rank variables of the DA vector 
      	  #to ~halve computational time
      	subsetting_genes = any(!genes_to_fit)
      	if(!subsetting_genes){
		vec_to_use = cur_da_c
	}else{
		vec_to_use = cur_da_c[genes_to_fit]
      	}
      	if (cor_type == "s"){
      		rank_vec = frank(vec_to_use)
      	}else{
      		rank_vec = c()
      	}
       	
	to_fit = colnames(cur_scaled)
	if (length(debug)>0) to_fit = debug
      ft = unlist(pbsapply(to_fit,function(cn){
      	  x = cur_scaled[,cn]
      	  if(!subsetting_genes){
      	  	x_to_use = x
      	  }else{
      	  	x_to_use= x[genes_to_fit]
      	  }
          if (sd(x_to_use) == 0) warning("Standard deviation of a sample is zero");
         
	start_locations = seq(-1,1,by=.5)
	p = sapply(start_locations,function(init){
		unlist(optim(par=init,fn=ns_min_corr,dat=x_to_use,vec=vec_to_use,rank_vec=rank_vec,lower=-5,upper=5,method="Brent"))
	})
	best_fits = which(p["value",] == min(p["value",]))
	r = min(p["par",best_fits])
	
	return(min(p["par",best_fits]))
	
          if(length(debug)>0)
          	browser()
          p$par
      }))

  }
  resid = matrix(NA,nrow=nrow(cur_scaled),ncol=ncol(cur_scaled))
  rownames(resid) = common_genes;
  colnames(resid) = colnames(log_counts);
  effect_v = matrix(NA,nrow=nrow(cur_scaled),ncol=ncol(cur_scaled))
  rownames(effect_v) = common_genes;
  colnames(effect_v) = colnames(log_counts);

  for (i in 1:ncol(cur_scaled)){
    resid[,i] = cur_scaled[,i]-ft[i]*cur_da_c
  }  
  for (i in 1:ncol(cur_scaled)){
    effect_v[,i] = ft[i]*cur_da_c
  }
  interquant_50 = function(x){
    unname(diff(quantile(x,c(.25,.75))) / 1.349)
  }
  interquant_90 = function(x){
    unname(diff(quantile(x,c(.05,.95))))
  }
 # resid_abs = (resid+log2_gene_means)
  if (log_transform){
    rfun = function(x)2^x;
  }else rfun = function(x)x
  
  #calculate stats for all genes
  orig_sd = apply(rfun(cur_scaled+log2_gene_means[common_genes]),1,var)
  resid_sd = apply(rfun(resid+log2_gene_means[common_genes]),1,var)
  effect_sd = apply(rfun(effect_v+log2_gene_means[common_genes]),1,var)
  
  orig_IQD50 = apply(rfun(cur_scaled+log2_gene_means[common_genes]),1,interquant_50)
  resid_IQD50 = apply(rfun(resid+log2_gene_means[common_genes]),1,interquant_50)
  effect_IQD50 = apply(rfun(effect_v+log2_gene_means[common_genes]),1,interquant_50)
  
  orig_IQD90 = apply(rfun(cur_scaled+log2_gene_means[common_genes]),1,interquant_90)
  resid_IQD90 = apply(rfun(resid+log2_gene_means[common_genes]),1,interquant_90)
  effect_IQD90 = apply(rfun(effect_v+log2_gene_means[common_genes]),1,interquant_90)
   
  
  if (allow_recursion_for_outliers > 0){
    ord = resid_IQD50[order(resid_IQD50,decreasing=T)]
    resid_to_remove = ord[1]
    if (resid_to_remove > 3*ord[allow_recursion_for_outliers+1]){
    #allow_recursion_for_outliers = 0
      warning("Found a large outlier.  Running again dropping that gene.")
      #browser()
      res = ns_fit_da_linear_model(counts=counts,
                                   differential_analysis=differential_analysis,
                                   differential_analysis_column_name=differential_analysis_column_name,
                                   RNAi_GeneName_to_exclude = c(names(resid_to_remove),RNAi_GeneName_to_exclude),
                                   min_mean_count_for_fitting=min_mean_count_for_fitting,
                                   min_min_count_for_fitting=min_min_count_for_fitting,
                                   allow_recursion_for_outliers=allow_recursion_for_outliers-1,
                                   keep_x=keep_x,
                                   keep_residuals=keep_residuals,
				   only_return_stats=only_return_stats);
      if (!only_return_stats)
        res$outliers_removed = c(names(resid_to_remove),res$outliers_removed)
      return(res);
    }
  } 
  good_IQD_ests = orig_IQD50 != 0
  good_sd = orig_sd != 0
  abundant_gene_lookup = common_genes %in% abundant_genes ;
  stats = c(all_genes_relative_total_variance = sum(resid_sd[good_sd])/sum(orig_sd[good_sd]),
            all_genes_relative_per_gene_variance = mean(resid_sd[good_sd]/orig_sd[good_sd]),
            all_genes_relative_total_IQD50 = sum(resid_IQD50[good_IQD_ests])/sum(orig_IQD50[good_IQD_ests]),
            all_genes_relative_per_gene_IQD50 = mean(resid_IQD50[good_IQD_ests]/orig_IQD50[good_IQD_ests]),
            all_genes_relative_total_IQD90 = sum(resid_IQD90[good_IQD_ests])/sum(orig_IQD90[good_IQD_ests]),
            all_genes_relative_per_gene_IQD90 = mean(resid_IQD90[good_IQD_ests]/orig_IQD90[good_IQD_ests]),
            all_genes_effect_relative_total_variance = sum(effect_sd[good_sd])/sum(orig_sd[good_sd]),
            all_genes_effect_relative_per_gene_variance = mean(effect_sd[good_sd]/orig_sd[good_sd]),
            all_genes_effect_relative_total_IQD50 = sum(effect_IQD50[good_IQD_ests])/sum(orig_IQD50[good_IQD_ests]),
            all_genes_effect_relative_per_gene_IQD50 = mean(effect_IQD50[good_IQD_ests]/orig_IQD50[good_IQD_ests]),
            all_genes_effect_relative_total_IQD90 = sum(effect_IQD90[good_IQD_ests])/sum(orig_IQD90[good_IQD_ests]),
            all_genes_effect_relative_per_gene_IQD90 = mean(effect_IQD90[good_IQD_ests]/orig_IQD90[good_IQD_ests]),
            all_genes_total_resid_variance=sum(resid_sd[good_sd]),
            all_genes_total_observation_variance=sum(orig_sd[good_sd]),
            all_genes_total_effect_variance=sum(effect_sd[good_sd]),
            
            relative_total_variance = sum(resid_sd[abundant_gene_lookup & good_sd])/sum(orig_sd[abundant_gene_lookup & good_sd]),
	    relative_per_gene_variance = mean(resid_sd[abundant_gene_lookup & good_sd]/orig_sd[abundant_gene_lookup & good_sd]),
		relative_total_IQD50 = sum(resid_IQD50[abundant_gene_lookup & good_IQD_ests])/sum(orig_IQD50[abundant_gene_lookup & good_IQD_ests]),
		relative_per_gene_IQD50 = mean(resid_IQD50[abundant_gene_lookup & good_IQD_ests]/orig_IQD50[abundant_gene_lookup & good_IQD_ests]),
		relative_total_IQD90 = sum(resid_IQD90[abundant_gene_lookup & good_IQD_ests])/sum(orig_IQD90[abundant_gene_lookup & good_IQD_ests]),
		relative_per_gene_IQD90 = mean(resid_IQD90[abundant_gene_lookup & good_IQD_ests]/orig_IQD90[abundant_gene_lookup & good_IQD_ests]),
		effect_relative_total_variance = sum(effect_sd[abundant_gene_lookup & good_sd])/sum(orig_sd[abundant_gene_lookup & good_sd]),
		effect_relative_per_gene_variance = mean(effect_sd[abundant_gene_lookup & good_sd]/orig_sd[abundant_gene_lookup & good_sd]),
		effect_relative_total_IQD50 = sum(effect_IQD50[abundant_gene_lookup & good_IQD_ests])/sum(orig_IQD50[abundant_gene_lookup & good_IQD_ests]),
		effect_relative_per_gene_IQD50 = mean(effect_IQD50[abundant_gene_lookup & good_IQD_ests]/orig_IQD50[abundant_gene_lookup & good_IQD_ests]),
		effect_relative_total_IQD90 = sum(effect_IQD90[abundant_gene_lookup & good_IQD_ests])/sum(orig_IQD90[abundant_gene_lookup & good_IQD_ests]),
		effect_relative_per_gene_IQD90 = mean(effect_IQD90[abundant_gene_lookup & good_IQD_ests]/orig_IQD90[abundant_gene_lookup & good_IQD_ests]),
		total_resid_variance=sum(resid_sd[abundant_gene_lookup & good_sd]),
		total_observation_variance=sum(orig_sd[abundant_gene_lookup & good_sd]),
            all_genes_total_effect_variance=sum(effect_sd[abundant_gene_lookup & good_sd])
    )
 # browser()
  if (any(is.infinite(stats))){
    cat("Found Inifinte!");
  }
  if (only_return_stats){
    return(stats)
  }
  
  x = NA;
  residuals = NA;
  log2_intercepts = NA
  if (keep_x){
    x = rfun(cur_scaled+log2_gene_means[common_genes]);
  }
  if (keep_residuals){
    residuals = resid;
    log2_intercepts = log2_gene_means[common_genes]
  }else{
    residuals = data.frame(residuals_var=resid_sd,
                           orig_var=orig_sd,
                           effect_var = effect_sd,
                           residuals_IQD50 =resid_IQD50,
                           orig_IQD50=orig_IQD50,
                           effect_IQD50=effect_IQD50,
                           residuals_IQD90 =resid_IQD90,
                           orig_IQD90=orig_IQD90 ,
                           effect_IQD90=effect_IQD90);
    rownames(residuals) = names(resid_sd);
  }
 
  
    return(list(coef=ft,
           log2_intercepts = log2_intercepts,
           residuals = residuals, 
           genes_used_to_fit = abundant_genes,
           all_genes_in_common = common_genes,
           x = x,
           stats = stats,
           offset=offset,
           min_mean_count_for_fitting = min_mean_count_for_fitting,
           min_min_count_for_fitting = min_min_count_for_fitting,
           outliers_removed = c()
    ))
  
}


ns_fit_da_linear_model_bootstrapped = function (bootstrap_N,counts,differential_analysis,differential_analysis_column_name = "log2FoldChange",
						RNAi_GeneName_to_exclude, keep_x=F,
						min_mean_count_for_fitting = 25,
						min_min_count_for_fitting = 0,
						allow_recursion_for_outliers=0,
						preset_log_means=NA,keep_residuals=F){
  print(Sys.time())
  cat("Fitting all observations...\n")
  cat(paste("Using mean gene expression threshold",min_mean_count_for_fitting,"\n"))
  cat(paste("Using min gene expression threshold",min_min_count_for_fitting,"\n"))
  mean_est = ns_fit_da_linear_model(counts,
                                 differential_analysis=differential_analysis,
                                 differential_analysis_column_name = differential_analysis_column_name,
                                 RNAi_GeneName_to_exclude = RNAi_GeneName_to_exclude, 
                                 keep_x=F,
                                 min_mean_count_for_fitting = min_mean_count_for_fitting,
                                 min_min_count_for_fitting=min_min_count_for_fitting,
                                 allow_recursion_for_outliers=allow_recursion_for_outliers,
                                 keep_residuals = keep_residuals,
                                 only_return_stats=F)
   #we only allow recursion and outlier removal in the first pass, so all bootstraps run on the same genes.          
   if (is.null(mean_est$outliers_removed)){
   	outliers_to_remove = c();
   }else{
  	outliers_to_remove = mean_est$outliers_removed        
  }

  if (is.na(RNAi_GeneName_to_exclude)){
  	  genes_to_exclude_during_bootstrapping = outliers_to_remove;
  }else genes_to_exclude_during_bootstrapping = c(outliers_to_remove,RNAi_GeneName_to_exclude)

  print(Sys.time())
  cat("\n")
  boot_data = sapply(1:bootstrap_N,function(x){
  	print(Sys.time())
    cat(paste("Boot",x,"..."))
    dat = counts[,sample.int(ncol(counts),replace=T)]
  
    ns_fit_da_linear_model(dat,
                            differential_analysis=differential_analysis,
                            differential_analysis_column_name = differential_analysis_column_name,
                            RNAi_GeneName_to_exclude = genes_to_exclude_during_bootstrapping, 
                            keep_x=F,
                            min_mean_count_for_fitting = min_mean_count_for_fitting,
                            min_min_count_for_fitting=min_min_count_for_fitting,
                            allow_recursion_for_outliers=0,
                            only_return_stats=T
    )
  })
  cat("\nDone\n")
  print(Sys.time())
  qt = apply(boot_data,1,FUN=quantile,probs=c(.025,.05,.95,.975))
  ci = 2*mean_est$stats-t(qt)
  colnames(ci) = rev(c("2.5%","5%","95%","97.5%"))
  pv = rowSums(boot_data>=1)/ncol(boot_data)
  r = t(boot_data);
  rownames(r) = NULL;
  list(dat = mean_est,
  	outliers_removed = outliers_to_remove,
       stats = data.frame(mean_est=mean_est$stats,ci=ci,pv=pv),
       boot_data = boot_data
  )
}