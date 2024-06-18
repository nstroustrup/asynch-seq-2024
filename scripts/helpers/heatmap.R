plotSplitCorrelationHeatmap <- function(...) {
   return(plotSplitCorrelationHeatmap_internal(...)[["heatmap"]]);                                
}

plotSplitCorrelationHeatmap_internal <- function(rho,
                                        split_list=NULL, 
                                        split_colors = NULL,
                                        beta = 2,name=NULL,
                                        difference=F,
                                        dend_spec=NULL,
                                        gene_labels=NULL,
                                        show_dendrogram=T,
                                        show_triangle="both",
                                        border=F) {
  # define distance function
  if (!difference){
    distfun <- function(x) {
      as.dist(1-abs(x[rownames(x), rownames(x)]))
    }
  }else{
    distfun <- function(x) {
      as.dist(x[rownames(x), rownames(x)])
    }
  }
  
  # make sure that genes are in both correlation matrix and split_list
  if (!is.null(split_list)){
    intersect_genes <- intersect(rownames(rho), unlist(split_list))
    rho <- rho[intersect_genes, intersect_genes]
	  split_list <- lapply(split_list, function(x) x[x %in% intersect_genes])

	  # convert split list to vector
	  if (! is.null(names(split_list))) {
	    iter <- names(split_list)
	  } else {
	    iter <- seq_along(split_list)
	  }

	  split <- unlist(lapply(iter, function(i) {
	    out <- rep_len(i, length(split_list[[i]]))
	    names(out) <- split_list[[i]]
	    out
	  }))
  
  	split <- split[rownames(rho)]
  }
  
  # raise correlation to some power
  rho <- sign(rho) * abs(rho)^beta
 # browser()
  # define colors
  if (is.null(split_colors)){
    if (is.null(split_list)){
      split_colors = "black";
    }else{ 
      split_colors <- colorspace::rainbow_hcl(length(split_list)) 
      names(split_colors) <- names(split_list)
    }
  }else  names(split_colors) <- names(split_list);
    
  if (!difference){
    col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("royalblue4", "white", "firebrick"))
  } else col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("royalblue4", "white", "firebrick"))
  
  # heatmap
 if (is.null(dend_spec)){
	  clustered_data =hclust(as.dist(distfun(rho)), method = "single")
	  gene_dendogram = as.dendrogram(clustered_data)
 }else{
	gene_dendogram = dend_spec
 }
 
 if (!is.null(gene_labels)){
 	label_groups = unique(gene_labels)
 	gene_label_color = rainbow(length(label_groups))
 	if(any(label_groups=="")){
 	  gene_label_color[which(label_groups=="")] = "black"
 	}else  if(any(label_groups=="0")){
 	  gene_label_color[which(label_groups==0)] = "black"
 	}else  if(any(label_groups==0))
 	  gene_label_color[which(label_groups==0)] = "black"
 	
 	column_annot = columnAnnotation(foo = anno_text(gene_labels, 
 	                                                location = 0.5, just = "center",
 	                                                gp = gpar(fontsize = 6, fill = gene_label_color[match(gene_labels,label_groups)],
 	                                                          col = "white"),height = max_text_height("1")*2,rot = 90))
 	
 	row_annot = rowAnnotation(foo = anno_text(gene_labels, 
 	                                          location = 0.5, just = "center",
 	                                          gp = gpar(fontsize = 6, fill = gene_label_color[match(gene_labels,label_groups)],
 	                                                    col = "white"),width = max_text_width("1")*2,rot = 90))
 }else{ 
   column_annot = row_annot = NULL
 }

if (show_triangle=="both"){
  cell_fun=NULL
  layer_fun=NULL
  rect_gp = gpar(col = NA)
}else if (show_triangle=="upper"){
  layer_fun = function(j, i, x, y, w, h, fill) {
    var1 = as.numeric(x) >= 1- as.numeric(y)
    grid.rect(x[var1],y[var1],w[var1],h[var1],gp = gpar(fill = fill[var1],col=fill[var1]))
    grid.rect(x[!var1],y[!var1],w[!var1],h[!var1],gp = gpar(fill = "white",col="white"))
  }
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(as.numeric(x) >= 1 - as.numeric(y)) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
    }
  } 
  rect_gp =  gpar(type = "none");
}else if (show_triangle=="lower"){
  layer_fun = function(j, i, x, y, w, h, fill) {
    var1 = as.numeric(x) < 1- as.numeric(y)
    grid.rect(x[var1],y[var1],w[var1],h[var1],gp = gpar(fill = fill[var1],col=fill[var1]))  
    grid.rect(x[!var1],y[!var1],w[!var1],h[!var1],gp = gpar(fill = "white",col="white"))
  }
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(as.numeric(x) <= 1 - as.numeric(y)) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
    }
  }
  rect_gp =  gpar(type = "none");
}else stop("Unknown show triangle option")
  #browser()
  hmp = ComplexHeatmap::Heatmap(rho, 
                          use_raster = TRUE,
                          col = col_fun,
                          layer_fun=layer_fun,
                          rect_gp=rect_gp,
                          row_title_gp = grid::gpar(col = split_colors), 
                          column_title_gp = grid::gpar(col = split_colors), 
                          
                          cluster_rows=gene_dendogram,
                          cluster_columns=gene_dendogram,
                          
                         # row_split = split, 
                         # column_split = split, 
                          top_annotation = column_annot,
                          left_annotation = row_annot,
                          cluster_row_slices = FALSE, 
                          cluster_column_slices = FALSE,
                          show_row_names = FALSE, 
                          show_column_names = FALSE,
                          show_column_dend = show_dendrogram,
                          show_row_dend = show_dendrogram,
                          name=name,row_title=name,column_title=name,
                         border=border)
   return(list(heatmap = hmp, dend = gene_dendogram))
}

#cluster both groups using rho_bl
plotSplitCorrelationHeatmap_dual <- function(rho_bl,
                                             rho_tr,
                                                 cluster_on = "bl",
                                                 split_list=NULL, 
                                                 split_colors = NULL,
                                                 beta = 2,
                                                 difference=F,
                                                 name_bl=NULL,
                                                 name_tr=NULL,
                                                 title=NULL,
                                                 dend_spec=NULL,
                                                 gene_labels=NULL,
                                                 show_dendrogram=T,
                                                 show_triangle="both",
                                                 border = F) {
  # define distance function
  if (!difference){
    distfun <- function(x) {
      as.dist(1-abs(x[rownames(x), rownames(x)]))
    }
  }else{
    distfun <- function(x) {
      as.dist(x[rownames(x), rownames(x)])
    }
  }
  if (any(rownames(rho_bl) != rownames(rho_tr))) stop("bottom and top matrixes must correspond to the same data!")
  
  # make sure that genes are in both correlation matrix and split_list
  if (!is.null(split_list)){
    
    
    if (! is.null(names(split_list))) {
      iter <- names(split_list)
    } else {
      iter <- seq_along(split_list)
    }
    #only plot genes in both split list and the rho vectors
    all_genes_in_split_list = do.call("c",lapply(iter,function(i)split_list[[i]]))
    
    if(cluster_on=="bl"){
      genes_to_plot <- intersect(rownames(rho_bl),all_genes_in_split_list)
    }else 
      genes_to_plot <- intersect(rownames(rho_bl),all_genes_in_split_list)
    
    rho_bl <- rho_bl[genes_to_plot, genes_to_plot]
    rho_tr <- rho_tr[genes_to_plot, genes_to_plot]
    
    split_list <- lapply(split_list, function(x) x[x %in% genes_to_plot])
    for (i in 1:length(split_list))
      for (j in i:length(split_list)){
        if (i==j) next;
        if (length(intersect(split_list[[i]],split_list[[j]]))>0)
          stop("Values are duplicated in the split list.")
      }
    
    # convert split list to vector
    split <- factor(unlist(lapply(iter, function(i) {
      if (length(split_list[[i]]) == 0) return(c());
      out <- rep_len(i, length(split_list[[i]]))
      names(out) <- split_list[[i]]
      out
    })))
    split = split[genes_to_plot]
    split = droplevels(split);
   
   
    
    fit_within_split_groups = T
    if (fit_within_split_groups){
      #get gene order such that members of each group will be nicely clustered
     # browser()
      distfun2 <- function(x) { as.dist(1-x[rownames(x), rownames(x)])}
      distfun3 <- function(x) { as.dist(1-cor(x[rownames(x),]))}
      
      #order the genes according to their final position in the split table--
      #this is important for the internal plotting of complex heatmap
      #(we order by the factor for each gene's entry in the split list)
     
      genes_reordered_by_clustering <- do.call("c",lapply(levels(split), function(cur_split_group) {
              x = split_list[[cur_split_group]] #genes in the current group
              if (length(x) < 2)
                return(x)
              if(cluster_on=="bl")
                group_dat = rho_bl[x,]
              else group_dat = rho_tr[x,]
              
              
              NA_rows = apply(group_dat,2,function(x)(any(!is.na(x))))
              if (any(!NA_rows))
                stop(paste0(length(which(NA_rows==F))," genes provided are NA across all of ",cur_split_group, " in the correlation matrix"))
            
              
              clusters =hclust(as.dist(distfun3(t(group_dat))), method = "complete")
              gene_dendogram = as.dendrogram(clusters)
              res = rownames(group_dat)[order.dendrogram(gene_dendogram)]
              res = res[res%in% x]
      }))
      ord = match(genes_reordered_by_clustering,rownames(rho_bl))
    }else ord = 1:nrow(rho_bl)
    
    
    #browser()
    
    split <- split[rownames(rho_bl)[ord]]
  }
  
 
  #browser()
  # define colors
  if (is.null(split_colors)){
    if (is.null(split_list)){
      split_colors = "black";
    }else{
      split_colors <- colorspace::rainbow_hcl(length(split_list)) 
      names(split_colors) <- names(split_list)
    }
  }else  names(split_colors) <- names(split_list);
 
  if (!difference){
    col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("royalblue4", "white", "firebrick"))
  } else col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("royalblue4", "white", "firebrick"))
  
  if (is.null(split_list)){#if we do not split we cluster across all genes

    if (is.null(dend_spec)){
      if (cluster_on == "bl"){
        clusters =hclust(as.dist(distfun(rho_bl)), method = "complete")
        gene_dendogram = as.dendrogram(clusters)

      }else{
        clusters =hclust(as.dist(distfun(rho_tr)), method = "complete")
        gene_dendogram = as.dendrogram(clusters)
      }
    }else{
      clusters = dend_spec
      gene_dendogram = as.dendogram(dend_spec);
    }
      ord = order.dendrogram(gene_dendogram)
  }
  else{
    gene_dendogram = F
   # ord = 1:nrow(rho_bl)
  }
  
  
  rho_ordered = matrix(NA,nrow=nrow(rho_bl),ncol=ncol(rho_bl))
  diag(rho_ordered) = 1;
  #browser()
  rho_bl_ordered = rho_bl[ord,ord]
  rho_tr_ordered = rho_tr[ord,ord]
                         
  rho_ordered[lower.tri(rho_ordered)] = rho_bl_ordered[lower.tri(rho_bl_ordered)]
  rho_ordered[upper.tri(rho_ordered)] = rho_tr_ordered[upper.tri(rho_tr_ordered)]
  # raise correlation to some power
  rho_ordered <- sign(rho_ordered) * abs(rho_ordered)^beta
  
  row_title = name_bl
  if (!is.null(split_list)){#if we split, we don't list names of individual columns
    row_title = "%s"
  }
  
      if (!is.null(gene_labels)){
        
       # browser()
        gene_labels = gene_labels[ord];
        
        gene_labels = as.character(gene_labels);
        gene_labels[gene_labels=="0"] = "";
        gene_labels[gene_labels=="NG"] = "";
        
        label_groups = unique(gene_labels)
        color = rainbow(length(label_groups))
        
        if(any(label_groups==""))
          color[which(label_groups=="")] = "black"
        
        fill_color = color[match(gene_labels,label_groups)]
        row_title = name_bl
        if (!is.null(split_list)){#if we split, we don't list names of individual columns
         
          gene_labels = rep("",length(gene_labels))
          row_title = "%s"
        }
        
        column_annot = columnAnnotation(foo = anno_text(gene_labels, 
                                              location = 0.5, just = "center",
                                              gp = gpar(fontsize = 4, fill = fill_color,
                                                        col = "white"),height = max_text_height("1")*2,rot = 0))
        row_annot = rowAnnotation(foo = anno_text(gene_labels, 
                                                  location = 0.5, just = "center",
                                                  gp = gpar(fontsize = 4, fill = fill_color,
                                                            col = "white"),width = max_text_width("1")*2,rot = 0))
      }else{
      column_annot = row_annot = NULL;
    }
  if (!is.null(split_list)){
  hmp = ComplexHeatmap::Heatmap(rho_ordered, 
                                use_raster = TRUE,
                                border = border,
                                col = col_fun,
                                #row_title_gp = grid::gpar(col = split_colors), 
                                #column_title_gp = grid::gpar(col = split_colors), 
                                
                                cluster_rows=F,
                                cluster_columns=F,
                                
                                row_split = split, 
                                column_split = split, 
                                top_annotation = column_annot,
                                left_annotation = row_annot,
                                cluster_row_slices = F, 
                                cluster_column_slices = F,
                                show_row_names = FALSE, 
                                show_column_names = FALSE,
                                show_column_dend = show_dendrogram,
                                show_row_dend = show_dendrogram,
                                name=title,
                                row_title=row_title,
                                column_title=row_title,
                                row_title_gp = gpar(fontsize = 10),
                                row_title_rot=0,
                                column_title_gp = gpar(fontsize = 10),
                                column_title_rot =90,
                                column_gap = unit(ifelse(border,c(0,0),c(2, 4)), "mm"),
                                row_gap = unit(ifelse(border,c(0,0),c(2, 4)), "mm")
                                )
  }else{
    hmp = ComplexHeatmap::Heatmap(rho_ordered, 
                                  use_raster = TRUE,
                                  border = border,
                                  col = col_fun,
                                  #row_title_gp = grid::gpar(col = split_colors), 
                                  #column_title_gp = grid::gpar(col = split_colors), 
                                  
                                  cluster_rows=F,
                                  cluster_columns=F,
                                  top_annotation = column_annot,
                                  left_annotation = row_annot,
                                  cluster_row_slices = F, 
                                  cluster_column_slices = F,
                                  show_row_names = FALSE, 
                                  show_column_names = FALSE,
                                  show_column_dend = show_dendrogram,
                                  show_row_dend = show_dendrogram,
                                  name=title,
                                  row_title=row_title,
                                  column_title=row_title,
                                  row_title_gp = gpar(fontsize = 10),
                                  row_title_rot=0,
                                  column_title_gp = gpar(fontsize = 10),
                                  column_title_rot =90,
                                  column_gap = unit(ifelse(border,c(0,0),c(2, 4)), "mm"),
                                  row_gap = unit(ifelse(border,c(0,0),c(2, 4)), "mm")
                                  )
  }
  return(list(heatmap = hmp, dend = gene_dendogram))
}

plotSplitCorrelationDualDiffHeatmap <- function(rho1,
					                                  rho2,
                                            name1=NULL,
                                            name2=NULL,
                                            show_dendrogram=F,
					                                  separate_plots,
					                                  gene_labels_1=NULL,
					                                  gene_labels_2=NULL,
					                                  split_list_1=NULL, 
					                                  split_list_2=NULL, 
					                                  split_colors_1=NULL, 
					                                  split_colors_2=NULL, 
                                            border=F,beta=2) { 
  if (separate_plots){

    layout(matrix(1:3,nrow=1))
    hm_1 = plotSplitCorrelationHeatmap_internal(rho=rho1,name=name1,difference=F,show_dendrogram=show_dendrogram,show_triangle="upper",gene_labels=gene_labels_1,split_list = split_list_1,split_colors=split_colors_1,border=border,beta=beta)
    dend_1 = hm_1[["dend"]]

    hm_2 = plotSplitCorrelationHeatmap_internal(rho=rho2,name=name2,difference=F,dend_spec = dend_1,show_dendrogram=show_dendrogram,show_triangle="lower",gene_labels=gene_labels_1,split_list = split_list_1,split_colors=split_colors_1,border=border,beta=beta)
    rho_diff = rho1-rho2
    hm_diff_1 = plotSplitCorrelationHeatmap_internal(rho=rho_diff,name=paste(name1,"-",name2),difference=T,dend_spec = dend_1,show_dendrogram=show_dendrogram,split_list = split_list_1,split_colors=split_colors_1,border=borde,beta=betar)
    return(hm_1[["heatmap"]]+hm_2[["heatmap"]]+hm_diff_1[["heatmap"]])
  }else{
    hm_1 = plotSplitCorrelationHeatmap_dual(rho_bl=rho1,rho_tr=rho2, cluster_on = "bl", title=paste("Fit on",name1),difference=F,name_bl=name1,name_tr=name2,
                                            show_dendrogram=show_dendrogram,gene_labels=gene_labels_1,split_list=split_list_1,split_colors=split_colors_1,border=border,beta=beta)
    hm_2 = plotSplitCorrelationHeatmap_dual(rho_bl=rho1,rho_tr=rho2, cluster_on = "tr", title=paste("Fit on",name2),difference=F,name_bl=name1,name_tr=name2,
                                            show_dendrogram=show_dendrogram,gene_labels=gene_labels_2,split_list=split_list_2,split_colors=split_colors_2,border=border,beta=beta)
    rho_diff_1 = rho1-rho2
    hm_diff_1 = plotSplitCorrelationHeatmap_internal(rho=rho_diff_1, name=paste(name1,"-",name2),difference=T,dend_spec = hm_1[["dend"]],
                                                     show_dendrogram=show_dendrogram,show_triangle="upper",gene_labels=gene_labels_1,split_list=split_list_1,split_colors=split_colors_1,border=border,beta=beta)
    hm_diff_2 = plotSplitCorrelationHeatmap_internal(rho=-rho_diff_1,name=paste(name2,"-",name1),difference=T,dend_spec = hm_2[["dend"]],
                                                     show_dendrogram=show_dendrogram,show_triangle="lower",gene_labels=gene_labels_2,split_list=split_list_2,split_colors=split_colors_2,border=border,beta=beta)
 
    return(list(rho1=hm_1[["heatmap"]],rho2=hm_2[["heatmap"]],rho1_diff=hm_diff_1[["heatmap"]],rho2_diff=hm_diff_2[["heatmap"]]))
    }
}
