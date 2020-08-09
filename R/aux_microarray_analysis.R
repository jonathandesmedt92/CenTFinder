###########################################################################################
###   CenTFinder                                                                        ###
###   Auxiliary function definitions: Microarray analysis and batch correction          ###
###   Author: Jonathan De Smedt                                                         ###
###########################################################################################

# Include
#' @include classes.R

get_microarrays_info<-function(path, platform_freq_cutoff=0.00258, filenames){
  if(is.null(path)){
    stop("Please provide a path.")
  }
  # Get file names and extract platform information
  files <- affy::list.celfiles(path=path)
  filenames<-files[files%in%filenames]
  f<-paste(path,filenames,sep = "")
  f<-data.frame(sapply(f, function(x) affyio::read.celfile.header(x)$cdfName),stringsAsFactors = F)
  f$File_name<-rownames(f)
  names(f)<-c("Platform","Path")
  f$File<-filenames
  f<-f%>%
    group_by(Platform)%>%
    mutate(Frequency = length(File))
  # Remove files belonging to infrequent platforms (<0.05%)
  f<-f[f$Frequency>ceiling(nrow(f)*platform_freq_cutoff),]
  # Return
  return(f)
}

annotate_platform<-function(CenTFinderObject,platform){
  if(platform%in%c("HG-U133_Plus_2","NuGO_Hs1a520180","HG-U219","HT_HG-U133A")){
    db<-switch(platform,
               "HG-U133_Plus_2" = list(hgu133plus2.db::hgu133plus2ACCNUM,hgu133plus2.db::hgu133plus2SYMBOL,hgu133plus2.db::hgu133plus2GENENAME),
               "NuGO_Hs1a520180"=list(nugohs1a520180.db::nugohs1a520180ACCNUM,nugohs1a520180.db::nugohs1a520180SYMBOL,nugohs1a520180.db::nugohs1a520180SYMBOL),
               "HG-U219"=list(hgu219.db::hgu219ACCNUM,hgu219.db::hgu219SYMBOL,hgu219.db::hgu219SYMBOL),
               "HT_HG-U133A"=list(hgu133a.db::hgu133aACCNUM,hgu133a.db::hgu133aSYMBOL,hgu133a.db::hgu133aSYMBOL)
    )
    annot<-data.frame(SYMBOL=sapply(contents(db[[2]]), paste, collapse=", "),
                      ACCNUM =0)
    annot$ACCNUM<-rownames(annot)
  } else if(platform%in%names(CenTFinderObject@probe_annotations)){
    annot<-CenTFinderObject@probe_annotations[[platform]]
  } else {
    annot<-NULL
  }
  return(annot)
}

gm_mean<-function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

analyse_microarray<-function(CenTFinderObject,arrays, platform){
  # Read in the arrays
  eset<-oligo::read.celfiles(unique(arrays))
  # RMA normalise
  affyRaw<-oligo::rma(eset)
  # Extract expression matrix
  data<-data.frame(Biobase::exprs(affyRaw),check.names = F)
  # Create probe annotation
  annot<-annotate_platform(CenTFinderObject,platform)
  # Merge annotation and data
  data <- merge(annot, data, by.x="ACCNUM", by.y=0, all=F)
  #Collapse probes to genes by means of geometric mean
  data<-data[data$SYMBOL!="NA",]
  data<-data[data$SYMBOL!="---",]
  if(any(grepl(data$SYMBOL, pattern = "//"))){
    data$SYMBOL<-trimws(sapply(strsplit(data$SYMBOL,split="//"),FUN=function(x){x[2]}))
  }
  data<-aggregate(.~data$SYMBOL,data=data[,-c(1:2)],FUN=gm_mean,na.rm=T)
  rownames(data)<-data$`data$SYMBOL`
  data$`data$SYMBOL`<-NULL
  return(data)
}

calculate_cumulative_gene_consensus<-function(ordered_platforms, data){
  gene_count<-NULL
  for(i in seq_along(ordered_platforms)){
    if(i==1){
      gene_count<-append(gene_count,nrow(data[[ordered_platforms[i]]]))
      genes<-rownames(data[[ordered_platforms[i]]])
    } else {
      gene_count<-append(gene_count,sum(genes%in%row.names(data[[ordered_platforms[i]]])))
      genes<-genes[genes%in%rownames(data[[ordered_platforms[i]]])]
    }
  }
  return(gene_count)
}

merge_expr<-function(x,y){
  res<-merge(x,y,by=0,all=F)
  rownames(res)<-res$Row.names
  res$Row.names<-NULL
  return(res)
}

mp_analyse_microarray<-function(CenTFinderObject, platform_keep=NULL){
  # Loop over the different platforms and analyse per platform
  platforms<-CenTFinderObject@platforms[CenTFinderObject@platforms$Annotation!="No annotation found.","Platform"]
  annotation<-CenTFinderObject@sample_annotation
  data<-list()
  for(i in seq_along(platforms)){
    data[[i]]<-tryCatch({analyse_microarray(arrays = as.character(as.matrix(annotation[annotation$Platform==platforms[i]&!is.na(annotation$Group),"Path"])),
                                  platform = platforms[i],
                                  CenTFinderObject = CenTFinderObject)},
                        error = function(e){
                          message(paste0("Could not analyse arrays for platform ",platforms[i],". Continuing without..."))
                          return(NA)
                        })
  }
  names(data)<-platforms
  data<-data[!is.logical(data)]
  # Weigh number of genes and number of arrays
  weighing<-data.frame(genes = sapply(data,nrow),
                       samples = sapply(data,ncol),
                       platforms = platforms,
                       stringsAsFactors = F)
  weighing<-weighing[order(weighing$genes, decreasing = T),]
  row.names(weighing)<-NULL
  weighing$samples_cumsum<-cumsum(weighing$samples)
  weighing$genecount_cumsum<-calculate_cumulative_gene_consensus(ordered_platforms = weighing$platforms, data=data)
  weighing$number_of_cells<-weighing$samples_cumsum*weighing$genecount_cumsum
  if(is.null(platform_keep)){
    keep_platforms<-weighing$platforms[1:which.max(weighing$number_of_cells)]
  } else {
    keep_platforms<-weighing$platforms[1:max(which.max(weighing$number_of_cells),which(weighing$platforms==platform_keep))]
  }
  # Remove data and change platform info
  CenTFinderObject@platforms$Used<-FALSE
  CenTFinderObject@platforms[CenTFinderObject@platforms$Platform%in%keep_platforms,"Used"]<-TRUE
  data<-data[keep_platforms]
  # Merge data
  data<-Reduce(function(x,y) merge_expr(x=x, y=y),data)
  # Perform prior quantile normalisation
  data<-limma::normalizeQuantiles(data)
  # Adjust batch effects with Combat
  if(length(keep_platforms)>1){
    modmat<-as.matrix(CenTFinderObject@sample_annotation[order(match(names(data),CenTFinderObject@sample_annotation$File),decreasing = F),])
    batch<-as.factor(modmat[,"Platform"])
    data<-sva::ComBat(dat=as.matrix(data), batch=batch, mod=NULL)
  }
  # Remove hyphens from gene names
  rownames(data)<-gsub(rownames(data), pattern = "-",replacement = "_")
  # Fill in the data
  CenTFinderObject@data<-as.matrix(data)
  # Return
  return(CenTFinderObject)
}

filter_genes<-function(data, tfs=humantfs, TF_range_cutoff=1, max_genes=9000){
  # We consider different cut-offs for TFs and non-TFs, so we split the data
  data_TFs<-data[rownames(data)%in%tfs,]
  data_noTFs<-data[!rownames(data)%in%tfs,]
  # For TFs, we require the range to be ideally more than 1. (Customisable with TF_range_cutoff)
  tfranges<-apply(data_TFs,1,FUN = function(x){max(x)-min(x)})
  data_TFs<-data_TFs[tfranges>TF_range_cutoff,]
  # To run things fast enough, we require the total number of genes to be max 9000. (Customisable with max_genes)
  # For non-TFs, we select the most variable genes.
  no_tf_vars<-apply(data_noTFs,1,FUN = function(x){var(x)})
  data_noTFs<-data_noTFs[order(no_tf_vars, decreasing = T),]
  data_noTFs<-data_noTFs[c(1:min(nrow(data_noTFs),max_genes-nrow(data_TFs))),]
  # Bind data
  data<-rbind(data_TFs,data_noTFs)
  genes<-unique(rownames(data))
  return(genes)
}

multipleTtest<-function(data, samples){
  data<-as.matrix(data)
  m1<-rowMeans(data[,colnames(data)%in%samples])
  m2<-rowMeans(data[,!colnames(data)%in%c(samples,"FC")])
  sd1<-matrixStats::rowMads(data[,colnames(data)%in%samples])
  sd2<-matrixStats::rowMads(data[,!colnames(data)%in%c(samples,"FC")])
  num1<-length(samples)
  num2<-ncol(data)-num1-1
  se <- sqrt(sd1*sd1/num1+sd2*sd2/num2)
  t <- (m1-m2)/se
  p<-pt(-abs(t),df=pmin(num1,num2)-1)*2
  return(p)
}

multiple_wilcoxon_rank_sum_test<-function(data, samples){
  data<-as.matrix(data)
  p<-NULL
  message("Calculating Wilcoxon rank sum tests ... ")
  pb <- txtProgressBar(char = '=', style = 3, file = stderr())
  for(i in 1:nrow(data)){
    p<-append(p,wilcox.test(x = as.numeric(data[i,colnames(data)%in%samples]), y = as.numeric(data[i,!colnames(data)%in%c(samples,"FC")]), paired = F)$p.value)
    setTxtProgressBar(pb = pb, value = i / nrow(data))
  }
  return(p)
}

calculate_OneSample_Pvalue<-function(data, sample){
  data<-as.matrix(data)
  m<-rowMeans(data[,!colnames(data)%in%c(sample,"FC")])
  sd<-matrixStats::rowMads(data[,!colnames(data)%in%c(sample,"FC")])
  z<-(as.numeric(data[,colnames(data)==sample])-m)/sd
  p<-pnorm(-abs(z))*2
  return(p)
}

microarrayDE_t_test<-function(CenTFinderObject, sig=0.05, pre_analysed=F, contrast=NULL){
  if(pre_analysed){
    groups<-unique(CenTFinderObject@sample_annotation$Group)
  } else {
    platforms<-CenTFinderObject@platforms[CenTFinderObject@platforms$Used,"Platform"]
    groups<-unique(as.character(as.matrix(CenTFinderObject@sample_annotation[CenTFinderObject@sample_annotation$Platform%in%platforms,"Group"])))
  }
  if(!all(contrast%in%groups)){
    stop("One or two of the contrast groups is absent or misspelled.")
  }
  if(!is.null(contrast)){
    groups<-contrast
  }
  DEs<-list()
  i=0
  for(group in groups){
    i=i+1
    message(paste0("Calculating differential expression of cell type ",i," of ",length(groups),"..."))
    samples<-as.character(as.matrix(CenTFinderObject@sample_annotation[CenTFinderObject@sample_annotation$Group==group,"Sample"]))
    samples<-gsub(samples,pattern = " |[+()/]",replacement = ".")
    expr<-data.frame(CenTFinderObject@data,stringsAsFactors = F)

    if(!is.null(contrast)){
      c_samples<-as.character(as.matrix(CenTFinderObject@sample_annotation[CenTFinderObject@sample_annotation$Group==groups[groups!=group],"Sample"]))
      c_samples<-gsub(c_samples,pattern = " |[+()/]",replacement = ".")
      expr<-expr[,names(expr)%in%append(samples,c_samples)]
    }

    if(length(samples)>1){
      expr$FC<-2^(as.numeric(rowSums(expr[,names(expr)%in%samples])/length(samples)-rowSums(expr[,!names(expr)%in%samples])/(ncol(expr)-length(samples))))
      expr$Pvalue<-multipleTtest(data=expr, samples = samples)
    } else {
      expr$FC<-2^(as.numeric(expr[,names(expr)==samples])-rowSums(expr[,!names(expr)%in%samples])/(ncol(expr)-1))
      expr$Pvalue<-calculate_OneSample_Pvalue(data=expr, sample = samples)
    }
    expr$AdjPvalue<-p.adjust(expr$Pvalue, method="BH")
    expr$Gene<-rownames(expr)
    DEs[[group]]<-expr[,c("Gene","FC","Pvalue","AdjPvalue")]
  }
  names(DEs)<-gsub(names(DEs),pattern = " |[+()/]",replacement = ".")

  return(DEs)
}

microarrayDE_Wilcoxon_rank_sum_test<-function(CenTFinderObject, sig=0.05, pre_analysed=F, contrast=NULL){
  if(pre_analysed){
    groups<-unique(CenTFinderObject@sample_annotation$Group)
  } else {
    platforms<-CenTFinderObject@platforms[CenTFinderObject@platforms$Used,"Platform"]
    groups<-unique(as.character(as.matrix(CenTFinderObject@sample_annotation[CenTFinderObject@sample_annotation$Platform%in%platforms,"Group"])))
  }
  if(!all(contrast%in%groups)){
    stop("One or two of the contrast groups is absent or misspelled.")
  }
  if(!is.null(contrast)){
    groups<-contrast
  }
  DEs<-list()
  i=0
  for(group in groups){
    i=i+1
    cat("\t\n")
    message(paste0("Calculating differential expression of cell type ",i," of ",length(groups),"..."))
    samples<-as.character(as.matrix(CenTFinderObject@sample_annotation[CenTFinderObject@sample_annotation$Group==group,"Sample"]))
    samples<-gsub(samples,pattern = " |[+()/]",replacement = ".")
    expr<-data.frame(CenTFinderObject@data,stringsAsFactors = F)

    if(!is.null(contrast)){
      c_samples<-as.character(as.matrix(CenTFinderObject@sample_annotation[CenTFinderObject@sample_annotation$Group==groups[groups!=group],"Sample"]))
      c_samples<-gsub(c_samples,pattern = " |[+()/]",replacement = ".")
      expr<-expr[,names(expr)%in%append(samples,c_samples)]
    }

    if(length(samples)>1){
      expr$FC<-2^(as.numeric(rowSums(expr[,names(expr)%in%samples])/length(samples)-rowSums(expr[,!names(expr)%in%samples])/(ncol(expr)-length(samples))))
      expr$Pvalue<-multiple_wilcoxon_rank_sum_test(data=expr,samples = samples)
    } else {
      expr$FC<-2^(as.numeric(expr[,names(expr)==samples])-rowSums(expr[,!names(expr)%in%samples])/(ncol(expr)-1))
      expr$Pvalue<-calculate_OneSample_Pvalue(data=expr, sample = samples)
    }
    expr$AdjPvalue<-p.adjust(expr$Pvalue, method="BH")
    expr$Gene<-rownames(expr)
    DEs[[group]]<-expr[,c("Gene","FC","Pvalue","AdjPvalue")]
  }
  names(DEs)<-gsub(names(DEs),pattern = " |[+()/]",replacement = ".")

  return(DEs)
}

calcZscore<-function(x){
  m<-mean(x)
  s<-sd(x)
  z<-(x-m)/s
  return(z)
}

calcZscoreOfMatrix<-function(expr){
  exprZ<-data.frame(t(apply(expr,MARGIN = 1, FUN=function(x){calcZscore(x)})),stringsAsFactors = F)
  return(exprZ)
}

microarrayDE_Z_Scores<-function(CenTFinderObject, pre_analysed=F){
  # Calculate Z scores
  Zmat<-calcZscoreOfMatrix(expr=data.frame(CenTFinderObject@data,stringsAsFactors = F))
  # Split per cell type
  if(pre_analysed){
    groups<-unique(CenTFinderObject@sample_annotation$Group)
  } else {
    platforms<-CenTFinderObject@platforms[CenTFinderObject@platforms$Used,"Platform"]
    groups<-unique(as.character(as.matrix(CenTFinderObject@sample_annotation[CenTFinderObject@sample_annotation$Platform%in%platforms,"Group"])))
  }
  groups<-groups[!is.na(groups)]
  Zs<-list()
  for(celltype in groups){
    samples<-as.character(as.matrix(CenTFinderObject@sample_annotation[CenTFinderObject@sample_annotation$Group==celltype,"Sample"]))
    samples<-gsub(samples,pattern = " |[+()/]",replacement = ".")

    if(length(samples)==1){
      Zs[[celltype]]<-data.frame(cbind(Zmat[,names(Zmat)%in%samples],Zmat[,names(Zmat)%in%samples]),stringsAsFactors = F)
      names(Zs[[celltype]])<-c(samples,"Avg_Z")
      rownames(Zs[[celltype]])<-Zs[[celltype]]$Gene<-rownames(Zmat)
    } else {
      Zs[[celltype]]<-Zmat[,names(Zmat)%in%samples]
      Zs[[celltype]]$Avg_Z<-rowSums(Zs[[celltype]])/ncol(Zs[[celltype]])
      rownames(Zs[[celltype]])<-Zs[[celltype]]$Gene<-rownames(Zmat)
    }
  }
  names(Zs)<-gsub(names(Zs),pattern = " |[+()/]",replacement = ".")
  return(Zs)
}




