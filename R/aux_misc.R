################################################################
###   CenTFinder                                             ###
###   Auxiliary function definitions: Miscellaneous          ###
###   Author: Jonathan De Smedt                              ###
################################################################

#' @include classes.R

arrayDistribution<-function(CenTFinderObject){
  if(all(dim(CenTFinderObject@platforms)==0)){
    stop("No microarray raw data was provided.")
  }
  if(any(names(CenTFinderObject@platforms)=="Used")){
    arrays<-CenTFinderObject@sample_annotation$Platform[CenTFinderObject@sample_annotation$Platform%in%CenTFinderObject@platforms$Platform[CenTFinderObject@platforms$Used]]
    title<-"Arrays per platform (filtered)"
  } else {
    arrays<-CenTFinderObject@sample_annotation$Platform
    title<-"Arrays per platform (unfiltered)"
  }
  arrays<-data.frame(table(arrays),stringsAsFactors = F)
  arrays$arrays<-factor(arrays$arrays,levels = arrays$arrays[order(arrays$Freq,decreasing = F)])
  ggplot2::ggplot(arrays,aes(x=arrays,y=Freq))+
    geom_bar(stat="identity",fill="#000000")+
    theme_bw()+
    theme(axis.text.x = element_text(size = 15, face="bold"),
          axis.text.y = element_text(size = 15, face="bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          plot.title = element_text(size=20, face = "bold"))+
    xlab("")+
    ylab("Number of arrays")+
    coord_flip()+
    ggtitle(title)
}

clusterDendrogram<-function(CenTFinderObject){
  if(length(CenTFinderObject@module_colors)==0){
    stop("No WGCNA analysis has been performed yet. Please run both applyWGCNA() and cutDendrogram().")
  }
  WGCNA::plotDendroAndColors(CenTFinderObject@tree, CenTFinderObject@module_colors,
                             "Module colors",rowText = CenTFinderObject@module_colors,
                             dendroLabels = FALSE, hang = 0.03,
                             addGuide = TRUE, guideHang = 0.05)
}


groupDistribution<-function(CenTFinderObject){
  if(all(dim(CenTFinderObject@sample_annotation)==0)){
    stop("No sample annotation has been loaded yet. Please run setSampleAnnotations().")
  }
  if(all(dim(CenTFinderObject@platforms)==0)&all(dim(CenTFinderObject@sample_annotation)==0)){
    stop("No platform information has been loaded yet. Please run getPlatformInfo().")
  }
  if(all(dim(CenTFinderObject@sample_annotation)!=0)&any(names(CenTFinderObject@platforms)=="Used")){
    groups<-CenTFinderObject@sample_annotation$Group[CenTFinderObject@sample_annotation$Platform%in%CenTFinderObject@platforms$Platform[CenTFinderObject@platforms$Used]]
  } else {
    groups<-CenTFinderObject@sample_annotation$Group
  }
  groups<-gsub(groups, pattern = "[.]",replacement = " ")
  groups<-data.frame(table(groups),stringsAsFactors = F)
  groups$groups<-factor(groups$groups,levels = groups$groups[order(groups$Freq,decreasing = F)])
  ggplot2::ggplot(groups,aes(x=groups,y=Freq))+
    geom_bar(stat="identity",fill="#000000")+
    theme_bw()+
    theme(axis.text.x = element_text(size = 15, face="bold"),
          axis.text.y = element_text(size = 10, face="bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          plot.title = element_text(size=20, face = "bold"))+
    xlab("")+
    ylab("Number of arrays")+
    coord_flip()+
    ggtitle("Arrays per group")
}

is.CenTFinder<-function(x){
  inherits(x,"CenTFinder")
}

meanConnectivityPlot<-function(CenTFinderObject){
  if(is.null(CenTFinderObject@plots$Mean_connectivity)){
    stop("No WGCNA analysis has been performed yet. Please run applyWGCNA().")
  }
  return(CenTFinderObject@plots$Mean_connectivity)
}

modulePCA<-function(CenTFinderObject, module=NULL, highlight=NULL){
  if(all(dim(CenTFinderObject@data)==0)){
    stop("No expression data has been loaded yet.")
  }
  if(all(is.na(CenTFinderObject@gene_cluster_annotation$Cluster_color))){
    stop("No WGCNA analysis has been performed yet on the data.")
  }
  if(is.null(module)){
    stop("Please provide a module of interest. Check available modules with modules().")
  }
  if(!any(names(CenTFinderObject@GO)%in%module)){
    message(paste0("The ",module," module is not one of the identified modules. Identified modules are: ",paste(names(CenTFinderObject@GO),collapse=", ")))
  } else {
    data<-CenTFinderObject@data[rownames(CenTFinderObject@data)%in%CenTFinderObject@gene_cluster_annotation$Gene[CenTFinderObject@gene_cluster_annotation$Cluster_color==module],]
    pcadata<-prcomp(t(data), scale. = T)
    dat.sum<-summary(pcadata)
    dat<-as.data.frame(pcadata$x)
    dat$Sample<-rownames(dat)
    dat<-merge(dat,CenTFinderObject@sample_annotation[,c("Sample","Group")], by="Sample",all=F)
    p<-ggplot2::ggplot(data = dat, aes(x = PC1, y = PC2, group=Group, col=Group)) +
      geom_point()+
      theme_bw()+
      geom_hline(yintercept = 0, colour = "gray65") +
      geom_vline(xintercept = 0, colour = "gray65")+
      xlab(paste("PC1",round(dat.sum$importance["Proportion of Variance",1]*100,2),"%",sep=", "))+
      ylab(paste("PC2",round(dat.sum$importance["Proportion of Variance",2]*100,2),"%",sep=", "))+
      theme(text = element_text(size=20),
            axis.text.x = element_text(angle=0, hjust=1))+
      ggtitle("PCA plot")+
      theme(legend.position="bottom",
            legend.text = element_text(size=5))+
      guides(fill=guide_legend(nrow=5, byrow=TRUE))
    if(!is.null(highlight)){
      dat2<-dat[dat$Group%in%highlight,]
      p<-p+geom_point(data=dat2,
                      pch=21, fill=NA, size=4, colour="red", stroke=1)
    }
    p
  }
}

plotGO<-function(CenTFinderObject, module, GO=c("BP","CC","MF"), GOID=F){
  if(length(CenTFinderObject@GO)==0){
    stop("No downstream analysis has been performed. Please run analyseClusters().")
  }
  dat<-CenTFinderObject@GO[[module]][[GO]]
  dat$classicFisher<-gsub(dat$classicFisher, pattern = "<",replacement = "")
  dat$classicFisher<-as.numeric(dat$classicFisher)
  if(GOID){
    dat$Term<-paste(dat$GO.ID, dat$Term, sep = " - ")
  } else {
    if(any(duplicated(dat$Term))){
      dat<-dat[!duplicated(dat$Term),]
    }
  }
  dat$Term<-factor(dat$Term, levels=dat$Term[order(dat$classicFisher,decreasing = T)])
  dat$logp<--log10(dat$classicFisher)
  title<-paste0("GO ",GO," of the ",module," module")
  ggplot2::ggplot(dat, aes(x=Term, y=logp))+
    theme_classic()+
    geom_bar(stat="identity",fill="#000000")+
    xlab("")+
    ylab("Fisher's exact test -log10(p value)")+
    scale_x_discrete(label = function(x) stringr::str_trunc(x, 35))+
    coord_flip()+
    ggtitle(title)+
    theme(plot.title = element_text(face="bold", size=10),
          axis.title.x = element_text(size=10))
}

plotExpression<-function(CenTFinderObject, genes, celltypes = "all", units){
  if(all(dim(CenTFinderObject@data)==0)){
    stop("No expression data has been loaded.")
  }
  if(all(dim(CenTFinderObject@sample_annotation)==0)){
    stop("No sample annotation has been loaded yet. Please run setSampleAnnotations().")
  }
  dat<-data.frame(CenTFinderObject@data,stringsAsFactors = F)
  if(celltypes!="all"){
    samples<-CenTFinderObject@sample_annotation$Sample[CenTFinderObject@sample_annotation$Group%in%celltypes]
    dat<-dat[,names(dat)%in%samples]
  }
  dat$Gene<-rownames(dat)
  dat<-dat[dat$Gene%in%genes,]
  dat<-gather(dat, key = "Sample", value = "Expression",1:(ncol(dat)-1))
  annot<-CenTFinderObject@sample_annotation[,c("Sample","Group")]
  annot$Sample<-gsub(annot$Sample,pattern = " |[+()/-]",replacement = ".")
  dat<-merge(dat,annot, by="Sample",all=F)
  if(length(genes)==1){
    title<-paste0("Expression profile of ",genes)
  } else {
    title<-paste0("Expression profiles of ",paste(genes[-length(genes)],collapse = ", "),", and ",genes[length(genes)])
  }
  ggplot2::ggplot(dat, aes(x = Group, y = Expression, group=interaction(Group,Gene), fill=Gene))+
    geom_boxplot(position = "dodge2")+
    theme_bw()+
    scale_y_continuous(breaks = seq(0,15,1), limits = c(0,15))+
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 40))+
    theme(axis.text.x = element_text(angle=60, hjust=1, size = 8))+
    xlab("")+
    ylab(units)+
    ggtitle(title)
}

plotGSVA<-function(CenTFinderObject,group, module_annot, labels , colors = c("#666666","#bfbfbf")){
  if(is.null(CenTFinderObject@integration$Uncompacted_GSVA)){
    stop("GSVA analysis has not yet been performed on this object.")
  }
  # Extract the GSVA results
  gsva_res<-CenTFinderObject@integration$Uncompacted_GSVA
  # Extract sample annotation
  annot<-CenTFinderObject@sample_annotation[,c("Sample","Group")]
  # Merge
  dat<-merge(annot,gsva_res, by.x="Sample",by.y="Cell")
  # Relabel groups
  dat$Group[dat$Group==group]<-labels[1]
  dat$Group[dat$Group!=labels[1]]<-labels[2]
  # Relabel modules
  for(module in unique(dat$Module)){
    dat$Module[dat$Module==module]<-module_annot[[module]]
  }
  # Generate plot
  dat$Group<-factor(dat$Group, levels = labels)
  ggplot2::ggplot(dat, aes(Module, gsva, fill=Group))+
    geom_hline(yintercept = 0, color="#FF0000", size=2)+
    geom_boxplot()+
    scale_fill_manual(values = colors, labels=labels)+
    theme_classic()+
    theme(axis.text.x = element_text(hjust=1, face="bold",size=16, angle = 70,color = "#000000"),
          axis.text.y = element_text(face="bold",size=16, color = "#000000"),
          axis.title = element_text(face = "bold", color = "#000000", size=20),
          legend.text = element_text(face = "bold", color = "#000000", size=16),
          legend.title = element_blank())+
    xlab("")+
    ylab("GSVA enrichment score")+
    coord_cartesian(xlim = c(1,length(unique(dat$Module))+1))+
    annotate(geom = "text", x = length(unique(dat$Module))+0.75, y= 0.5, label="Higher activity", color="#FF0000", angle=90, fontface="bold", size = 5)+
    annotate(geom = "text", x = length(unique(dat$Module))+0.75, y= -0.5, label="Lower activity", color="#FF0000", angle=90, fontface="bold", size = 5)
}

plotIntegration<-function(CenTFinderObject, minFC=NULL, maxFC=NULL, minKME=0.5, top=15, module_annot, modules_keep=names(module_annot), celltype=NULL, add_dataset_to_plot=NULL, topTFs_to_highlight=0){
  if(all(dim(CenTFinderObject@motif_enrichment)==0)){
    stop("No RcisTarget analysis has been performed yet.")
  }
  if(length(CenTFinderObject@additional_DE)==0){
    if(length(CenTFinderObject@DE)==0){
      stop("No differential expression analysis had been performed.")
    }
    warning("No additional differential expression data was provided. Output will make use of fold changes from the coexpression dataset.")
    if(is.null(celltype)){
      stop("Please provide a celltype for which to plot the fold changes.")
    }
  }
  # Extract the RcisTarget results
  rct_res<-CenTFinderObject@motif_enrichment
  # Extract the high-confidence TF from this frame
  rct_res$TF<-gsub(sapply(strsplit(rct_res$TF_highConf, split = " "), FUN = function(x){x[1]}), pattern = ";", replacement = "")
  # Remove NA values
  rct_res<-rct_res[!is.na(rct_res$TF),]
  # Select cols
  rct_res<-rct_res[,c("NES","TF","geneSet","enrichedGenes")]
  # Select only module enrichments
  rct_res<-rct_res[rct_res$geneSet!="markers",]
  # Reformat
  rct_new<-rct_res%>%
    group_by(geneSet,TF)%>%
    summarise(maxNES = max(NES),
              NrTarget = length(unique(strsplit(paste(enrichedGenes,collapse=";"),split = ";")[[1]])))
  # Merge with genes_per_cluster
  genes_per_cluster<-CenTFinderObject@gene_cluster_annotation%>%
    group_by(Cluster_color)%>%
    summarise(NrGenes = length(Gene))
  rct_new<-merge(rct_new,genes_per_cluster,by.x="geneSet",by.y="Cluster_color",all.x=T, all.y=F)
  rct_new$Fraction<-rct_new$NrTarget/rct_new$NrGenes
  # Relabel modules
  rct_new$Color<-rct_new$geneSet
  for(module in unique(rct_new$geneSet)){
    rct_new$geneSet[rct_new$geneSet==module]<-module_annot[[module]]
  }
  rct_new$geneSet<-gsub(pattern = " - ", replacement = "\n", rct_new$geneSet)
  # Select labels to plot
  rct_new<-rct_new%>%
    group_by(geneSet)%>%
    mutate(Plot_label = Fraction>sort(Fraction, decreasing = T)[top+1])
  # Generate a plot
  rct_new$geneSet<-paste0(rct_new$geneSet, "\nCluster size: ",rct_new$NrGenes," genes")
  # Merge with DE
  if(length(CenTFinderObject@additional_DE)==0){
    if(!is.null(celltype)){
      de<-CenTFinderObject@DE[[celltype]]
      de$logFC<-log2(de$FC)
    } else {
      de=NULL
    }
  } else {
    if(is.null(add_dataset_to_plot)){
      de<-CenTFinderObject@additional_DE[[1]]
    } else {
      de<-CenTFinderObject@additional_DE[[add_dataset_to_plot]]
    }
  }

  if(!is.null(de)){
    wrd<-merge(rct_new, de, by.x="TF",by.y="Gene",all.x=T, all.y=F)
    wrd<-wrd[!is.na(wrd$FC),]
    wrd<-wrd[wrd$Color!="grey",]
    kme<-CenTFinderObject@kMEs

    wrd$kME<-0

    for(i in 1:nrow(wrd)){
      wrd[i,"kME"]<-kme[wrd$TF[i],paste0("kME",wrd$Color[i])]
    }
    wrd<-wrd[!is.na(wrd$kME),]
    wrd$TF<-paste0(wrd$TF,"\nkME: ",round(wrd$kME,2),"\nNES: ",wrd$maxNES)
    wrd<-wrd[wrd$kME>minKME,]

    if(!is.null(minFC)){
      wrd<-wrd[wrd$FC>minFC,]
      if(nrow(wrd)==0){
        stop("No identified TFs are higher expressed than the cut-off in the provided additional data set.")
      }
    }
    if(!is.null(maxFC)){
      wrd<-wrd[wrd$FC<maxFC,]
      stop("No identified TFs are lower expressed than the cut-off in the provided additional data set.")
    }
    wrd$fontcolour="black"
    wrd<-wrd[,!is.na(names(wrd))]
    if(topTFs_to_highlight>0){
      wrd<-wrd%>%
        group_by(geneSet)%>%
        mutate(fontcolour = ifelse(TF%in%TF[order(Fraction, decreasing=T)][1:min(topTFs_to_highlight,length(TF))],"red","black"))
    }

    ggplot2::ggplot(wrd[wrd$Color%in%modules_keep,], aes(FC, Fraction))+
      geom_point(size=0.5, alpha=0.5)+
      theme_classic()+
      theme(axis.text.x = element_text(angle=70, hjust = 1, face="bold", color = "#000000"),
            axis.text.y = element_text(face = "bold", color = "#000000"),
            axis.title = element_text(face = "bold", color = "#000000"),
            strip.text = element_text(face = "bold", color = "#000000"),
            legend.position = "none")+
      scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))+
      scale_x_continuous(trans = "log2", breaks = c(1,10,100,1000,10000), limits = c(1,10000))+
      ylab("Fraction of module regulated")+
      xlab("Fold change")+
      facet_grid(.~geneSet)+
      ggrepel::geom_text_repel(aes(label=TF, color=fontcolour), force = 3, size=2, fontface="bold", segment.alpha = 0.3)+
      scale_color_manual(values = c("red"="#FF0000","black"="#000000"))
  } else {
    plotRcisTarget(CenTFinderObject = CenTFinderObject, module_annot = module_annot, top=top, topTFs_to_highlight = topTFs_to_highlight)
  }
}

plotKEGG<-function(CenTFinderObject, module){
  if(length(CenTFinderObject@KEGG)==0){
    stop("No downstream analysis has been performed. Please run analyseClusters().")
  }
  dat<-CenTFinderObject@KEGG[[module]]
  if(nrow(dat)==0){
    message("No KEGG pathways were enriched in this module.")
  } else {
    dat$Pathway<-paste(dat$ID, dat$Description, sep = "-")
    dat$qvalue<-as.numeric(dat$qvalue)
    dat$Pathway<-factor(dat$Pathway, levels = dat$Pathway[order(dat$qvalue,decreasing = T)])
    dat$logq<--log10(dat$qvalue)
    dat<-dat[order(dat$logq,decreasing = T),]
    dat<-dat[1:min(nrow(dat),15),]
    title<-paste0("KEGG enrichment of the ",module," module.")
    ggplot2::ggplot(dat, aes(x=Pathway, y=logq))+
      theme_bw()+
      geom_bar(stat="identity",fill="#000000")+
      xlab("")+
      ylab("Hypergeometric test -log10(q value)")+
      coord_flip()+
      ggtitle(title)
  }
}

plotMarkers<-function(CenTFinderObject, celltype, n, labels=NULL, units=c("log2int","fc","logfc"),sig=0.05, colours=c("#666666","#bfbfbf"), minimal=T, rotate=T){
  if(minimal){
    labels<-NULL
  }
  # Extract genes that are significantly differentially expressed
  genes<-CenTFinderObject@DE[[celltype]]
  genes<-genes[genes$AdjPvalue<sig,]
  genes<-genes[order(genes$FC, decreasing = T),][1:n,]
  genes<-rownames(genes)
  # Extract expression data of the selected genes
  expr<-CenTFinderObject@data[rownames(CenTFinderObject@data)%in%genes,]
  # Extract annotation information
  annot<-CenTFinderObject@sample_annotation[,c("Sample","Group")]
  # Transform expression data
  expr<-t(expr)
  expr<-data.frame(expr, stringsAsFactors = F)
  expr$Sample<-rownames(expr)
  expr<-merge(annot,expr,by="Sample")
  expr<-gather(expr, key = "Gene", value = "Expression", 3:ncol(expr))
  # Whichever is not the intended cell type is set to 'other'.
  expr$Sample<-NULL
  if(all(!is.null(labels))){
    expr$Group[expr$Group!=celltype]<-labels[2]
    expr$Group[expr$Group==celltype]<-labels[1]
  }
  # Calculate fold changes
  fc<-expr
  fc<-fc%>%
    group_by(Gene)%>%
    mutate(FoldChange = 2^(Expression - mean(Expression[Group!=ifelse(is.null(labels[1]),celltype,labels[1])])),
           LogFC = log2(FoldChange))
  gene_levels<-fc%>%
    filter(Group==ifelse(is.null(labels[1]),celltype,labels[1]))%>%
    group_by(Gene)%>%
    summarise(meanlogfc = mean(LogFC))
  gene_levels<-as.character(gene_levels$Gene[order(gene_levels$meanlogfc, decreasing = T)])
  fc$Gene<-factor(fc$Gene, levels = if(rotate){rev(gene_levels)}else{gene_levels})
  # Generate plot

  if(minimal){
    expr<-expr[expr$Group==celltype,]
    fc<-fc[fc$Group==celltype,]
  } else {
    expr$Group<-factor(expr$Group, levels = labels)
    fc$Group<-factor(fc$Group, levels = labels)
  }
  res<-ggplot2::ggplot(fc, aes(Gene, switch(units, "log2int"=Expression, "fc"=FoldChange, "log2fc"=LogFC), fill=Group))+
    scale_fill_manual(values = colours)+
    theme_classic()+
    theme(axis.text.x = if(rotate){element_text(face="bold",size=12, color = "#000000")}else{element_text(hjust=1, face="bold.italic",size=12, angle = 70, color = "#000000")},
          axis.text.y = if(rotate){element_text(hjust=1, face="bold.italic",size=12, color = "#000000")}else{element_text(face="bold",size=12, color = "#000000")},
          legend.title = element_blank(),
          axis.title = element_text(face="bold", color = "#000000", size=14),
          title = element_text(face = "bold", size=14))+
    scale_y_continuous(trans=switch(units, "log2int"="identity", "fc"="log10", "log2fc"="identity"))+
    xlab("")+
    ylab(switch(units, "log2int"="Log2(intensity)", "fc"="Fold change", "log2fc"="Log2(Fold change)"))+
    ggtitle(paste0("Putative ", gsub(celltype,pattern="[.]", replacement = " ")," markers"))

  if(minimal){
    res<-res+
      stat_summary(fun = mean,
                   fun.min = function(x) quantile(x,0.25),
                   fun.max = function(x) quantile(x,0.75),
                   geom = "pointrange",
                   fatten=3)+
      theme(legend.position = "none")
    } else {
    res<-res+
      geom_boxplot(outlier.size = 1)
  }
  if(rotate){
    res<-res+
      coord_flip()
  }
  res
}

# plotMarkers(lsec,"LSEC",40,units = "fc")
# plotMarkers(lsec,"LSEC",40,labels=c("LSECs","Other ECs"),units = "fc", colours = c("#FF0000","#00FF00"), minimal = F)
# plotMarkers(lsec,"LSEC",40,labels=c("LSECs","Other ECs"),units = "fc", colours = c("#FF0000","#00FF00"), minimal = T)
# plotMarkers(lsec,"LSEC",40,labels=c("LSECs","Other ECs"),units = "fc", colours = c("#FF0000","#00FF00"), minimal = F, rotate = F)


plotModuleScores<-function(CenTFinderObject, module_annot=NULL, traces=c("column","row","both","none")){
  if(length(CenTFinderObject@integration)==0){
    message("Please run integrateAnalyses in order to get module scores.")
  } else {
    mat<-as.matrix(CenTFinderObject@integration$Module_scores)
    if(!is.null(module_annot)){
      for(i in 1:ncol(mat)){
        colnames(mat)[i]<-module_annot[[colnames(mat)[i]]]
      }
    }
    gplots::heatmap.2(mat, col=gplots::bluered, margins = c(15,20), cexRow = 0.5, trace = traces, key.title = "Module activity")
  }
}



plotOverlapWithOtherDatasets<-function(CenTFinderObject, x_axis, y_axis, size=NULL, colour=NULL, titles=NULL, sig=0.05, FC_cutoff=2, label_FC_cutoff=7){
  if(missing(x_axis)){
    stop("Please provide a data set for which to plot the fold changes on the x-axis.")
  }
  if(missing(y_axis)){
    stop("Please provide a data set for which to plot the fold changes on the y-axis.")
  }

  # Extract the data sets
  if(substr(x_axis,1,3)=="DE:"){
    if(length(CenTFinderObject@DE)==0){
      stop("No differential expression was performed on the gene co-expression data set.")
    }
    if(!any(names(CenTFinderObject@DE)==substr(x_axis,4,nchar(x_axis)))){
      stop("x_axis data was not found. Please provide the exact group name. Group names can be extracted with the 'groups' function.")
    }
    x_axis<-CenTFinderObject@DE[[substr(x_axis,4,nchar(x_axis))]]
  } else {
    if(length(CenTFinderObject@additional_DE)==0){
      stop("No additional data sets were provided.")
    }
    if(!is.numeric(x_axis)){
      if(!any(names(CenTFinderObject@additional_DE)==x_axis)){
        stop(paste0("No additional data set has the name '",x_axis,"'."))
      }
    } else {
      if(length(CenTFinderObject@additional_DE)<x_axis){
        stop("Index of x_axis data set is out of bounds.")
      }
    }
    x_axis<-CenTFinderObject@additional_DE[[x_axis]]
  }
  if(substr(y_axis,1,3)=="DE:"){
    if(length(CenTFinderObject@DE)==0){
      stop("No differential expression was performed on the gene co-expression data set.")
    }
    if(!any(names(CenTFinderObject@DE)==substr(y_axis,4,nchar(y_axis)))){
      stop("y_axis data was not found. Please provide the exact group name. Group names can be extracted with the 'groups' function.")
    }
    y_axis<-CenTFinderObject@DE[[substr(y_axis,4,nchar(y_axis))]]
  } else {
    if(length(CenTFinderObject@additional_DE)==0){
      stop("No additional data sets were provided.")
    }
    if(!is.numeric(y_axis)){
      if(!any(names(CenTFinderObject@additional_DE)==y_axis)){
        stop(paste0("No additional data set has the name '",y_axis,"'."))
      }
    } else {
      if(length(CenTFinderObject@additional_DE)<y_axis){
        stop("Index of y_axis data set is out of bounds.")
      }
    }
    y_axis<-CenTFinderObject@additional_DE[[y_axis]]
  }
  if(!is.null(size)){
    if(substr(size,1,3)=="DE:"){
      if(length(CenTFinderObject@DE)==0){
        stop("No differential expression was performed on the gene co-expression data set.")
      }
      if(!any(names(CenTFinderObject@DE)==substr(size,4,nchar(size)))){
        stop("size data was not found. Please provide the exact group name. Group names can be extracted with the 'groups' function.")
      }
      size<-CenTFinderObject@DE[[substr(size,4,nchar(size))]]
    } else {
      if(length(CenTFinderObject@additional_DE)==0){
        stop("No additional data sets were provided.")
      }
      if(!is.numeric(size)){
        if(!any(names(CenTFinderObject@additional_DE)==size)){
          stop(paste0("No additional data set has the name '",size,"'."))
        }
      } else {
        if(length(CenTFinderObject@additional_DE)<size){
          stop("Index of size data set is out of bounds.")
        }
      }
      size<-CenTFinderObject@additional_DE[[size]]
    }
  }
  if(!is.null(colour)){
    if(substr(colour,1,3)=="DE:"){
      if(length(CenTFinderObject@DE)==0){
        stop("No differential expression was performed on the gene co-expression data set.")
      }
      if(!any(names(CenTFinderObject@DE)==substr(colour,4,nchar(colour)))){
        stop("colour data was not found. Please provide the exact group name. Group names can be extracted with the 'groups' function.")
      }
      colour<-CenTFinderObject@DE[[substr(colour,4,nchar(colour))]]
    } else {
      if(length(CenTFinderObject@additional_DE)==0){
        stop("No additional data sets were provided.")
      }
      if(!is.numeric(colour)){
        if(!any(names(CenTFinderObject@additional_DE)==colour)){
          stop(paste0("No additional data set has the name '",colour,"'."))
        }
      } else {
        if(length(CenTFinderObject@additional_DE)<colour){
          stop("Index of colour data set is out of bounds.")
        }
      }
      colour<-CenTFinderObject@additional_DE[[colour]]
    }
  }
  # Merge the datasets into one
  xy<-merge(x_axis, y_axis, by="Gene", all=F, suffixes = c("_x","_y"))

  if(!is.null(size)){
    names(size)[names(size)!="Gene"]<-paste0(names(size)[names(size)!="Gene"],"_size")
    xys<-merge(xy, size, by="Gene", all=F)
  }

  if(!is.null(colour)){
    names(colour)[names(colour)!="Gene"]<-paste0(names(colour)[names(colour)!="Gene"],"_colour")
    xysc<-merge(xys, colour, by="Gene", all=F)
  }

  # Keep only the significantly differentially expressed genes
  xysc<-xysc[xysc$AdjPvalue_x<sig&xysc$AdjPvalue_y<sig,]

  if(!is.null(size)){
    xysc<-xysc[xysc$AdjPvalue_size<sig,]
  }

  if(!is.null(colour)){
    xysc<-xysc[xysc$AdjPvalue_colour<sig,]
  }

  # Keep only the genes with a fold change higher than FC_cutoff
  xysc<-xysc[xysc$FC_x>FC_cutoff&xysc$FC_y>FC_cutoff,]

  if(!is.null(size)){
    xysc<-xysc[xysc$FC_size>FC_cutoff,]
  }

  if(!is.null(colour)){
    xysc<-xysc[xysc$FC_colour>FC_cutoff,]
  }

  # Plot
  if(is.null(size)&is.null(colour)){
    p<-ggplot(xysc, aes(FC_x,FC_y))+
      geom_point()+
      theme_classic()+
      scale_x_log10()+
      scale_y_log10()+
      xlab("scRNA MacParland et al. fold change")+
      ylab("scRNA Aizarani et al. fold change")+
      theme(axis.text = element_text(face = "bold"),
            axis.title = element_text(face = "bold"),
            legend.text = element_text(face = "bold"),
            legend.title = element_text(face = "bold"))+
      geom_text_repel(aes(label=ifelse(FC_x>label_FC_cutoff&FC_y>label_FC_cutoff,Gene,"")), force = 2, fontface = "bold.italic")
    if(!is.null(titles)){
      if(length(titles)==2){
        p<-p+
          xlab(titles[1])+
          ylab(titles[2])
      } else {
        warning("One x-axis title and one y-axis title need to be provided.")
      }
    }
  }
  if(is.null(size)&!is.null(colour)){
    p<-ggplot(xysc, aes(FC_x,FC_y, colour=FC_colour))+
      geom_point()+
      theme_classic()+
      scale_x_log10()+
      scale_y_log10()+
      xlab("scRNA MacParland et al. fold change")+
      ylab("scRNA Aizarani et al. fold change")+
      theme(axis.text = element_text(face = "bold"),
            axis.title = element_text(face = "bold"),
            legend.text = element_text(face = "bold"),
            legend.title = element_text(face = "bold"))+
      geom_text_repel(aes(label=ifelse(FC_x>label_FC_cutoff&FC_y>label_FC_cutoff&FC_colour>label_FC_cutoff,Gene,"")), force = 2, fontface = "bold.italic")
    if(!is.null(titles)){
      if(length(titles)==3){
        p<-p+
          xlab(titles[1])+
          ylab(titles[2])+
          labs(colour = titles[3])
      } else {
        warning("One x-axis title, one y-axis title, and one colour legend title need to be provided.")
      }
    }
  }
  if(!is.null(size)&is.null(colour)){
    p<-ggplot(xysc, aes(FC_x,FC_y, size=FC_size))+
      geom_point()+
      theme_classic()+
      scale_x_log10()+
      scale_y_log10()+
      xlab("scRNA MacParland et al. fold change")+
      ylab("scRNA Aizarani et al. fold change")+
      theme(axis.text = element_text(face = "bold"),
            axis.title = element_text(face = "bold"),
            legend.text = element_text(face = "bold"),
            legend.title = element_text(face = "bold"))+
      geom_text_repel(aes(label=ifelse(FC_x>label_FC_cutoff&FC_y>label_FC_cutoff&FC_size>label_FC_cutoff,Gene,"")), force = 2, fontface = "bold.italic")
    if(!is.null(titles)){
      if(length(titles)==3){
        p<-p+
          xlab(titles[1])+
          ylab(titles[2])+
          labs(size = titles[3])
      } else {
        warning("One x-axis title, one y-axis title, and one size legend title need to be provided.")
      }
    }
  }
  if(!is.null(size)&!is.null(colour)){
    p<-ggplot(xysc, aes(FC_x,FC_y, size=FC_size, colour=FC_colour))+
      geom_point()+
      theme_classic()+
      scale_x_log10()+
      scale_y_log10()+
      xlab("scRNA MacParland et al. fold change")+
      ylab("scRNA Aizarani et al. fold change")+
      theme(axis.text = element_text(face = "bold"),
            axis.title = element_text(face = "bold"),
            legend.text = element_text(face = "bold"),
            legend.title = element_text(face = "bold"))+
      geom_text_repel(aes(label=ifelse(FC_x>label_FC_cutoff&FC_y>label_FC_cutoff&FC_size>label_FC_cutoff&FC_colour>label_FC_cutoff,Gene,"")), force = 2, fontface = "bold.italic")
    if(!is.null(titles)){
      if(length(titles)==4){
        p<-p+
          xlab(titles[1])+
          ylab(titles[2])+
          labs(size = titles[3])+
          labs(colour = titles[4])
      } else {
        warning("One x-axis title, one y-axis title, one size legend title, and one colour legend title need to be provided.")
      }
    }
  }
  p
}



plotRcisTarget<-function(CenTFinderObject, module_annot, top, topTFs_to_highlight=0){
  if(all(dim(CenTFinderObject@motif_enrichment)==0)){
    stop("No RcisTarget analysis has been performed yet.")
  }
  # Extract the RcisTarget results
  rct_res<-CenTFinderObject@motif_enrichment
  # Extract the high-confidence TF from this frame
  rct_res$TF<-gsub(sapply(strsplit(rct_res$TF_highConf, split = " "), FUN = function(x){x[1]}), pattern = ";", replacement = "")
  # Remove NA values
  rct_res<-rct_res[!is.na(rct_res$TF),]
  # Select cols
  rct_res<-rct_res[,c("NES","TF","geneSet","enrichedGenes")]
  # Select only module enrichments
  rct_res<-rct_res[rct_res$geneSet!="markers",]
  # Reformat
  rct_new<-rct_res%>%
    group_by(geneSet,TF)%>%
    summarise(maxNES = max(NES),
              NrTarget = length(unique(strsplit(paste(enrichedGenes,collapse=";"),split = ";")[[1]])))
  # Merge with genes_per_cluster
  genes_per_cluster<-CenTFinderObject@gene_cluster_annotation%>%
    group_by(Cluster_color)%>%
    summarise(NrGenes = length(Gene))
  rct_new<-merge(rct_new,genes_per_cluster,by.x="geneSet",by.y="Cluster_color",all.x=T, all.y=F)
  rct_new$Fraction<-rct_new$NrTarget/rct_new$NrGenes
  # Relabel modules
  rct_new$Color<-rct_new$geneSet
  for(module in unique(rct_new$geneSet)){
    rct_new$geneSet[rct_new$geneSet==module]<-module_annot[[module]]
  }
  rct_new$geneSet<-gsub(pattern = " - ", replacement = "\n", rct_new$geneSet)
  # Select labels to plot
  rct_new<-rct_new%>%
    group_by(geneSet)%>%
    mutate(Plot_label = Fraction>sort(Fraction, decreasing = T)[top+1])
  # Generate a plot
  rct_new$geneSet<-paste0(rct_new$geneSet, "\nCluster size: ",rct_new$NrGenes," genes")

  rct_new<-rct_new[,!is.na(names(rct_new))]
  rct_new$fontcolour="black"
  if(topTFs_to_highlight>0){
    rct_new<-rct_new%>%
      group_by(geneSet, Plot_label)%>%
      mutate(fontcolour = ifelse(TF%in%TF[order(Fraction, decreasing=T)][1:min(topTFs_to_highlight,length(TF))],"red","black"))
  }

  p<-ggplot2::ggplot(rct_new[!grepl(rct_new$geneSet,pattern="Unclustered genes"),], aes(maxNES, Fraction))+
    geom_point(size=0.5, alpha=0.5)+
    theme_classic()+
    scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1))+
    scale_x_log10(limits=c(3,40), breaks = c(3,5,10,20), labels=c(3,5,10,20))+
    xlab("Maximal network enrichment score (NES)")+
    ylab("Fraction of module regulated")+
    theme(axis.text = element_text(face="bold", color = "#000000"),
          axis.title = element_text(face = "bold", color = "#000000"),
          strip.text = element_text(face = "bold", color = "#000000"),
          legend.position = "none")
  if(length(unique(rct_new$geneSet))>7){
    p<-p+
      facet_wrap(facets="geneSet", nrow=2)+
      ggrepel::geom_text_repel(aes(label=ifelse(Plot_label,TF,""),color=fontcolour), force = 3, size=2, fontface="bold.italic", segment.alpha = 0.3)+
      scale_color_manual(values = c("red"="#FF0000","black"="#000000"))
  } else {
    p<-p+
      facet_grid(.~geneSet)+
      ggrepel::geom_text_repel(aes(label=ifelse(Plot_label,TF,""), color=fontcolour), force = 3, size=2, fontface="bold.italic", segment.alpha = 0.3)+
      scale_color_manual(values = c("red"="#FF0000","black"="#000000"))
  }
  p
}

plotTFs<-function(CenTFinderObject, tfs, celltype, n, labels=NULL, units=c("log2int","fc","logfc"),sig=0.05, colours=c("#666666","#bfbfbf"), minimal=T, rotate=T){
  if(minimal){
    labels<-NULL
  }
  # Extract genes that are significantly differentially expressed
  genes<-CenTFinderObject@DE[[celltype]]
  genes<-genes[genes$AdjPvalue<sig,]
  genes<-genes[genes$Gene%in%tfs,]
  genes<-genes[order(genes$FC, decreasing = T),][1:n,]
  genes<-rownames(genes)
  # Extract expression data of the selected genes
  expr<-CenTFinderObject@data[rownames(CenTFinderObject@data)%in%genes,]
  # Extract annotation information
  annot<-CenTFinderObject@sample_annotation[,c("Sample","Group")]
  # Transform expression data
  expr<-t(expr)
  expr<-data.frame(expr, stringsAsFactors = F)
  expr$Sample<-rownames(expr)
  expr<-merge(annot,expr,by="Sample")
  expr<-gather(expr, key = "Gene", value = "Expression", 3:ncol(expr))
  # Whichever is not the intended cell type is set to 'other'.
  expr$Sample<-NULL
  if(all(!is.null(labels))){
    expr$Group[expr$Group!=celltype]<-labels[2]
    expr$Group[expr$Group==celltype]<-labels[1]
  }
  # Calculate fold changes
  fc<-expr
  fc<-fc%>%
    group_by(Gene)%>%
    mutate(FoldChange = 2^(Expression - mean(Expression[Group!=ifelse(is.null(labels[1]),celltype,labels[1])])),
           LogFC = log2(FoldChange))
  gene_levels<-fc%>%
    filter(Group==ifelse(is.null(labels[1]),celltype,labels[1]))%>%
    group_by(Gene)%>%
    summarise(meanlogfc = mean(LogFC))
  gene_levels<-as.character(gene_levels$Gene[order(gene_levels$meanlogfc, decreasing = T)])
  fc$Gene<-factor(fc$Gene, levels = if(rotate){rev(gene_levels)}else{gene_levels})
  # Generate plot

  if(minimal){
    expr<-expr[expr$Group==celltype,]
    fc<-fc[fc$Group==celltype,]
  } else {
    expr$Group<-factor(expr$Group, levels = labels)
    fc$Group<-factor(fc$Group, levels = labels)
  }
  res<-ggplot2::ggplot(fc, aes(Gene, switch(units, "log2int"=Expression, "fc"=FoldChange, "log2fc"=LogFC), fill=Group))+
    scale_fill_manual(values = colours)+
    theme_classic()+
    theme(axis.text.x = if(rotate){element_text(face="bold",size=12, color = "#000000")}else{element_text(hjust=1, face="bold.italic",size=12, angle = 70, color = "#000000")},
          axis.text.y = if(rotate){element_text(hjust=1, face="bold.italic",size=12, color = "#000000")}else{element_text(face="bold",size=12, color = "#000000")},
          legend.title = element_blank(),
          axis.title = element_text(face="bold", color = "#000000", size=14),
          title = element_text(face = "bold", size=14))+
    scale_y_continuous(trans=switch(units, "log2int"="identity", "fc"="log10", "log2fc"="identity"))+
    xlab("")+
    ylab(switch(units, "log2int"="Log2(intensity)", "fc"="Fold change", "log2fc"="Log2(Fold change)"))+
    ggtitle(paste0("Putative ", gsub(celltype,pattern="[.]", replacement = " "),"-specific transcription factors"))

  if(minimal){
    res<-res+
      stat_summary(fun = mean,
                   fun.min = function(x) quantile(x,0.25),
                   fun.max = function(x) quantile(x,0.75),
                   geom = "pointrange",
                   fatten=3)+
      theme(legend.position = "none")
  } else {
    res<-res+
      geom_boxplot(outlier.size = 1)
  }
  if(rotate){
    res<-res+
      coord_flip()
  }
  res
}




plotVariance<-function(CenTFinderObject, genes=NULL, top=15){
  if(all(dim(CenTFinderObject@data)==0)){
    stop("No expression data has been loaded.")
  }
  dat<-data.frame(cbind(rownames(CenTFinderObject@data),matrixStats::rowVars(CenTFinderObject@data)), stringsAsFactors = F)
  names(dat)<-c("Gene","Variance")
  if(is.null(genes)){
    dat<-dat[order(dat$Variance, decreasing = T),][c(1:top,(nrow(dat)-top+1):nrow(dat)),]
    dat$Variance<-as.numeric(dat$Variance)
    dat$Gene<-factor(dat$Gene, levels = dat$Gene[order(dat$Variance,decreasing = F)])
    ggplot2::ggplot(dat,aes(x=Gene, y=Variance))+
      theme_bw()+
      geom_bar(stat="identity",fill="#000080")+
      xlab("")+
      ylab("Variance")+
      coord_flip()+
      ggtitle(paste0("Top ",top," most and top ",top," least variable genes"))
  } else {
    dat<-dat[dat$Gene%in%genes,]
    dat$Variance<-as.numeric(dat$Variance)
    dat$Gene<-factor(dat$Gene, levels = dat$Gene[order(dat$Variance,decreasing = F)])
    ggplot2::ggplot(dat,aes(x=Gene, y=Variance))+
      theme_bw()+
      geom_bar(stat="identity",fill="#000080")+
      xlab("")+
      ylab("Variance")+
      coord_flip()+
      ggtitle("Variance of selected genes")
  }
}


scaleIndependencePlot<-function(CenTFinderObject){
  if(is.null(CenTFinderObject@plots$Scale_independence)){
    stop("No WGCNA analysis has been performed yet. Please run applyWGCNA().")
  }
  return(CenTFinderObject@plots$Scale_independence)
}


tree<-function(CenTFinderObject){
  if(length(CenTFinderObject@module_colors)==0){
    stop("No WGCNA analysis has been performed yet. Please run applyWGCNA().")
  }
  plot(CenTFinderObject@tree)
}
