################################################################
###   CenTFinder                                             ###
###   Auxiliary function definitions: core analyses          ###
###   Author: Jonathan De Smedt                              ###
################################################################

# Include
#' @include classes.R

calc_TOM<-function(data, network_type="signed",powers=1:30){
  # Transpose data
  data<-data.frame(t(data),stringsAsFactors = F)
  # Picking soft tresholding power
  message("Picking soft thresholding power...")
  sft <- WGCNA::pickSoftThreshold(data, powerVector = powers, networkType = network_type)
  # Create plots
  p<-list()
  p[["Scale_independence"]]<-ggplot(sft$fitIndices,aes(x=Power,y=-sign(slope)*SFT.R.sq))+
    geom_text(aes(label=Power))+
    theme_bw()+
    xlab("Soft Treshold (power)")+
    ylab(paste0("Scale Free Topology Model Fit,",network_type,"R^2"))+
    ggtitle("Scale independence")+
    geom_hline(yintercept = 0.9, col="red")
  p[["Mean_connectivity"]]<-ggplot(sft$fitIndices,aes(x=Power,y=mean.k.))+
    geom_text(aes(label=Power))+
    theme_bw()+
    xlab("Soft Treshold (power)")+
    ylab("Mean Connectivity")+
    ggtitle("Mean connectivity")
  # Set power
  if(is.na(sft$powerEstimate)|any(sft$fitIndices[sft$fitIndices$Power==sft$powerEstimate,"slope"]>0)){
    power<-sft$fitIndices[sft$fitIndices$SFT.R.sq==max(sft$fitIndices$SFT.R.sq[sft$fitIndices$slope<0]),"Power"]
  } else {
    power<-sft$powerEstimate
  }
  # Calculate adjacency matrix
  message("Calculating adjacency matrix...")
  A <- WGCNA::adjacency(datExpr =data, power = power, type = network_type)
  # Calculate topological overlap matrix TOM
  message("Calculating topological overlap matrix...")
  dissTOM <- WGCNA::TOMdist(A, TOMType = network_type)
  # Return
  res<-list(TOM=dissTOM, sft_power=power, plots=p)
  return(res)
}

define_modules<-function(data,TOM, cut_height=0.3){
  # Hierarchical clustering
  geneTree <-flashClust::flashClust(as.dist(TOM), method = "average")
  # Dynamic tree cut
  moduleLabelsManual <- dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = TOM, method = "hybrid",
                                                      deepSplit = 2, pamRespectsDendro = F, minClusterSize = 30)#, cutHeight = cut_height
  moduleColorsManual <- WGCNA::labels2colors(moduleLabelsManual)
  # Merge highly correlated clusters
  data<-data.frame(t(data),stringsAsFactors = F)
  merge<-WGCNA::mergeCloseModules(data,moduleColorsManual, cutHeight = cut_height)
  moduleColorsManual<-merge$colors
  # Get genes associated with each module
  gene_cluster_annot<-data.frame(cbind(names(data),moduleColorsManual),stringsAsFactors = F)
  names(gene_cluster_annot)<-c("Gene","Cluster_color")
  # Calculate kMEs
  MElist <- WGCNA::moduleEigengenes(data, colors = moduleColorsManual,nPC = 6)
  varexp<-MElist$varExplained
  MEs <- MElist$eigengenes
  MEs = WGCNA::orderMEs(MEs)
  datKME <- WGCNA::signedKME(data,MEs)
  # Return
  res<-list(cluster_genes=gene_cluster_annot, module_colors = moduleColorsManual, tree=geneTree, KMEs=datKME)
  return(res)
}

getGO<-function(geneset,universe,species=c("Homo sapiens","Mus musculus")){
  all_genes<-append(rep(0,length(geneset)),rep(1,length(universe[!universe%in%geneset])))
  names(all_genes)<-append(geneset,universe[!universe%in%geneset])
  GOBPdata <- new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(p) p == 0, description = "Test", annot = topGO::annFUN.org, mapping=if(species=="Homo sapiens"){"org.Hs.eg.db"} else {if(species=="Mus musculus"){"org.Mm.eg.db"}}, ID="SYMBOL")
  GOCCdata <- new("topGOdata", ontology = "CC", allGenes = all_genes, geneSel = function(p) p == 0, description = "Test", annot = topGO::annFUN.org, mapping=if(species=="Homo sapiens"){"org.Hs.eg.db"} else {if(species=="Mus musculus"){"org.Mm.eg.db"}}, ID="SYMBOL")
  GOMFdata <- new("topGOdata", ontology = "MF", allGenes = all_genes, geneSel = function(p) p == 0, description = "Test", annot = topGO::annFUN.org, mapping=if(species=="Homo sapiens"){"org.Hs.eg.db"} else {if(species=="Mus musculus"){"org.Mm.eg.db"}}, ID="SYMBOL")

  resultFisherBP <- topGO::runTest(GOBPdata, algorithm = "classic", statistic = "fisher")
  resultFisherCC <- topGO::runTest(GOCCdata, algorithm = "classic", statistic = "fisher")
  resultFisherMF <- topGO::runTest(GOMFdata, algorithm = "classic", statistic = "fisher")

  BP<-data.frame(topGO::GenTable(GOBPdata, classicFisher = resultFisherBP, topNodes = 15, numChar=100),stringsAsFactors = F)
  CC<-data.frame(topGO::GenTable(GOCCdata, classicFisher = resultFisherCC, topNodes = 15, numChar=100),stringsAsFactors = F)
  MF<-data.frame(topGO::GenTable(GOMFdata, classicFisher = resultFisherMF, topNodes = 15, numChar=100),stringsAsFactors = F)
  res<-list(BP,CC,MF)
  names(res)<-c("BP","CC","MF")
  return(res)
}

getKEGG<-function(geneset,geneset_IDtype, species){
  if(species=="Homo sapiens"){organism="hsa"}
  if(species=="Mus musculus"){organism="mmu"}
  first<-function(x){x[[1]]}

  geneset<-as.character(AnnotationDbi::mapIds(if(species=="Homo sapiens"){org.Hs.eg.db}else{if(species=="Mus musculus"){org.Mm.eg.db}},
                                              keys = geneset,
                                              column = "ENTREZID",
                                              keytype = geneset_IDtype,
                                              multiVals = first))
  kk<-clusterProfiler::enrichKEGG(gene = geneset,
                                  organism = organism,
                                  pvalueCutoff = 0.05)
  res<-head(kk,100)
  return(res)
}

calc_regulons<-function(gene_lists, motif_db_path, cores=detectCores()-1){
  # Import motif rankings
  motifRankings<-importRankings(motif_db_path)
  # 1. Calculate AUC

  gene_lists<-lapply(gene_lists, FUN=function(x){
    x<-x[x%in%names(motifRankings@rankings)]
  })


  motifs_AUC <- RcisTarget::calcAUC(gene_lists, motifRankings)

  # 2. Select significant motifs, add TF annotation & format as table
  motifEnrichmentTable <- RcisTarget::addMotifAnnotation(motifs_AUC,
                                                         motifAnnot=motifAnnotations_hgnc)

  # 3. Identify significant genes for each motif
  # (i.e. genes from the gene set in the top of the ranking)
  # Note: Method 'iCisTarget' instead of 'aprox' is more accurate, but slower
  motifEnrichmentTable_wGenes <- RcisTarget::addSignificantGenes(motifEnrichmentTable,
                                                                 geneSets=gene_lists,
                                                                 rankings=motifRankings,
                                                                 nCores=cores,
                                                                 method="aprox")
  # Return
  return(motifEnrichmentTable_wGenes)
}

createGeneSetfromGeneClusters<-function(gene_cluster_annotation){
  # Create geneset list
  gset.idx.list<-list()
  modules<-unique(gene_cluster_annotation$Cluster_color)
  for(module in modules){
    gset.idx.list[[module]]<-as.character(gene_cluster_annotation[gene_cluster_annotation$Cluster_color==module,"Gene"])
  }
  return(gset.idx.list)
}

scoreGenes<-function(exprDat, geneList, method = "gsva"){
  # Perform scoring
  if(method=="gsva"){
    gsva_res<-GSVA::gsva(expr = exprDat, gset.idx.list = geneList, method = "gsva")
  }
  # Format scoring
  gsva_t<-data.frame(gsva_res,stringsAsFactors = F)
  gsva_t$Module<-rownames(gsva_t)
  gsva_t<-tidyr::gather(data.frame(gsva_t,stringsAsFactors = F),key="Cell",value = "gsva",1:(ncol(gsva_t)-1))
  gsva_scores<-gsva_t
  #gsva_t$Cell<-gsub("^\\d+|\\d+$", "", gsva_t$Cell)
  gsva_t<-gsva_t%>%
    group_by(Cell,Module)%>%
    summarise(Avggsva = mean(gsva))
  gsva_t<-spread(gsva_t, key = "Module",value = "Avggsva")
  gsva_t<-as.matrix(gsva_t)
  rownames(gsva_t)<-gsva_t[,"Cell"]
  gsva_t<-gsva_t[,-1]
  gsva_mat<-matrix(as.numeric(gsva_t), ncol = ncol(gsva_t),nrow=nrow(gsva_t),dimnames=list(rownames(gsva_t), colnames(gsva_t)))
  # Return
  res<-list(Uncompacted_GSVA = gsva_scores, GSVA_matrix = gsva_mat)
  return(res)
}

calculateVariableImportance<-function(gsva_mat){
  y<-rownames(gsva_mat)
  module_var_imps<-list()
  for(celltype in y){
    tmp<-y
    tmp[tmp!=celltype]<-"Other"
    module_var_imps[[celltype]]<-caret::filterVarImp(x=data.frame(gsva_mat,stringsAsFactors = F),y=as.factor(tmp))
  }
  return(module_var_imps)
}

generateSummary<-function(CenTFinderObject, tfs, kME_cutoff=0.7){
  # Extract TF kMEs
  tf_kMEs<-CenTFinderObject@kMEs[rownames(CenTFinderObject@kMEs)%in%tfs,]
  # Reshape
  tf_kMEs$TF<-rownames(tf_kMEs)
  tf_kMEs<-tidyr::gather(tf_kMEs, key = "Cluster",value = "kME",1:(ncol(tf_kMEs)-1))
  tf_kMEs$Cluster<-gsub(tf_kMEs$Cluster, pattern = "kME", replacement = "")
  # Keep only TFs for which |kME|>kME_cutoff
  tf_kMEs<-tf_kMEs[tf_kMEs$kME>abs(kME_cutoff),]
  # Filter motif enrichment table
  motifs<-CenTFinderObject@motif_enrichment[CenTFinderObject@motif_enrichment$TF_highConf!="",]
  motifs$TF<-sapply(strsplit(motifs$TF_highConf,split=" |;"),FUN=function(x){x[1]})
  names(motifs)[1]<-"Cluster"
  # Merge WGCNA and RCisTarget
  WR_merged<-merge(tf_kMEs,motifs,by=c("TF","Cluster"),all=T)
  # Check for additional ChipSeq
  return(WR_merged)
}

make_PCA<-function(data){
  if(!is.matrix(data)){
    return(NULL)
  }
  pcadata<-prcomp(t(data), scale. = T)
  dat.sum<-summary(pcadata)
  dat<-as.data.frame(pcadata$x)
  res<-ggplot(data = dat, aes(x = PC1, y = PC2, label = rownames(dat))) +
    geom_point()+
    geom_hline(yintercept = 0, colour = "gray65") +
    geom_vline(xintercept = 0, colour = "gray65")+
    geom_text(colour = "black", alpha = 0.8, size = 2, nudge_x = 1)+
    xlab(paste("PC1",round(dat.sum$importance["Proportion of Variance",1]*100,2),"%",sep=", "))+
    ylab(paste("PC2",round(dat.sum$importance["Proportion of Variance",2]*100,2),"%",sep=", "))+
    theme(text = element_text(size=20),
          axis.text.x = element_text(angle=0, hjust=1))+
    ggtitle("PCA plot")
  return(res)
}

