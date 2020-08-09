########################################################
###   CenTFinder                                     ###
###   Main Function definitions: Integration         ###
###   Author: Jonathan De Smedt                      ###
########################################################

#' @include classes.R
#' @include aux_microarray_retrieval.R
#' @include aux_microarray_analysis.R
#' @include aux_WGCNA_core_analyses.R
#' @include aux_misc.R
#' @include main_functions_preprocessing.R
#' @include main_functions_WGCNA.R
#' @include main_functions_module_level_analyses.R

#' Integrate analyses.
#'
#' \code{integrateAnalyses} integrates gene co-expression analysis, transcription factor binding motif enrichment analysis, and differential gene expression analysis. Additionally, it creates regulons per cluster per cell type, which can be exported and visualised in Cytoscape.
#'
#' @name integrateAnalyses
#' @title integrateAnalyses
#' @param  CenTFinderObject An instance of the CenTFinder class.
#' @param  tfs A vector of transcription factor gene symbols.
#' @param  kME_cutoff A numeric cut-off for cluster centrality. Cluster centrality is measured with the module eigengene. A kME of 1 indicates complete centrality while a kME of 0 indicates no cluster membership.
#' @return An updated CenTFinderObject with added integrative summary, module GSVA scores, and module variable importances.
#' @examples
#' CenTFinderObject<-integrateAnalyses(CenTFinderObject, tfs=humantfs, target_cell = "LSEC")
#' @export
#' @importFrom tidyr gather
#' @importFrom tidyr spread
#' @importFrom GSVA gsva
#' @importFrom caret filterVarImp
integrateAnalyses<-function(CenTFinderObject, tfs=NULL, kME_cutoff=0.7){
  if(is.null(tfs)){
    stop("Please provide a vector transcription factor gene symbols of the species you are working with.")
  }
  # Generate summary table
  integrative_summary<-generateSummary(CenTFinderObject = CenTFinderObject, tfs=tfs, kME_cutoff=kME_cutoff)
  if(length(CenTFinderObject@additional_DE)==0){
    integrative_summary$TF_diff_expressed<-FALSE
  } else {
    integrative_summary$TF_diff_expressed<-integrative_summary$TF%in%CenTFinderObject@additional_DE$Gene
  }
  integrative_summary$Class<-4-(as.numeric(!is.na(integrative_summary$kME))+as.numeric(!is.na(integrative_summary$NES))+as.numeric(integrative_summary$TF_diff_expressed))
  # Score each sample for each module with GSVA
  geneset<-createGeneSetfromGeneClusters(gene_cluster_annotation = CenTFinderObject@gene_cluster_annotation)
  moduleScores_res<-scoreGenes(exprDat = CenTFinderObject@data[rownames(CenTFinderObject@data)%in%CenTFinderObject@gene_cluster_annotation$Gene,], geneList = geneset, method = "gsva")
  moduleScores<-data.frame(moduleScores_res$GSVA_matrix, stringsAsFactors = F)
  moduleScores<-merge(CenTFinderObject@sample_annotation[,c("Sample","Group")],moduleScores,by.x="Sample",by.y=0,all=F)
  moduleScores$Sample<-NULL
  moduleScores<-moduleScores%>%
    group_by(Group)%>%
    summarise_all(.funs=mean)
  moduleScores<-data.frame(moduleScores, stringsAsFactors = F)
  rownames(moduleScores)<-moduleScores$Group
  moduleScores$Group<-NULL
  # Calculate variable importance of each module
  module_var_imps<-calculateVariableImportance(gsva_mat = moduleScores)
  module_var_imps<-data.frame(bind_cols(module_var_imps),stringsAsFactors = F, row.names = rownames(module_var_imps[[1]]))
  module_var_imps<-module_var_imps[,!grepl(names(module_var_imps),pattern = "Other")]
  # Return results
  CenTFinderObject@integration<-list(Modules_and_regulons = integrative_summary,
                                     Module_scores = moduleScores,
                                     Variable_importances = module_var_imps,
                                     Uncompacted_GSVA = moduleScores_res$Uncompacted_GSVA)
  return(CenTFinderObject)
}



#' Extract network for Cytoscape visualisation
#'
#' \code{extractNetwork} extracts Cytoscape networks from a CenTFinder object for specified modules.
#'
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @param modules Modules of interest.
#' @param class Integers specifying the transcription factor classes to be displayed. Class 1 transcription factors are differentially expressed, are WGCNA cluster-central, and have RcisTarget enriched binding motifs in the promoters of the genes of their respective clusters. Can only be 1, 2, and/or 3.
#' @param DEset Dataset from which to use the fold changes.
#' @param TFs_diff_expr Logical indicating whether or not differentially expressed transcription factors of the module are selected. Transcription factors that are not differentially expressed might play a posttranscriptional role. Defaults to TRUE.
#' @param Nr_diff_targets Integer specifying the maximum number of top differentially expressed targets are included in the Cytoscape network. Defaults to 20.
#' @param target_cell Cell type in which the TFs and/or targets should be differentially expressed. If target_cell and additional_DE are NULL, TFs_diff_expr and Nr_diff_targets will be ignored.
#' @param z_score_cutoff In case target_cell is NULL and additional_DE is FALSE and Z-scores were calculated for the expression data, then this parameter specifies the cutoff to apply.
#' @param p_value_cutoff In case target_cell is NULL and additional_DE is FALSE and t-tests were calculated for the expression data, then this parameter specifies the cutoff to apply.
#' @return A dataframe with the desired Cytoscape network
#' @export
extractNetwork<-function(CenTFinderObject, modules=NULL, class=c(1,2,3), DEset=NULL, TFs_diff_expr=T, Nr_diff_targets=20, target_cell=NULL, z_score_cutoff=2, p_value_cutoff=0.05){
  if(is.null(modules)){
    stop("Please provide modules of interest.")
  }
  if(!is.logical(TFs_diff_expr)){
    stop("TFs_diff_expr should be a logical indicating whether or not to select differentially expressed transcription factors.")
  }
  if(length(CenTFinderObject@integration)==0){
    stop("Please run integrateAnalyses in order to extract networks.")
  }
  if(any(class==1)&!TFs_diff_expr){
    stop("If TFs_diff_expr=F then no TFs can be of class 1.")
  }
  if(is.null(DEset)&is.null(target_cell)){
    if(length(CenTFinderObject@DE)==0){
      stop("Please provide a target_cell.")
    }
  }
  # Extract modules of interest from the integrative summary
  regulons<-CenTFinderObject@integration$Modules_and_regulons[CenTFinderObject@integration$Modules_and_regulons$Cluster%in%modules&CenTFinderObject@integration$Modules_and_regulons$Class%in%class,]
  regulons<-regulons[!is.na(regulons$enrichedGenes),]
  # Reformat regulons
  regulons<-regulons%>%
    group_by(TF,Cluster)%>%
    summarise(Targets = paste(enrichedGenes, collapse=";"))
  regulons$ID<-paste(regulons$TF,regulons$Cluster,sep="-")
  regulons_reformatted<-list()
  for(id in unique(regulons$ID)){
    tmp<-regulons[regulons$ID==id,]
    regulons_reformatted[[id]]<-data.frame(TF = tmp$TF,
                                           Cluster = tmp$Cluster,
                                           Target = strsplit(paste(tmp$Targets,sep=";"),split=";")[[1]],
                                           stringsAsFactors = F)
  }
  regulons<-dplyr::bind_rows(regulons_reformatted)
  regulons<-regulons[!duplicated(regulons),]
  # Extract differentially expressed genes
  if(!is.null(DEset)){
    diffgenes<-CenTFinderObject@additional_DE[[DEset]][,c("Gene","FC")]
  } else if(is.null(target_cell)){
    if(length(CenTFinderObject@additional_DE)==0){
      diffgenes<-data.frame(Gene=rownames(CenTFinderObject@data),
                            FC=0,
                            stringsAsFactors = F)
    } else {
      diffgenes<-CenTFinderObject@additional_DE[[1]][,c("Gene","FC")]
    }
  } else if(any(names(CenTFinderObject@DE[[target_cell]])=="Avg_Z")){
    diffgenes<-CenTFinderObject@DE[[target_cell]][abs(CenTFinderObject@DE[[target_cell]]$Avg_Z)>z_score_cutoff,c("Gene","Avg_Z")]
  } else {
    diffgenes<-CenTFinderObject@DE[[target_cell]][CenTFinderObject@DE[[target_cell]]$AdjPvalue<p_value_cutoff,c("Gene","AdjPvalue")]
  }
  # Filter for differentially expressed targets
  regulons<-regulons[regulons$Target%in%diffgenes$Gene,]
  # Filter for differentially expressed TFs
  if(TFs_diff_expr){
    regulons<-regulons[regulons$TF%in%diffgenes$Gene,]
  } else {
    if(!is.null(DEset)){
      regulons<-regulons[!regulons$TF%in%diffgenes$Gene,]
    }
  }
  # Extend the regulons with fold changes or Z-scores
  regulons<-regulons%>%
    group_by(TF)%>%
    mutate(nTargets = length(unique(Target)))
  if(TFs_diff_expr){
    regulons<-merge(regulons,diffgenes, by.x="TF",by.y="Gene",all=F)
    names(regulons)[5]<-"logFC_TF_or_Z_score"
    regulons<-merge(regulons,diffgenes, by.x="Target",by.y="Gene",all=F)
    names(regulons)[6]<-"logFC_Target_Z_score"
  } else {
    regulons<-merge(regulons,diffgenes, by.x="Target",by.y="Gene",all=F)
    names(regulons)[5]<-"logFC_Target_Z_score"
  }
  # Assign TFs to the cluster with the highest number of their respective targets
  TF_assign<-data.frame(as.matrix(table(regulons[,c("TF","Cluster")])),stringsAsFactors = F)
  TF_assign<-TF_assign%>%
    group_by(TF)%>%
    summarise(Cluster = Cluster[which.max(Freq)])
  TF_assign$Target<-TF_assign$TF


  regulons<-regulons%>%
    group_by(TF, Cluster)%>%
    filter(abs(logFC_Target_Z_score)>abs(logFC_Target_Z_score)[order(abs(logFC_Target_Z_score),decreasing = T)][min(length(logFC_Target_Z_score),20)])


  regulons<-suppressWarnings(bind_rows(regulons,TF_assign))
  regulons<-regulons%>%
    group_by(TF)%>%
    filter(length(Target)>1)%>%
    mutate(nTargets = unique(nTargets[!is.na(nTargets)]),
           logFC_TF_or_Z_score = unique(logFC_TF_or_Z_score[!is.na(logFC_TF_or_Z_score)]))
  #regulons[is.na(regulons)]<-0
  regulons<-regulons%>%
    group_by(TF,Cluster)%>%
    mutate(nTargets = max(nTargets))
  regulons$logFC_Target_Z_score[is.na(regulons$logFC_Target_Z_score)]<-regulons$logFC_TF_or_Z_score[is.na(regulons$logFC_Target_Z_score)]
  return(regulons)
}

#' Compare two CenTFinder analyses
#'
#' \code{compareCTF} allows one to compare two separate CenTFinder analyses.
#'
#' @name compareCTF
#' @title compareCTF
#' @param  ctf1 An instance of the CenTFinder class.
#' @param  ctf2 An instance of the CenTFinder class.
#' @param  module_annot1 List with module colour-to-name annotations.
#' @param  module_annot2 List with module colour-to-name annotations.
#' @return An alluvial plot displaying the module memberships. Following CenTFinder releases may extend this comparison.
#' @export
#' @import dplyr
#' @import alluvial
compareCTF<-function(ctf1, ctf2, module_annot1, module_annot2){
  # Extract the gene cluster annotations
  annot1<-geneAnnotations(ctf1)
  annot2<-geneAnnotations(ctf2)

  # Change module annotations
  annot1$Module<-unlist(module_annot1[match(annot1$Cluster_color, names(module_annot1))])
  annot2$Module<-unlist(module_annot2[match(annot2$Cluster_color, names(module_annot2))])
  annot1$Cluster_color<-NULL
  annot2$Cluster_color<-NULL

  # Merge annotations
  annot<-merge(annot1,annot2,by="Gene", all = T,suffixes = c("_1","_2"))
  annot$Module_1[is.na(annot$Module_1)]<-"Not in panel"
  annot$Module_2[is.na(annot$Module_2)]<-"Not in panel"

  # Aggregate
  annot_sum<-annot%>%
    group_by(Module_1, Module_2)%>%
    summarise(Freq = length(unique(Gene)))

  # Add missing combos
  missing<-data.frame(expand.grid(unique(annot_sum$Module_1), unique(annot_sum$Module_2)), stringsAsFactors = F)
  names(missing)<-c("Module_1","Module_2")
  missing$Freq<-0
  annot_sum<-data.frame(bind_rows(annot_sum, missing), stringsAsFactors = F)
  annot_sum<-annot_sum%>%
    group_by(Module_1, Module_2)%>%
    summarise(Freq = max(Freq))

  return(annot_sum)
}


