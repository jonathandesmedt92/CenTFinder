########################################################
###   CenTFinder                                     ###
###   Main Function definitions: WGCNA               ###
###   Author: Jonathan De Smedt                      ###
########################################################

#' @include classes.R
#' @include aux_microarray_retrieval.R
#' @include aux_microarray_analysis.R
#' @include aux_WGCNA_core_analyses.R
#' @include aux_misc.R
#' @include main_functions_preprocessing.R
#' @include main_functions_module_level_analyses.R
#' @include main_functions_integration.R


#' Apply Weighted Gene Correlation Network Analysis (WGCNA).
#'
#' \code{applyWGCNA} runs the entire WGCNA pipeline: 1) calculation of the similarity matrix, 2) calculation of the adjacency matrix, 3) calculation of the toplogical overlap matrix (TOM), and 4) the assignment of genes to clusters. Please check https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/ for more information.
#'
#' @name applyWGCNA
#' @title applyWGCNA
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @param network_type This argument specifies whether WGCNA should make use of signed or unsigned networks.
#' @param powers An integer vector of powers to be evaluated in order to obtain a scale-free topology.
#' @param cut_height A numeric specifying the cutting height of the gene tree in order to cluster genes.
#' @return A CenTFinderObject with added TOM, power, connectivity plots, gene cluster annotation, gene tree, and kMEs.
#' @examples
#' CenTFinderObject<-applyWGCNA(CenTFinderObject, network_type="signed")
#' @export
#' @importFrom WGCNA pickSoftThreshold
#' @importFrom WGCNA adjacency
#' @importFrom WGCNA TOMdist
#' @importFrom flashClust flashClust
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom WGCNA labels2colors
#' @importFrom WGCNA mergeCloseModules
#' @importFrom WGCNA moduleEigengenes
#' @importFrom WGCNA signedKME
#' @importFrom WGCNA orderMEs
applyWGCNA<-function(CenTFinderObject, network_type=c("signed","unsigned"),powers=1:30, cut_height=0.3){
  # Check if data is present
  if(all(dim(CenTFinderObject@data)==0)){
    stop("No data is provided.")
  }
  # Calculate TOM
  res<-calc_TOM(data = CenTFinderObject@data[rownames(CenTFinderObject@data)%in%CenTFinderObject@gene_cluster_annotation$Gene,], network_type = network_type,powers=powers)
  CenTFinderObject@TOM<-as.matrix(res$TOM)
  CenTFinderObject@power<-res$sft_power
  CenTFinderObject@plots<-append(CenTFinderObject@plots,res$plots)
  # Return
  return(CenTFinderObject)
}

#' Cut WGCNA tree
#'
#' \code{cutDendrogram} allows custom cutting of the WGCNA dendrogram, as to set the number of clusters.
#'
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @param cut_height Numeric treshold to cut the dendrogram. Default to 0.2.
#' @return A CenTFinder object with a dendrogram and defined modules.
#' @export
cutDendrogram<-function(CenTFinderObject, cut_height=0.2){
  # Define modules
  res<-define_modules(data = CenTFinderObject@data[rownames(CenTFinderObject@data)%in%CenTFinderObject@gene_cluster_annotation$Gene,], TOM = CenTFinderObject@TOM, cut_height=cut_height)
  CenTFinderObject@gene_cluster_annotation<-res$cluster_genes
  CenTFinderObject@module_colors<-res$module_colors
  CenTFinderObject@tree<-res$tree
  CenTFinderObject@kMEs<-res$KMEs
  res$tree
  # Return
  return(CenTFinderObject)
}


#' Change module colours
#'
#' \code{changeModuleColours} allows to change the cluster colours after WGCNA analysis.
#'
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @param colours Vector containing the new module colours. The first colour should be grey, the subsequent colours should be ranked by decreasing module size.
#' @return A CenTFinder object with updated colour scheme.
#' @importFrom WGCNA standardColors
#' @export
changeModuleColours<-function(CenTFinderObject, colours=NULL){
  if(is.null(colours)){
    message("Colours will default to the order of the standard colours.")
    colours<-c("grey",standardColors(100))
  }
  if(!is.null(colours)){
    if(colours[1]!="grey"){
      stop("The first colour should always be grey. Grey is used for genes that could not be classified in any cluster.")
    }
  }
  colours_used<-as.character(names(sort(table(CenTFinderObject@module_colors),decreasing = T)))
  colours_used<-append("grey",colours_used[colours_used!="grey"])
  CenTFinderObject@gene_cluster_annotation$Cluster_color2<-"grey"
  module_colors2<-rep("grey",length(CenTFinderObject@module_colors))
  kme_names<-rep("grey",ncol(CenTFinderObject@kMEs))

  for(i in 2:length(colours_used)){
    CenTFinderObject@gene_cluster_annotation$Cluster_color2[CenTFinderObject@gene_cluster_annotation$Cluster_color==colours_used[i]]<-colours[i]
    module_colors2[CenTFinderObject@module_colors==colours_used[i]]<-colours[i]
    kme_names[gsub(names(CenTFinderObject@kMEs), pattern = "kME", replacement = "")==colours_used[i]]<-paste0("kME",colours[i])
  }
  CenTFinderObject@gene_cluster_annotation$Cluster_color<-CenTFinderObject@gene_cluster_annotation$Cluster_color2
  CenTFinderObject@gene_cluster_annotation$Cluster_color2<-NULL
  CenTFinderObject@module_colors<-module_colors2
  names(CenTFinderObject@kMEs)<-kme_names
  return(CenTFinderObject)
}
