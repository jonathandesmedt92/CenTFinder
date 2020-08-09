############################################################
###   CenTFinder                                         ###
###   Main Function definitions: Module-level analyses   ###
###   Author: Jonathan De Smedt                          ###
############################################################

#' @include classes.R
#' @include aux_microarray_retrieval.R
#' @include aux_microarray_analysis.R
#' @include aux_WGCNA_core_analyses.R
#' @include aux_misc.R
#' @include main_functions_preprocessing.R
#' @include main_functions_WGCNA.R
#' @include main_functions_integration.R


#' Analyse gene co-expression clusters.
#'
#' \code{analyseClusters} performs several downstream analyses on WGCNA gene co-expression clusters, including Gene Ontology analysis (using topGO), KEGG pathway enrichment analysis (using ClusterProfiler), principal component analysis (PCA, using prcomp of stats package), and transcription factor binding motif enrichment analysis (using RCisTarget).
#'
#' @name analyseClusters
#' @title analyseClusters
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @param GO A logical indicating whether Gene Ontology analysis should be performed.
#' @param KEGG A logical indicating whether KEGG pathway enrichment should be performed.
#' @param PCA A logical indicating whether principal component analysis should be performed.
#' @param RCisTarget A logical indicating whether transcription factor binding motif enrichment analysis should be performed. This function is a wrapper for usage of the RCisTarget package.
#' @param species Species of the samples. Only Latin names are supported (e.g. Homo sapiens, Mus musculus, ...)
#' @param cores Number of cores to be used for RCisTarget analysis. Defaults to total number of cores minus one.
#' @param motif_db_path Path to the motif db needed for RCisTarget analysis. For more info, please check http://bioconductor.org/packages/release/bioc/vignettes/RcisTarget/inst/doc/RcisTarget.html.
#' @return An updated CenTFinderObject with cluster analyses.
#' @examples
#' motif_db_path<-"hg19-tss-centered-10kb-7species.mc9nr.feather"
#' CenTFinderObject<-analyseClusters(CenTFinderObject, GO=T, KEGG=T, PCA=T, RCisTarget=T,species="Homo sapiens", cores=10, motif_db_path=motif_db_path)
#' @export
#' @importClassesFrom topGO topGOdata
#' @import topGO
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom AnnotationDbi mapIds
#' @importFrom parallel detectCores
#' @import RcisTarget
analyseClusters<-function(CenTFinderObject, GO=T, KEGG=T, PCA=T, RCisTarget=T,species=c("Homo sapiens","Mus musculus"), cores=detectCores()-1, motif_db_path=NULL){
  if(is.null(motif_db_path)){
    stop("Please provide the path where the motif database is stored.")
  }
  GOs<-list()
  KEGGs<-list()
  PCAs<-list()
  RCisTargets<-NULL
  # Perform GO analyses
  if(GO){
    message("Performing Gene Ontology analyses on the identified WGCNA modules...")
    for(module in unique(CenTFinderObject@gene_cluster_annotation$Cluster_color)){
      GOs[[module]]<-getGO(geneset = as.character(CenTFinderObject@gene_cluster_annotation[CenTFinderObject@gene_cluster_annotation$Cluster_color==module,"Gene"]),
                           universe = as.character(unique(CenTFinderObject@gene_cluster_annotation$Gene)),
                           species = species)
    }
    CenTFinderObject@GO<-GOs
  }
  # Perform KEGG analyses
  if(KEGG){
    message("Performing KEGG analyses on the identified WGCNA modules...")
    for(module in unique(CenTFinderObject@gene_cluster_annotation$Cluster_color)){
      KEGGs[[module]]<-getKEGG(geneset = as.character(CenTFinderObject@gene_cluster_annotation[CenTFinderObject@gene_cluster_annotation$Cluster_color==module,"Gene"]),
                               geneset_IDtype = "SYMBOL",
                               species = species)
    }
    CenTFinderObject@KEGG<-KEGGs
  }
  # Perform PCA analyses
  if(PCA){
    message("Creating PCA plots per cluster...")
    for(module in unique(CenTFinderObject@gene_cluster_annotation$Cluster_color)){
      genes<-as.character(CenTFinderObject@gene_cluster_annotation[CenTFinderObject@gene_cluster_annotation$Cluster_color==module,"Gene"])
      PCAs[[module]]<-make_PCA(data = CenTFinderObject@data[rownames(CenTFinderObject@data)%in%genes,])
    }
    CenTFinderObject@plots$PCA_plots_per_cluster<-PCAs
  }
  # Perform RCisTarget analysis
  if(RCisTarget){
    message("Running RCisTarget analyses on the identified WGCNA modules...")
    gene_lists<-list()
    for(color in unique(CenTFinderObject@gene_cluster_annotation$Cluster_color)){
      gene_lists[[color]]<-as.character(CenTFinderObject@gene_cluster_annotation[CenTFinderObject@gene_cluster_annotation$Cluster_color==color,"Gene"])
    }
    RCisTargets<-calc_regulons(gene_lists = gene_lists,
                               motif_db_path = motif_db_path,
                               cores = cores)
    RCisTargets$Analysis<-"modules"
    RCisTargets$Level<-"primary"
    CenTFinderObject@motif_enrichment<-RCisTargets
  }
  # Return
  return(CenTFinderObject)
}


#' Add markers.
#'
#' \code{addMarkers} allows to add RCisTarget analysis on a set of provided markers for the cell type of interest.
#'
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @param markers A vector of gene symbols for a set of markers.
#' @param cores Number of cores to be used for RCisTarget analysis. Defaults to total number of cores minus one.
#' @param motif_db_path Path to the motif db needed for RCisTarget analysis. For more info, please check http://bioconductor.org/packages/release/bioc/vignettes/RcisTarget/inst/doc/RcisTarget.html.
#' @return An updated CenTFinderObject with RCisTarget analysis of provided markers.
#' @examples
#' motif_db_path<-"hg19-tss-centered-10kb-7species.mc9nr.feather"
#' CenTFinderObject<-addMarkers(CenTFinderObject,markers=c("FCGR2B","PECAM1","STAB1","STAB2","LYVE1","CLEC4G","CLEC4M"),cores=10,motif_db_path=motif_db_path)
#' @export
addMarkers<-function(CenTFinderObject, markers=NULL, cores=detectCores()-1, motif_db_path=NULL){
  if(is.null(markers)){
    stop("Please provide a vector marker gene symbols.")
  }
  if(is.null(motif_db_path)){
    stop("Please provide the path where the motif database is stored.")
  }
  RCisTargets<-calc_regulons(gene_lists = list(markers=markers),
                             motif_db_path = motif_db_path,
                             cores = cores)
  RCisTargets$Analysis<-"markers"
  RCisTargets$Level<-"primary"
  CenTFinderObject@motif_enrichment<-bind_rows(CenTFinderObject@motif_enrichment,RCisTargets)
  return(CenTFinderObject)
}


