########################################################
###   CenTFinder                                     ###
###   Main Function definitions: Miscellaneous       ###
###   Author: Jonathan De Smedt                      ###
########################################################

#' @include classes.R
#' @include aux_microarray_retrieval.R
#' @include aux_microarray_analysis.R
#' @include aux_WGCNA_core_analyses.R
#' @include aux_misc.R
#' @include main_functions_WGCNA.R
#' @include main_functions_module_level_analyses.R
#' @include main_functions_integration.R

#' Retrieve the differential expression analysis for a given cell type
#'
#' \code{DE} allows one to extract the differential expression analysis for a given cell type (contrasting with the other cell types)
#'
#' @name DE
#' @title DE
#' @param CenTFinderObject An instance of the CenTFinder class
#' @param celltype Cell type for which the differential expression analysis should be extracted.
#' @return A dataframe with the differential expression results
#' @export
DE<-function(CenTFinderObject, celltype){
  if(length(CenTFinderObject@DE)==0){
    stop("No microarray or RNAseq analysis has been performed.")
  }
  if(!any(names(CenTFinderObject@DE)%in%celltype)){
    stop("Provided cell type is not present in the dataset.")
  }
  return(CenTFinderObject@DE[[celltype]])
}

#' Extract the optionally provided differential expression data
#'
#' \code{DEset} allows one to extract the additional and optional differential expression data.
#'
#' @name DEset
#' @title DEset
#' @param CenTFinderObject An instance of the CenTFinder class
#' @return A dataframe with the differential expression data.
#' @export
DEset<-function(CenTFinderObject){
  if(length(CenTFinderObject@additional_DE)==0){
    stop("No additional expression data has been loaded.")
  }
  return(CenTFinderObject@additional_DE)
}

#' Load an additional differential expression dataset
#'
#' \code{'DEset<-'} allows one to load a differential expression analysis.
#'
#' @name DEset<-
#' @title DEset<-
#' @param x An instance of the CenTFinder class
#' @param value A list containing one or more differential expression analyses, with 'Gene', 'FC', and 'Pvalue' columns in each of the dataframes.
#' @return The CenTFinder object updated with a differential expression set. P-values are automatically Benjamini-Hochberg adjusted.
#' @export
'DEset<-'<-function(x, value){
  if(!all(names(value)%in%c("Gene","FC","Pvalue"))){
    stop("The first three column names of each provided DEset should be Gene, FC, and Pvalue.")
  }
  value<-value[,c("Gene","FC","Pvalue")]
  value$AdjPvalue<-p.adjust(value$Pvalue, method = "BH")
  x@additional_DE[[length(x@additional_DE)+1]]<-value
  validObject(x)
  return(x)
}

#' Extract the gene expression matrix
#'
#' \code{expressionData} allows one to extract the gene expression matrix used for analysis (i.e. any filtering will be taken into account).
#'
#' @name expressionData
#' @title expressionData
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @return The gene expression matrix used for CenTFinder analysis.
#' @export
expressionData<-function(CenTFinderObject){
  if(all(dim(CenTFinderObject@data)==0)){
    stop("No expression data has been loaded.")
  }
  if(all(dim(CenTFinderObject@gene_cluster_annotation)==0)){
    res<-CenTFinderObject@data
  } else {
    res<-CenTFinderObject@data[rownames(CenTFinderObject@data)%in%CenTFinderObject@gene_cluster_annotation$Gene,]
  }
  return(res)
}

#' Extract the gene-module annotations
#'
#' \code{geneAnnotations} allows one to extract the gene-module annotations from a CenTFinder object.
#'
#' @name geneAnnotations
#' @title geneAnnotations
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @return A dataframe with the gene cluster annotation.
#' @export
geneAnnotations<-function(CenTFinderObject){
  if(all(dim(CenTFinderObject@gene_cluster_annotation)==c(0,0))){
    stop("No WGCNA analysis has been performed yet. Please run applyWGCNA().")
  }
  return(CenTFinderObject@gene_cluster_annotation)
}

#' Extract the variances of all genes
#'
#' \code{getAllVariances} extracts all gene expression variances.
#'
#' @name getAllVariances
#' @title getAllVariances
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @return A dataframe with the variance of each gene.
#' @export
getAllVariances<-function(CenTFinderObject){
  if(is.null(dim(CenTFinderObject@data))){
    stop("No expression data has been loaded.")
  }
  dat<-data.frame(cbind(rownames(CenTFinderObject@data),matrixStats::rowVars(CenTFinderObject@data)), stringsAsFactors = F)
  names(dat)<-c("Gene","Variance")
  return(dat)
}


#' Extract groups
#'
#' \code{groups} extracts the groups within a CenTFinder object.
#'
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @return A vector containing all unique group names.
#' @export
groups<-function(CenTFinderObject){
  if(all(dim(CenTFinderObject@sample_annotation)==0)){
    stop("No sample annotation has been loaded yet. Please run setSampleAnnotations().")
  }
  return(unique(CenTFinderObject@sample_annotation$Group))
}

#' Extract the GSVA module scores of the CenTFinder object
#'
#' \code{getModuleScores} lets one extract the GSVA module activities from a CenTFinder object.
#'
#' @param CenTFinderObject An instance of the CenTFinder class
#' @return A dataframe with containing the GSVA activity
#' @export
getModuleScores<-function(CenTFinderObject){
  if(length(CenTFinderObject@integration)==0){
    stop("No integrative analysis has been performed yet.")
  }
  return(CenTFinderObject@integration$Module_scores)
}

#' Extract the overlap between WGCNA, RcisTarget, and any optional differential expression sets
#'
#' \code{getOverlap} extracts the overlap or integration of all analyses.
#'
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @return A dataframe containing the overlap of all analyses.
#' @export
getOverlap<-function(CenTFinderObject){
  if(length(CenTFinderObject@integration)==0){
    message("Please run integrateAnalyses in order to get an overlap.")
  } else {
    return(CenTFinderObject@integration$Modules_and_regulons)
  }
}

#' Extract the sample names
#'
#' \code{getSampleNames} allows to extract the sample names of samples used in the CenTFinder analysis.
#'
#' @param CenTFinderObject An instance of the CenTFinder class
#' @return A character vector with the sample names used.
#' @export
getSampleNames<-function(CenTFinderObject){
  if(is.null(dim(CenTFinderObject@data))){
    stop("No expression data has been loaded.")
  }
  return(colnames(CenTFinderObject@data))
}


#' Extract the variable importances for each WWGCNA module
#'
#' \code{getVariableImportances} extracts the variable importance for each module. This is a wrapper function for the 'filterVarImp' function of the 'caret' package.
#'
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @return A dataframe with the variable importance of each module.
#' @export
getVariableImportances<-function(CenTFinderObject){
  if(length(CenTFinderObject@integration)==0){
    stop("No integrative analysis has been performed yet.")
  }
  return(CenTFinderObject@integration$Variable_importances)
}

#' Extract the module colours
#'
#' \code{modules} allows one to extract the module colours of a CenTFinder object.
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @return A character vector containing the module colours in the CenTFinder analysis.
#' @export
modules<-function(CenTFinderObject){
  if(all(dim(CenTFinderObject@gene_cluster_annotation)==0)){
    stop("No WGCNA analysis has been performed yet. Please run applyWGCNA().")
  }
  return(unique(as.character(CenTFinderObject@gene_cluster_annotation$Cluster_color)))
}

#' Extract microarray platforms
#'
#' \code{platforms} returns the microarray platform types that were used in the CenTFinder analysis.
#'
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @return A character vector containing the platform types used.
#' @export
platforms<-function(CenTFinderObject){
  if(all(dim(CenTFinderObject@platforms)==0)){
    stop("No platform information has been loaded yet. Please run getPlatformInfo().")
  }
  return(CenTFinderObject@platforms)
}


#' Remove samples
#'
#' \code{removeSamples} allows one to remove samples from a CenTFinder object.
#'
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @param samples Names or indices of the samples to be removed.
#' @return An updated CenTFinder object in which samples were removed.
#' @export
removeSamples<-function(CenTFinderObject, samples=NULL){
  if(is.null(samples)){
    stop("No samples were provided.")
  }
  if(class(samples)=="character"){
    if(sum(colnames(CenTFinderObject@data)%in%samples)!=length(samples)){
      stop("Not all provided sample names are currently present in this CenTFinder object.")
    }
    CenTFinderObject@data<-CenTFinderObject@data[,!colnames(CenTFinderObject@data)%in%samples]
  }
  if(class(samples)=="integer"){
    if(max(samples)>ncol(CenTFinderObject@data)){
      stop("Some provided indices are out of bounds.")
    }
    CenTFinderObject@data<-CenTFinderObject@data[,-samples]
  }
  CenTFinderObject@sample_annotation<-CenTFinderObject@sample_annotation[CenTFinderObject@sample_annotation$Sample%in%colnames(CenTFinderObject@data),]
  # Create a subset of the SDRF
  sdrf<-read.table("SDRF.txt", sep = "\t", header = T, stringsAsFactors = F)
  sdrf_used<-sdrf[sdrf$filenames%in%CenTFinderObject@sample_annotation$File,]
  write.table(sdrf_used, "SDRF_used.txt",sep="\t",row.names = F)
  return(CenTFinderObject)
}

#' Extract sample annotation
#'
#' \code{sampleAnnotations} extracts the sample annotation of a CenTFinder object.
#'
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @return A dataframe containing the sample annotation of the samples in a CenTFinder object.
#' @export
sampleAnnotations<-function(CenTFinderObject){
  if(all(dim(CenTFinderObject@sample_annotation)==0)){
    stop("No sample annotation has been loaded yet. Please run setSampleAnnotations().")
  }
  return(CenTFinderObject@sample_annotation)
}

#' Extract thresholding power
#'
#' \code{tresholdingPower} extracts the soft-thresholding power used in the CenTFinder analysis.
#'
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @return The applied soft-thresholding power.
#' @export
thresholdingPower<-function(CenTFinderObject){
  if(length(CenTFinderObject@power)==0){
    stop("No WGCNA analysis has been performed yet. Please run applyWGCNA().")
  }
  return(CenTFinderObject@power)
}
















