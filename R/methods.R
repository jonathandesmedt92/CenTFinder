########################################################
###   CenTFinder                                     ###
###   Class definitions                              ###
###   version 21032018                               ###
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
#' @include main_functions_integration.R


# Make S3 generics S4-aware
setOldClass("plot")
setOldClass("summary")
setOldClass("show")

# Plot needs to be used for..
# plotModuleScores
# scaleIndependencePlot
# MeanConnectivityPlot
# modulePCA
# tree
# plotVariance
# plotExpression
# clusterdendrogram
# plotgo
# plotkegg
# plotmodulescores

# groupDistribution


#' Plot an aspect of a CenTFinder object
#'
#' Plot allows one to visualise the different aspects of the CenTFinder analysis.
#'
#' @param object An instance of the CenTFinder class.
#' @param aspect An aspect of the CenTFinder class. Aspect can be either of the following: "array distribution" (plots the microarray platform distribution), "group distribution" (plots the distribution of the different groups (cell types)), "expression" (plots the expression of one or more genes in the meta-analysis data), "markers" (plots the expression of top identified markers in the meta-analysis data), "tfs" (plots the expression of top expressed tfs in the cell type of interest in the meta-analysis data), "variance" (plots the least and most variable genes in the meta-analysis data), "scale independence" (plots the correlation of each soft-thresholding connectivity distribution with the scale-free topology connectivity distribution), "mean connectivity" (plots the mean connectivity in the network in function of the soft-thresholding power used), "tree" (plots the WGCNA dendrogram tree), "cluster dendrogram" (plots the WGCNA cluster dendrogram), "GSVA heatmap" (plots a heatmap of GSVA activity scores for each module and each cell type/group of cells), "GSVA boxplot" (plots a boxplot of the GSVA scores of the cell type/group of interest versus the other cell types/groups), "module PCA" (plots a PCA plot for a given module), "GO" (plots the gene ontology enerichment for a given module), "KEGG" (plots the KEGG enrichment for a given module),"RcisTarget" (plots the RcisTarget results), and "integrative" (plots the overlap and final ranking of all analyses).
#' @import ggplot2
#' @import gplots
#' @export
setMethod(f = "plot",
          signature = "CenTFinder",
          definition = function(x, aspect, ...){
            switch(aspect,
                   "array distribution" = arrayDistribution(x),
                   "group distribution" = groupDistribution(x),
                   "expression" = plotExpression(x,...),
                   "markers" = plotMarkers(x,...),
                   "tfs" = plotTFs(x,...),
                   "variance" = plotVariance(x,...),
                   "scale independence" = scaleIndependencePlot(x),
                   "mean connectivity" = meanConnectivityPlot(x),
                   "tree" = tree(x),
                   "cluster dendrogram" = clusterDendrogram(x),
                   "GSVA heatmap" = plotModuleScores(x,...),
                   "GSVA boxplot" = plotGSVA(x,...),
                   "module PCA" = modulePCA(x,...),
                   "GO" = plotGO(x,...),
                   "KEGG" = plotKEGG(x,...),
                   "RcisTarget" = plotRcisTarget(x,...),
                   "integrative" = plotIntegration(x,...))
          })


#' Get a summary of a CenTFinder object
#'
#' Summary reveals the main features of the CenTFinder object, including analyses performed.
#'
#' @param object An instance of the CenTFinder class.
#' @export
setMethod(f = "summary",
          signature="CenTFinder",
          definition = function(object,...){
            CenTFinderObject<-object
            if(all(dim(CenTFinderObject@motif_enrichment)!=0)){
              tfs<-sapply(strsplit(CenTFinderObject@motif_enrichment$TF_highConf,split = " \\("),FUN=function(x){x[1]})
              tfs<-tfs[!is.na(tfs)]
              tfs<-unlist(strsplit(tfs,split = "; "))
              tfs<-unique(tfs)
            }
            cat("\t\n",
                "CenTFinder object\n",
                sprintf("Expression data: %s\n", ifelse(all(dim(CenTFinderObject@data)==0),"Empty",paste0(nrow(CenTFinderObject@data)," genes x ",ncol(CenTFinderObject@data)," samples."))),
                sprintf("Filtered expression data: %s\n", ifelse(all(dim(CenTFinderObject@data)==0),"Empty",paste0(sum(rownames(CenTFinderObject@data)%in%CenTFinderObject@gene_cluster_annotation$Gene)," genes x ",ncol(CenTFinderObject@data)," samples."))),
                sprintf("Distinct platforms: %s\n", ifelse(all(dim(CenTFinderObject@platforms)==0),"None", as.character(length(unique(CenTFinderObject@platforms$Platform[CenTFinderObject@platforms$Used]))))),
                sprintf("Distinct cell types: %s\n", ifelse(all(dim(CenTFinderObject@sample_annotation)==0),"None",as.character(length(unique(CenTFinderObject@sample_annotation$Group))))),
                sprintf("Additional differential expression data: %s\n", ifelse(length(CenTFinderObject@additional_DE)==0,"None",ifelse(length(CenTFinderObject@additional_DE)==1,paste0(1," dataset loaded."),paste0(length(CenTFinderObject@additional_DE)," datasets loaded.")))),
                "\t\n",
                "-----------------------------\n",
                "WGCNA\n",
                "-----------------------------\n",
                sprintf("TOM: %s\n", ifelse(all(dim(CenTFinderObject@TOM)==0),"None","Calculated")),
                sprintf("Soft-thresholding power: %s\n", ifelse(length(CenTFinderObject@power)==0,"None",as.character(CenTFinderObject@power))),
                sprintf("Number of modules detected: %i\n",ifelse(length(CenTFinderObject@module_colors)==0,0,length(unique(CenTFinderObject@module_colors)))),
                "\t\n",
                "-----------------------------\n",
                "RCisTarget\n",
                "-----------------------------\n",
                sprintf("Motifs detected: %s\n", ifelse(all(dim(CenTFinderObject@motif_enrichment)==0),"None",as.character(nrow(CenTFinderObject@motif_enrichment)))),
                sprintf("Number of TFs with enriched binding motifs: %s\n", ifelse(all(dim(CenTFinderObject@motif_enrichment)==0),"None",as.character(length(tfs))))

            )
          })

#' Show a summary of a CenTFinder object
#'
#' Show reveals the main features of the CenTFinder object, including analyses performed.
#'
#' @param object An instance of the CenTFinder class.
#' @export
setMethod(f = "show",
          signature="CenTFinder",
          definition = function(object){
            CenTFinderObject<-object
            if(all(dim(CenTFinderObject@motif_enrichment)!=0)){
              tfs<-sapply(strsplit(CenTFinderObject@motif_enrichment$TF_highConf,split = " \\("),FUN=function(x){x[1]})
              tfs<-tfs[!is.na(tfs)]
              tfs<-unlist(strsplit(tfs,split = "; "))
              tfs<-unique(tfs)
            }
            cat("\t\n",
                "CenTFinder object\n",
                sprintf("Expression data: %s\n", ifelse(all(dim(CenTFinderObject@data)==0),"Empty",paste0(nrow(CenTFinderObject@data)," genes x ",ncol(CenTFinderObject@data)," samples."))),
                sprintf("Filtered expression data: %s\n", ifelse(all(dim(CenTFinderObject@data)==0),"Empty",paste0(sum(rownames(CenTFinderObject@data)%in%CenTFinderObject@gene_cluster_annotation$Gene)," genes x ",ncol(CenTFinderObject@data)," samples."))),
                sprintf("Distinct platforms: %s\n", ifelse(all(dim(CenTFinderObject@platforms)==0),"None", as.character(length(unique(CenTFinderObject@platforms$Platform[CenTFinderObject@platforms$Used]))))),
                sprintf("Distinct cell types: %s\n", ifelse(all(dim(CenTFinderObject@sample_annotation)==0),"None",as.character(length(unique(CenTFinderObject@sample_annotation$Group))))),
                sprintf("Additional differential expression data: %s\n", ifelse(length(CenTFinderObject@additional_DE)==0,"None",ifelse(length(CenTFinderObject@additional_DE)==1,paste0(1," dataset loaded."),paste0(length(CenTFinderObject@additional_DE)," datasets loaded.")))),
                "\t\n",
                "-----------------------------\n",
                "WGCNA\n",
                "-----------------------------\n",
                sprintf("TOM: %s\n", ifelse(all(dim(CenTFinderObject@TOM)==0),"None","Calculated")),
                sprintf("Soft-thresholding power: %s\n", ifelse(length(CenTFinderObject@power)==0,"None",as.character(CenTFinderObject@power))),
                sprintf("Number of modules detected: %i\n",ifelse(length(CenTFinderObject@module_colors)==0,0,length(unique(CenTFinderObject@module_colors)))),
                "\t\n",
                "-----------------------------\n",
                "RCisTarget\n",
                "-----------------------------\n",
                sprintf("Motifs detected: %s\n", ifelse(all(dim(CenTFinderObject@motif_enrichment)==0),"None",as.character(nrow(CenTFinderObject@motif_enrichment)))),
                sprintf("Number of TFs with enriched binding motifs: %s\n", ifelse(all(dim(CenTFinderObject@motif_enrichment)==0),"None",as.character(length(tfs))))

            )
          })



