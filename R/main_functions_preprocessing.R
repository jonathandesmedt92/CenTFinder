########################################################
###   CenTFinder                                     ###
###   Main Function definitions: preprocessing       ###
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



#' Get ArrayExpress accession codes.
#'
#' \code{getAE_accession_codes} retrieves the ArrayExpress accession codes of studies containing microarrays matching the queried keywords.
#'
#' @name getAE_accession_codes
#' @title getAE_accession_codes
#' @param query One or more keywords for querying the ArrayExpress database.
#' @return A vector containing the retrieved ArrayExpress accession codes.
#' @examples
#' getAE_accession_codes(c("endothelial cell AND Homo sapiens","LSEC AND Homo sapiens"))
#' @export
#' @importFrom RCurl getURL
getAE_accession_codes<-function(query){
  base_url = "https://www.ebi.ac.uk/arrayexpress/xml/v3/experiments?keywords="
  responses<-list()
  for(i in 1:length(query)){
    query_tmp<-paste(strsplit(query[i],split=" ")[[1]], collapse = "+")
    query_tmp<-paste(base_url, query_tmp, sep = "")
    responses[[i]]<-retrieveURLdata(URL=query_tmp,FUN = RCurl::getURL)
    responses[[i]]<-strsplit(responses[[i]], split = "<experiment><id>")
    responses[[i]]<-lapply(responses[[i]],FUN = strsplit, split = "</id><accession>")
    responses[[i]]<-sapply(responses[[i]][[1]], FUN = function(x){strsplit(x[2],split = "</accession>")[[1]][1]})
  }
  response<-unlist(responses)
  response<-response[!is.na(response)]
  response<-response[!duplicated(response)]
  return(response)
}



#' Get sample-data relationship data.
#'
#' \code{getBatchSDRF} retrieves the ArrayExpress Sample Data Relationships Files (SDRF) for the provided accession codes.
#'
#' After using \code{\link{getAE_accession_codes}} this function should be used to retrieve the SDRF data of the acquired accession codes. This function further filters the datasets for the desired assay type (i.e. sequencing or array).
#'
#' @param AEcodes A vector containing ArrayExpress accession codes.
#' @param technology Has to be either "array" or "sequencing".
#' @param label Is NULL in case 'technology = sequencing', otherwise it defaults to 'biotin'. Here, different array types can be specified (i.e. biotin arrays, dual arrays, ...).
#' @param cores Integer specifying the number of cores to be used. Defaults to the number of available cores minus one.
#' @return A dataframe containing the merged SDRF files.
#' @examples
#' getBatchSDRF(AEcodes=c("E-MTAB-4025","E-MTAB-2495"), technology = "array", label = "biotin")
#' @export
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom data.table fread
#' @importFrom parallel mclapply
getBatchSDRF<-function(AEcodes, technology=c("array","sequencing"),label=if(technology=="array"){"biotin"}else{NULL}, cores=detectCores()-1){
  if(is.null(technology)){
    stop("Please specify whether you want arrays or sequencing assays.")
  }
  if(!class(AEcodes)=="character"){
    stop("Please provide ArrayExpress accession codes as a character vector.")
  }
  # Get batch sdrf files
  AEcodes_list<-as.list(AEcodes)
  names(AEcodes_list)<-AEcodes
  if(.Platform$OS.type == "unix"){
    sdrf<-suppressWarnings(parallel::mclapply(AEcodes_list, FUN = getSDRF, mc.cores = cores))
  } else {
    sdrf<-suppressWarnings(lapply(AEcodes_list, FUN = getSDRF))
  }
  # Bind the files
  columns<-unique(unlist(sapply(sdrf, colnames)))
  sdrf<-lapply(sdrf, FUN = function(x){
    missing_columns<-columns[!columns%in%colnames(x)]
    if(length(missing_columns)!=0&length(missing_columns)!=length(columns)){
      missing_column_data<-data.frame(matrix(ncol=length(missing_columns), nrow = nrow(x)), stringsAsFactors = F)
      names(missing_column_data)<-missing_columns
      x<-data.frame(cbind(x,missing_column_data),stringsAsFactors = F)
    }
    if(length(missing_columns)!=length(columns)){
      x<-x[,columns]
    }
    if(length(missing_columns)==length(columns)){
      x<-NA
    }
    return(x)
  })
  sdrf<-sdrf[sapply(sdrf,FUN=function(x){class(x)})!="logical"]
  sdrf<-dplyr::bind_rows(sdrf)
  # Select technology
  if(any(names(sdrf)=="technologytype")){
    sdrf<-switch (technology,
                  "array" = sdrf[sdrf$technologytype=="array assay"|is.na(sdrf$technologytype),],
                  "sequencing" = sdrf[sdrf$technologytype=="sequencing assay"|is.na(sdrf$technologytype),]
    )
  }
  # Select labelling
  if(!is.null(label)){
    if(any(names(sdrf)=="label")){
      sdrf<-switch (label,
                    "biotin" = sdrf[sdrf$label=="biotin",]
      )
    }
  }
  sdrf<-sdrf[!is.na(sdrf$AEcode),]
  return(sdrf)
}

#' Download expression data.
#'
#' \code{downloadAE} downloads the expression data associated with provided ArrayExpress accession codes.
#'
#' @param  AEcodes A vector of ArrayExpress accession codes.
#' @param path Path of the directory where the downloaded expression data should be stored. In case this directory does not yet exist, it will be created automatically.
#' @param processed_data Logical indicating whether or not to retrieve pre-processed expression data in case the raw data is missing. Defaults to FALSE.
#' @return None.
#' @examples
#' downloadAE(AEcodes = c("E-MTAB-4025","E-MTAB-2495"), path="data.downloads")
#' @export
downloadAE<-function(AEcodes, path="data.downloads", processed_data=F){
  if(!class(AEcodes)=="character"){
    stop("Please provide ArrayExpress accession codes as a character vector.")
  }
  # Base URL for data retrieval
  base<-"ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/"
  # Check and/make directory
  if(!dir.exists(path)){
    dir.create(path)
  }
  # Loop over the AEcodes
  for(i in 1:length(AEcodes)){
    url<-paste(base,substr(AEcodes[i],3,6),"/",AEcodes[i],"/",sep = "")
    url<-tryCatch(RCurl::getURL(url), error = function(e){cat('In error handler\n'); print(e); e})
    if(any(class(url)=="error")){
      next
    }
    url<-unlist(strsplit(url,split=" ")[[1]])

    if(length(grep(pattern = "raw",url))!=0){
      url<-url[grep(pattern = "raw",url)]
      for(j in 1:length(url)){
        url[j]<-paste(strsplit(url[j], split = ".zip")[[1]][1],".zip",sep = "")
        url[j]<-paste(base,"/",substr(AEcodes[i],3,6),"/",AEcodes[i],"/",url[j],sep = "")
        downloadURLdata(URL=url[j],destfile = paste(path,"/",AEcodes[i],"_",j,".zip",sep=""))
      }
    } else {
      if(processed_data&length(grep(pattern = "processed",url))!=0){
        url<-url[grep(pattern = "processed",url)]
        for(j in 1:length(url)){
          url[j]<-paste(strsplit(url[j], split = ".zip")[[1]][1],".zip",sep = "")
          url[j]<-paste(base,"/",substr(AEcodes[i],3,6),"/",AEcodes[i],"/",url[j],sep = "")
          downloadURLdata(URL=url[j],destfile = paste(path,"/",AEcodes[i],"_",j,".zip",sep=""))
        }
      }
    }
    flush.console()
    message(paste("Data from ",AEcodes[i]," downloaded. - ",round(i/length(AEcodes)*100,2),"% completed.",sep = ""))
  }
  message("All data downloaded.")
  # Unzip data
  files<-list.files(path=path, pattern = ".zip")
  for(i in 1:length(files)){
    unzip(zipfile = paste(path,"/",files[i],sep = ""), exdir=paste0("./",path), overwrite = F)
    message(paste0("Zip file ",i," of ",length(files)," unzipped."))
  }
  message("All data unpacked.")
}

#' Unpack and filter downloaded expression data.
#'
#' \code{filterAEdata} unzips all zip files present in the path provided in \code{\link{downloadAE}}. Only expression data of the desired species are retained. Subsequently, the expression data is filtered for keywords.
#'
#' @param sdrf The SDRF data describing the expression data in the provided path. This should be the output of \code{\link{getBatchSDRF}}.
#' @param path Path to directory where expression data is stored.
#' @param species Desired species of the expression data. Latin names are supported (e.g. 'Homo sapiens').
#' @param keywords A character vector containing keywords for files of interest. Expression data for which at least one of the annotation fields contains at least one keyword are kept.
#' @param exclusion_keywords A character vector containing keywords for expression data that (for user-specificshould be excluded from analysis (e.g. cancer data). Defaults to NULL.
#' @param allowed_extensions A character vector containing the extensions of data types allowed for analysis.
#' @examples
#' sdrf<-getBatchSDRF(AEcodes=c("E-MTAB-4025","E-MTAB-2495"), technology = "array", label = "biotin")
#' filterAEdata(sdrf = sdrf, path="data.downloads",species = "Homo sapiens", keywords = keywords,exclusion_keywords=c("oma","tumor","tumour","cancer","malignant"))
#' @export
filterAEdata<-function(sdrf, path="data.downloads", species="Homo sapiens", keywords=NULL, exclusion_keywords=NULL, allowed_extensions=c(".cel",".txt")){
  # if(sum(names(sdrf)%in%c("characteristicsorganism","characteristicsorganismpart","characteristicscelltype","characteristicscellline","filenames"))<5){
  #   stop("The SDRF file should contain columns named \'characteristicsorganism\',\'characteristicsorganismpart\',\'characteristicscelltype\',\'characteristicscellline\',\'filenames\'")
  # }
  if(!any(names(sdrf)=="characteristicsorganism")){
    stop("The SDRF file should contain a column named \'characteristicsorganism\'.")
  }

  if(is.null(path)){
    stop("Please provide a path.")
  }
  if(is.null(keywords)){
    stop("Please provide keywords.")
  }
  # Filter sdrf for species of interest
  sdrf<-sdrf[sdrf$characteristicsorganism%in%species,]
  # Filter sdrf for keywords
  keywords<-paste0(keywords,collapse = "|")
  keep<-apply(sdrf, MARGIN = 1, FUN=function(x){
    any(grepl(pattern = keywords, x=x))
  })
  sdrf<-sdrf[keep,]
  # Filter sdrf for exclusion keywords
  if(!is.null(exclusion_keywords)){
    exclusion_keywords<-paste0(exclusion_keywords,collapse = "|")
    keep<-apply(sdrf, MARGIN = 1, FUN=function(x){
      !any(grepl(pattern = exclusion_keywords, x=x))
    })
    sdrf<-sdrf[keep,]
  }
  # Add celltype
  sdrf$celltype<-apply(sdrf[,colnames(sdrf)%in%c("characteristicsorganismpart","characteristicscelltype","characteristicscellline")], MARGIN = 1, FUN = function(x){paste(x, collapse = ";")})
  sdrf$celltype<-gsub(pattern = ";NA|NA;|;  ",replacement = "",sdrf$celltype)
  # Filter for files that were actually downloaded
  file_annot<-data.frame(Filename = list.files(path = path),
                         Filename_match = list.files(path = path),
                         stringsAsFactors = F)
  file_annot<-file_annot[grepl(pattern = paste0(allowed_extensions, collapse = "|"),file_annot$Filename_match, ignore.case = T),]
  file_annot$Filename_match<-gsub(pattern = "_sample_table", replacement = "", file_annot$Filename_match)
  file_annot$Filename_match<-gsub(pattern = paste0(allowed_extensions, collapse = "|"), replacement = "", file_annot$Filename_match, ignore.case = T)

  sdrf$File_id<-apply(sdrf[,colnames(sdrf)%in%c("sourcename","filenames","arraydatafile","derivedarraydatamatrixfile")], MARGIN = 1, FUN = function(x){paste(x, collapse = " ")})
  sdrf$Rank<-1:nrow(sdrf)
  file_pos<-sapply(file_annot$Filename_match, function(i) grep(i,sdrf$File_id , fixed=TRUE))
  file_pos<-file_pos[sapply(file_pos,length)!=0]
  for(i in 1:length(file_pos)){
    file_pos[[i]]<-data.frame(Filename_match = names(file_pos)[i],
                              Rank = file_pos[[i]],
                              stringsAsFactors = F)
  }
  file_pos<-dplyr::bind_rows(file_pos)
  file_annot<-merge(file_annot,file_pos,by="Filename_match",all=F)
  sdrf<-merge(sdrf, file_annot, by="Rank", all=F)

  # Delete unused files
  files<-list.files(path = path)
  unused<-files[!files%in%sdrf$Filename]
  unused<-paste0(path,"/",unused,sep="")
  file.remove(unused)
  #sdrf<-sdrf[sdrf$filenames%in%list.celfiles(path=path,ignore.case = T),]
  #sdrf<-sdrf[!duplicated(sdrf$filenames),]
  # Write sdrf file
  sdrf$Rank<-NULL
  sdrf$File_id<-NULL
  sdrf$Filename_match<-NULL

  write.table(sdrf,"SDRF.txt",sep = "\t")
}

#' Gather information about the downloaded platforms.
#'
#' \code{getPlatformInfo} evaluates the platform type and frequency of the downloaded expression data.
#'
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @param path The path pointing to the downloaded CEL files.
#' @param platform_freq_cutoff Minimal percentage of arrays allowed of any given platform. Defaults to zero.
#' @return A CenTFinderObject with updated platform annotation.
#' @examples
#' CenTFinderObject<-getPlatformInfo(CenTFinderObject,path = "data.downloads", platform_freq_cutoff = 0.00258)
#' @export
#' @importFrom affy list.celfiles
#' @importFrom affyio read.celfile.header
getPlatformInfo<-function(CenTFinderObject,filenames=NULL, path, platform_freq_cutoff=0){
  if(is.null(filenames)){
    stop("Please provide file names.")
  }
  if(is.null(path)){
    stop("Please provide the data download path.")
  }
  # Get and set platform information
  result<-get_microarrays_info(path,platform_freq_cutoff, filenames)
  CenTFinderObject@sample_annotation<-result
  # Notify in case platform annotations are missing
  platform_annots<-data.frame(Platform = unique(result$Platform),
                              Annotation = "No annotation found.",
                              stringsAsFactors = F)
  platform_annots[platform_annots$Platform%in%c("HG-U133_Plus_2","NuGO_Hs1a520180","HG-U219","HT_HG-U133A"),"Annotation"]<-"From package"
  CenTFinderObject@platforms<-platform_annots
  if(any(platform_annots$Annotation=="No annotation found.")){
    msg<-paste0("Please provide an annotation file for :", paste(as.character(platform_annots[platform_annots$Annotation=="No annotation found.","Platform"]),collapse = ", "),".\nArrays for which no probe annotation file is provided will automatically be excluded from analysis.")
    message(msg)
  }
  return(CenTFinderObject)
}


#' Add platform probe annotation data.
#'
#' \code{loadProbeAnnotation} allows to add probe annotation for platforms other than HG-U133_Plus_2, NuGO_Hs1a520180 , HG-U219, or HT_HG-U133A. Arrays not belonging to aforementioned four platform types and for which no probe annotation is provided through \code{loadProbeAnnotation} will be excluded from further analysis.
#'
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @param platform Name of the platform.
#' @param probes A vector of probe identifiers.
#' @param symbols A vector of gene symbols corresponding to the respective probes.
#' @return A CenTFinderObject with additional probe annotation data.
#' @export
loadProbeAnnotation<-function(CenTFinderObject, platform=NULL, probes=NULL, symbols=NULL){
  if(all(dim(CenTFinderObject@platforms)==0)){
    stop("No platform information is loaded. Please make sure to use this function only when you provide raw CEL files. Please run 'getPlatformInfo' first.")
  }
  if(is.null(platform)){
    stop("Please provide platform name.")
  }
  if(is.null(probes)){
    stop("Please provide probe IDs.")
  }
  if(is.null(symbols)){
    stop("Please provide gene symbols.")
  }
  CenTFinderObject@probe_annotations[[platform]]<-data.frame(ACCNUM = probes,
                                                             SYMBOL = symbols,
                                                             stringsAsFactors = F)
  CenTFinderObject@platforms[CenTFinderObject@platforms$Platform==platform,"Annotation"]<-"Annotation loaded."
  return(CenTFinderObject)
}


#' Set array annotations.
#'
#' \code{setArrayAnnotations} sets for each sample the provided cell type and accession code. This information is easily put together based on the SDRF information.
#'
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @param files A vector of CEL file names in the download path.
#' @param groups Groups corresponding to the provided CEL file names. Groups could indicate the cell type, but any other grouping is possible.
#' @param AEcodes Accession codes corresponding to the provided CEL file names.
#' @return A CenTFinderObject with updated sample annotations.
#' @export
setArrayAnnotations<-function(CenTFinderObject, files, groups, AEcodes){
  groups<-gsub(groups, pattern = " ",replacement = ".")
  annot<-data.frame(File = files,
                    Group = groups,
                    AEcodes = AEcodes,
                   stringsAsFactors = F)
  if(!any(names(CenTFinderObject@sample_annotation)=="Group")){
    CenTFinderObject@sample_annotation<-merge(CenTFinderObject@sample_annotation,annot,by="File",all.x=T,all.y=F)
    CenTFinderObject@sample_annotation<-CenTFinderObject@sample_annotation[complete.cases(CenTFinderObject@sample_annotation),]
    CenTFinderObject@sample_annotation<-CenTFinderObject@sample_annotation%>%
      group_by(Group)%>%
      mutate(Sample = paste0(Group,"_",1:length(Group)))
  } else {
    message("Sample annotations have already been set.")
  }
  return(CenTFinderObject)
}




#' Analyse microarrays.
#'
#' \code{analyseMicroarrays} performs RMA normalisation per platform type, subsequently runs Combat to remove batch effects, calculates which platforms to keep to maximise both number of samples and number of genes, and finally calculates either Z scores or T test p-values for all genes and for each cell type.
#'
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @param platforms_keep Platform types that should be kept for analyses. Normally, arrays of platforms with too few genes are removed from analysis. However, one can overrule this by providing the platforms in this argument.
#' @param sig Significance level for T-tests. Defaults to 0.05.
#' @param statistic Statistical method for gene ranking for each cell type. 'Wilcoxon' calculates Wilcoxon rank sum tests for all genes for each cell type with all other cell types as contrast. 'T-test' calculates t-tests for all genes for each cell type with all other cell types as contrast . 'Z-scores' merely calculates the Z-score of expression for each gene across all samples.
#' @return A CenTFinderObject updated with a normalised expression data matrix as well as differential expression analysis.
#' @examples
#' CenTFinderObject<-analyseMicroarrays(CenTFinderObject, platforms_keep = "HG-U219", statistic = "Zscores")
#' @export
#' @importFrom oligo read.celfiles
#' @importFrom Biobase exprs
#' @importFrom oligo rma
#' @importFrom limma normalizeQuantiles
#' @importFrom sva ComBat
#' @importFrom matrixStats rowMads
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import nugohs1a520180.db
#' @import hgu133a.db
#' @import hgu133plus2.db
#' @import hgu219.db
#' @import hgu133plus2.db
analyseMicroarrays<-function(CenTFinderObject, platforms_keep=NULL, sig=0.05, statistic = c("Wilcoxon","T-test","Zscores","none"), contrast = NULL){
  if(all(dim(CenTFinderObject@platforms)==0)){
    stop("No arrays were detected. Please run first 'getPlatformInfo' and 'loadProbeAnnotation'.")
  }
  if(!is.null(contrast)){
    if(length(contrast)!=2){
      stop("Contrast assessed by Student's T-test requires two sample groups.")
    }
  }
  # Analyse arrays from multiple platforms
  CenTFinderObject<-mp_analyse_microarray(CenTFinderObject, platform_keep=platforms_keep)
  # Set sample names
  colnames(CenTFinderObject@data)<-CenTFinderObject@sample_annotation$Sample[match(colnames(CenTFinderObject@data),CenTFinderObject@sample_annotation$File)]
  CenTFinderObject@sample_annotation<-CenTFinderObject@sample_annotation[CenTFinderObject@sample_annotation$Sample%in%colnames(CenTFinderObject@data),]
  # Initialise gene_cluster_annotation
  CenTFinderObject@gene_cluster_annotation<-data.frame(Gene = rownames(CenTFinderObject@data),
                                                       Cluster_color = NA,
                                                       stringsAsFactors = F)
  # Create a subset of the SDRF
  if(file.exists("SDRF.txt")){
    sdrf<-read.table("SDRF.txt", sep = "\t", header = T, stringsAsFactors = F)
    sdrf_used<-sdrf[sdrf$filenames%in%CenTFinderObject@sample_annotation$File,]
    write.table(sdrf_used, "SDRF_used.txt",sep="\t",row.names = F)
  } else {
    message("No SDRF.txt found.")
  }
  # Perform DE analysis
  CenTFinderObject@DE<-switch(statistic,
                              "Wilcoxon" = microarrayDE_Wilcoxon_rank_sum_test(CenTFinderObject = CenTFinderObject, sig=sig, contrast = contrast),
                              "T-test" = microarrayDE_t_test(CenTFinderObject = CenTFinderObject, sig=sig, contrast = contrast),
                              "Zscores" = microarrayDE_Z_Scores(CenTFinderObject = CenTFinderObject),
                              "none" = list())
  return(CenTFinderObject)
}


#' Load expression matrix.
#'
#' \code{loadExpression} loads a clean expression matrix into a CenTFinder object. Expression values can be log2(intensity) values from microarray data, log2(TPM+1) counts from bulk RNA sequencing data or log2(CPM+1) values from TempO-Seq data.
#'
#' @param CenTFinderObject An instance of the CeNTFinder class.
#' @param expr The expression matrix, either microarray, bulkRNAseq, or TempO-Seq data.
#' @param colData An annotation dataframe for the samples. Should have two columns: File and Group. File names should match the column names in the expression data.
#' @param ID_type Type of gene identifier. Currently, 'SYMBOL' and 'ENSEMBL' ID types are supported.
#' @param species Species of the samples. Currently, only 'Homo sapiens' and 'Mus musculus' are supported.
#' @return A CenTFinder object updated with the expression matrix.
#' @export
#' @import dplyr
loadExpression<-function(CenTFinderObject, expr, colData, ID_type=c("SYMBOL","ENSEMBL"), species = c("Homo sapiens","Mus musculus"), statistic = c("Wilcoxon","T-test","Zscores","none"), contrast = NULL){
  if(missing(CenTFinderObject)){
    stop("Please provide a CenTFinder object to update.")
  }
  if(missing(expr)){
    stop("Please provide an expression matrix.")
  }
  if(missing(colData)){
    stop("Please provide colData")
  }
  if(length(ID_type)==2){
    ID_type<-ID_type[1]
  }
  if(length(species)==2){
    species<-"Homo sapiens"
  }
  # Format gene IDs
  if(ID_type == "ENSEMBL"){
    geneset<-as.character(AnnotationDbi::mapIds(if(species=="Homo sapiens"){org.Hs.eg.db}else{if(species=="Mus musculus"){org.Mm.eg.db}},
                                                keys = rownames(expr),
                                                column = "SYMBOL",
                                                keytype = ID_type,
                                                multiVals = function(x){x[[1]]}))
    expr$Symbol<-geneset
    expr<-expr[!is.na(expr$Symbol),]
    expr<-expr[!duplicated(expr$Symbol),]
    rownames(expr)<-expr$Symbol
    expr$Symbol<-NULL
  }
  # Remove hyphens from genes
  rownames(expr)<-gsub(rownames(expr), pattern = "-", replacement = "_")
  # Normalise and store
  CenTFinderObject@data<-as.matrix(expr)
  # Format the annotation
  colData<-colData[complete.cases(colData[,c("File","Group")]),]
  colData<-colData%>%
    group_by(Group)%>%
    mutate(Sample = paste0(Group,"_",1:length(Group)))
  CenTFinderObject@sample_annotation<-colData
  # Rename the expression data with the sample IDs instead of the file IDs
  for(i in 1:length(colData$File)){
    colnames(CenTFinderObject@data)[colnames(CenTFinderObject@data)==colData$File[i]]<-colData$Sample[i]
  }
  # Initialise gene_cluster_annotation
  CenTFinderObject@gene_cluster_annotation<-data.frame(Gene = rownames(CenTFinderObject@data),
                                                       Cluster_color = NA,
                                                       stringsAsFactors = F)
  # Perform DE analysis
  CenTFinderObject@DE<-switch(statistic,
                              "Wilcoxon" = microarrayDE_Wilcoxon_rank_sum_test(CenTFinderObject = CenTFinderObject, sig=sig, pre_analysed=T),
                              "T-test" = microarrayDE_t_test(CenTFinderObject = CenTFinderObject, sig=sig, contrast = contrast, pre_analysed=T),
                              "Zscores" = microarrayDE_Z_Scores(CenTFinderObject = CenTFinderObject, pre_analysed=T),
                              "none" = list())
  return(CenTFinderObject)
}


#' Filter genes.
#'
#' \code{filterGenes} filters a CenTFinder expression matrix for the most variable genes.
#'
#' This function applies different filtering criteria for transcription factors and the other genes. Transcription factors are typically retained if they vary more than two-fold across all samples. This behaviour can be regulated with the 'TF_range_cutoff' argument. The most variant genes that are not transcriptional regulators are retained.
#' @param CenTFinderObject An instance of the CenTFinder class.
#' @param tfs A vector of transcription factor gene symbols.
#' @param TF_range_cutoff An numeric specifying the minimum log2 fold change of a transcription factor between any two samples. Defaults to 1.
#' @param max_genes An integer specifying the maximum total number of genes to be retained for analysis. Increasing this parameter will increase computation time. Computational complexity is of the order O(n^2). Defaults to 9000.
#' @return  A CenTFinderObject with a filtered expression matrix.
#' @examples
#' humantfs<-as.character(read.table("allhumantfs.txt",sep = "\t")[,1])
#' CenTFinderObject<-filterGenes(CenTFinderObject,tfs=humantfs, TF_range_cutoff=1, max_genes=9000)
#' @export
filterGenes<-function(CenTFinderObject, tfs=NULL, TF_range_cutoff=1, max_genes=9000){
  if(is.null(tfs)){
    stop("Please provide a vector of transcription factors for the species you are working with.")
  }
  CenTFinderObject@gene_cluster_annotation<-CenTFinderObject@gene_cluster_annotation[CenTFinderObject@gene_cluster_annotation$Gene%in%filter_genes(CenTFinderObject@data, tfs=tfs, TF_range_cutoff = TF_range_cutoff, max_genes = max_genes),]
  return(CenTFinderObject)
}


