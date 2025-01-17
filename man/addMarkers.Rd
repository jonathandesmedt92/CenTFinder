% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_functions_module_level_analyses.R
\name{addMarkers}
\alias{addMarkers}
\title{Add markers.}
\usage{
addMarkers(
  CenTFinderObject,
  markers = NULL,
  cores = detectCores() - 1,
  motif_db_path = NULL
)
}
\arguments{
\item{CenTFinderObject}{An instance of the CenTFinder class.}

\item{markers}{A vector of gene symbols for a set of markers.}

\item{cores}{Number of cores to be used for RCisTarget analysis. Defaults to total number of cores minus one.}

\item{motif_db_path}{Path to the motif db needed for RCisTarget analysis. For more info, please check http://bioconductor.org/packages/release/bioc/vignettes/RcisTarget/inst/doc/RcisTarget.html.}
}
\value{
An updated CenTFinderObject with RCisTarget analysis of provided markers.
}
\description{
\code{addMarkers} allows to add RCisTarget analysis on a set of provided markers for the cell type of interest.
}
\examples{
motif_db_path<-"hg19-tss-centered-10kb-7species.mc9nr.feather"
CenTFinderObject<-addMarkers(CenTFinderObject,markers=c("FCGR2B","PECAM1","STAB1","STAB2","LYVE1","CLEC4G","CLEC4M"),cores=10,motif_db_path=motif_db_path)
}
