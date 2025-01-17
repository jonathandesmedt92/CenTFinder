% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_functions_preprocessing.R
\name{analyseMicroarrays}
\alias{analyseMicroarrays}
\title{Analyse microarrays.}
\usage{
analyseMicroarrays(
  CenTFinderObject,
  platforms_keep = NULL,
  sig = 0.05,
  statistic = c("Wilcoxon", "T-test", "Zscores", "none"),
  contrast = NULL
)
}
\arguments{
\item{CenTFinderObject}{An instance of the CenTFinder class.}

\item{platforms_keep}{Platform types that should be kept for analyses. Normally, arrays of platforms with too few genes are removed from analysis. However, one can overrule this by providing the platforms in this argument.}

\item{sig}{Significance level for T-tests. Defaults to 0.05.}

\item{statistic}{Statistical method for gene ranking for each cell type. 'Wilcoxon' calculates Wilcoxon rank sum tests for all genes for each cell type with all other cell types as contrast. 'T-test' calculates t-tests for all genes for each cell type with all other cell types as contrast . 'Z-scores' merely calculates the Z-score of expression for each gene across all samples.}
}
\value{
A CenTFinderObject updated with a normalised expression data matrix as well as differential expression analysis.
}
\description{
\code{analyseMicroarrays} performs RMA normalisation per platform type, subsequently runs Combat to remove batch effects, calculates which platforms to keep to maximise both number of samples and number of genes, and finally calculates either Z scores or T test p-values for all genes and for each cell type.
}
\examples{
CenTFinderObject<-analyseMicroarrays(CenTFinderObject, platforms_keep = "HG-U219", statistic = "Zscores")
}
