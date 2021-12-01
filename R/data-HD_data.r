#' Mice models with Huntington’s disease
#'
#' This study explored differences in microbial composition between Huntington's
#' disease (HD) and wild-type (WT) mice. The samples were collected from the
#' mice with different genotypes and housed in different cages. This study
#' includes a nested batch x treatment design.
#'
#'
#' @name HD_data
#' @docType data
#' @usage data('HD_data')
#' @format A list containing two sub-lists of data \code{FullData},
#' \code{EgData}:
#' \describe{
#' \item{FullData}{A list containing two data sets: \code{X.count}
#' which are the raw data represented as a matrix with 30 samples and 368 OTUs;
#' \code{metadata} which is a data frame containing the meta data information
#' of samples in \code{X.count}.}
#' \item{EgData}{A list containing three data sets: \code{X.clr} which are the
#' filtered and centered log ratio transformed data including 30 samples and 269
#' OTUs from the raw data in the \code{FullData} list; \code{Y.trt} which is a
#' factor of genotypes for each sample that is the effect of interest in the HD
#' study; \code{Y.bat} which is a factor of cages where the mice were housed
#' for each sample treated as the batch effect.}}
#'
#' @return None.
#' @references
#' \insertRef{kong2020microbiome}{PLSDAbatch}
#' @source The raw data were provided by Geraldine Kong and published at
#' the referenced article. Filtering and normalisation described in our package
#' vignette.
#' @keywords datasets
NULL
