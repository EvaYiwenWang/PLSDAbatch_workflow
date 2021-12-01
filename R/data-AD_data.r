#' Anaerobic digestion study
#'
#' This study explored the microbial indicators that could improve the efficacy
#' of anaerobic digestion (AD) bioprocess and prevent its failure. The samples
#' were treated with two different ranges of phenol concentration (effect of
#' interest) and processed at five different dates (batch effect). This study
#' includes a clear and strong batch effect with an approx. balanced
#' batch x treatment design.
#'
#'
#' @name AD_data
#' @docType data
#' @usage data('AD_data')
#' @format A list containing three sub-lists of data \code{FullData},
#' \code{EgData} and \code{CorrectData}:
#' \describe{
#' \item{FullData}{A list containing three data sets: \code{X.count}
#' which are the raw data represented as a matrix with 75 samples and 567 OTUs;
#' \code{metadata} which is a data frame containing the meta data information
#' of samples in \code{X.count}; \code{taxa} which is a data frame containing
#' the taxonomy of each OTU in \code{X.count}.}
#' \item{EgData}{A list containing three data sets: \code{X.clr} which are the
#' filtered and centered log ratio transformed data including 75 samples and
#' 231 OTUs from the raw data in the \code{FullData} list; \code{Y.trt} which
#' is a factor of phenol concentrations for each sample that is the effect of
#' interest in the AD study; \code{Y.bat} which is a factor of sample
#' processing dates for each sample treated as the batch effect.}
#' \item{CorrectData}{A list containing seven data sets before or after batch
#' effect correction using different methods. Each data set includes 75 samples
#' and 231 OTUs.}}
#'
#' @return None.
#' @references
#' \insertRef{chapleur2016increasing}{PLSDAbatch}
#' @source The raw data were provided by Dr. Olivier Chapleur and published at
#' the referenced article. Filtering and normalisation described in our package
#' vignette.
#' @keywords datasets
NULL
