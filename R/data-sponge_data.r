#' Sponge \emph{A. aerophoba} study
#'
#' This study investigated the relationship between metabolite concentration
#' and microbial abundance of specific sponge tissues. The samples were
#' collected from two types of tissues (Ectosome vs. Choanosome) and processed
#' on two separate denaturing gradient gels in electrophoresis. This study
#' includes relative abundance data only and a completely balanced
#' batch x treatment design.
#'
#'
#' @name sponge_data
#' @docType data
#' @usage data('sponge_data')
#' @format A list containing four data sets:
#' \describe{
#' \item{X.tss}{The relative abundance data represented as a data
#' frame with 32 samples and 24 OTUs.}
#' \item{X.clr}{The centered log ratio transformed data including 32 samples
#' and 24 OTUs from the \code{X.tss}.}
#' \item{Y.trt}{A factor of sponge tissues for each sample that is the effect
#' of interest in the sponge study.}
#' \item{Y.bat}{A factor of electrophoresis gels where the sample processed for
#' each sample treated as the batch effect.}}
#'
#'
#' @return None.
#' @references
#' \insertRef{sacristan2011exploring}{PLSDAbatch}
#' @source The raw data were downloaded from the referenced article. Filtering
#' and normalisation described in our package vignette.
#' @keywords datasets
NULL
