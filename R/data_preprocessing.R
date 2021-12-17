####################################################
#' Prefiltering for Microbiome Data
#'
#' This function prefilters the data to remove samples or microbial variables
#' with excess zeroes.
#'
#' @param data The data to be prefiltered. The samples in rows and variables
#' in columns.
#' @param keep.spl The minimum counts of a sample to be kept.
#' @param keep.var The minimum proportion of counts of a variable to be kept.
#'
#' @return \code{PreFL} returns a list that contains the following components:
#' \item{data.filter}{The filtered data matrix.}
#' \item{sample.idx}{The indexs of samples kept.}
#' \item{var.idx}{The indexs of variables kept.}
#' \item{zero.prob}{The proportion of zeros of the input data.}
#'
#' @author Yiwen Wang, Kim-Anh Lê Cao
#'
#'
#' @references
#' \insertRef{le2016mixmc}{PLSDAbatch}
#'
#' @export
#'
#' @examples
#' data('AD_data')
#' ad.count <- AD_data$FullData$X.count # microbial count data
#' ad.filter.res <- PreFL(data = ad.count)
#'
#' # The proportion of zeroes of the AD count data
#' ad.zero.prob <- ad.filter.res$zero.prob
#'
#' # The filtered AD count data
#' ad.filter <- ad.filter.res$data.filter
#'
#'
PreFL <- function(data,
                  keep.spl = 10,
                  keep.var = 0.01){

  # zero prob
  zero.prob <- mean(data == 0)

  # sample
  keep.sample = which(rowSums(data) > keep.spl)
  data <- data[keep.sample, ]

  # variable
  keep.var = which(colSums(data)*100/(sum(colSums(data))) > keep.var)
  data.out <- as.matrix(data[, keep.var])

  result = list(data.filter = data.out,
                sample.idx = keep.sample,
                var.idx = keep.var,
                zero.prob = zero.prob)

  return(invisible(result))
}


