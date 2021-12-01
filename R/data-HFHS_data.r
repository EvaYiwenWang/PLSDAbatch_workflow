#' High fat high sugar diet study
#'
#' This study aimed to investigate the effect of high fat high sugar (HFHS)
#' diet on the mouse microbiome. The samples were collected at different days
#' from the mice treated with two types of diets (HFHS vs. normal) and housed
#' in different cages and different facilities. This study includes weak batch
#' effects with a extremely unbalanced (nested) batch x treatment design between
#' cages and diets and an approx. balanced design between days and diets or
#' facilities and diets.
#'
#'
#' @name HFHS_data
#' @docType data
#' @usage data('HFHS_data')
#' @format A list containing three sub-lists of data \code{FullData},
#' \code{EgData} and \code{CorrectData}:
#' \describe{
#' \item{FullData}{A list containing three data sets: \code{X.count}
#' which are the raw data represented as a matrix with 250 samples and 4524
#' OTUs; \code{metadata} which is a data frame containing the meta data
#' information of samples in \code{X.count}; \code{taxa} which is a matrix
#' containing the taxonomy of each OTU in \code{X.count}.}
#' \item{EgData}{A list containing four data sets: \code{X.clr} which are the
#' filtered and centered log ratio transformed data including 54 samples and 419
#' OTUs from a subset of the raw data in the \code{FullData} list; \code{Y.trt}
#' which is a factor of diet types for each sample that is the effect of
#' interest in the HFHS study; \code{Y.bat} which is a factor of facilities
#' where the mice were housed for each sample treated as the batch effect;
#' \code{taxa} which is a matrix containing the taxonomy of each OTU in
#' \code{X.clr}.}
#' \item{CorrectData}{A list containing seven data sets before or after batch
#' effect correction using different methods. Each data set includes 149 samples
#' and 515 OTUs. In addition, the list also contains \code{Y.trt} which is a
#' factor of diet types for each sample that is the effect of interest
#' in the HFHS study; \code{Y.bat} which is a factor of collection dates for
#' each sample treated as the batch effect.}}
#'
#' @return None.
#' @references
#' \insertRef{susin2020variable}{PLSDAbatch}
#' @source The raw data were provided by Prof. Kim-Anh Le Cao and published at
#' the referenced article. Filtering and normalisation described in our package
#' vignette.
#' @keywords datasets
NULL
