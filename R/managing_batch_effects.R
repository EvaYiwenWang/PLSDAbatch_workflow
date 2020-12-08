###################################
#' Linear Regression
#'
#' This function fits linear regression on each variable of the input data including linear model,
#' linear mixed model and adjusts p-values for multiple comparisons using one of several methods.
#'
#' @importFrom lmerTest lmer
#'
#' @param data A data frame that contains the response variables in the model.
#' Samples as rows and variables as columns.
#' @param trt A factor or a class vector for the treatment grouping information (categorical outcome variable).
#' @param batch A factor or a class vector for the batch grouping information (categorical outcome variable).
#' @param model The model to be used for fitting, either 'linear model' or 'linear mixed model'.
#' @param p.adjust.method The method to be used for p-value adjustment, either "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr" or "none".
#'
#' @return \code{linear_regres} returns a list that contains the following components:
#' \item{p_adjusted}{The adjusted p-values for each variable from the input data.}
#' \item{model}{The model used for fitting.}
#' \item{p.adjust.method}{The method used for p-value adjustment.}
#'
#'
#' @export
#'
#' @examples
#' data('AD_data')
#' ad.clr <- AD_data$EgData$X.clr
#' ad.batch = AD_data$EgData$Y.bat
#' ad.trt = AD_data$EgData$Y.trt
#' ad.lm <- linear_regres(data = ad.clr, trt = ad.trt, batch = ad.batch, model = 'linear model')
#' ad.p.adj <- ad.lm$p_adjusted
#'
#'
linear_regres <- function(data,
                          trt,
                          batch = NULL,
                          model = 'linear model',
                          p.adjust.method = 'fdr'){

  if(!is.null(batch)){
    if(model == 'linear model'){
      p <- apply(data, 2, FUN = function(x){
        res.lm <- lm(x ~ trt + batch)
        summary.res <- summary(res.lm)
        p <- summary.res$coefficients[2,4]
      })
      p.adj <- p.adjust(p, method = p.adjust.method)
    }

    if(model == 'linear mixed model'){
      p <- apply(data, 2, FUN = function(x){
        res.lmm <- lmer(x ~ trt + (1|batch))
        summary.res <- summary(res.lmm)
        p <- summary.res$coefficients[2,4]
      })
      p.adj <- p.adjust(p, method = p.adjust.method)
    }}else{
      p <- apply(data, 2, FUN = function(x){
        res.lm <- lm(x ~ trt)
        summary.res <- summary(res.lm)
        p <- summary.res$coefficients[2,4]
      })
      p.adj <- p.adjust(p, method = p.adjust.method)
    }
  result = list(p_adjusted = p.adj,
                model = model,
                p.adjust.method = p.adjust.method)

  return(invisible(result))

}


#################################
#' Percentile score
#'
#' This function converts the relative abundance of microbial variables (i.e. bacterial taxa) in case (i.e. disease) samples
#' to percentiles of the equivalent variables in control (i.e. healthy) samples.
#'
#' @importFrom Rdpack reprompt
#' @param df A data frame that contains the microbial variables and required to be converted into percentile scores.
#' Samples as rows and variables as columns.
#' @param control.index A numeric vector that contains the indexes of control samples.
#'
#' @return A data frame of percentile scores for each microbial variable and each sample.
#' @keywords Internal
#'
#' @references
#' \insertRef{gibbons2018correcting}{PLSDAbatch}
#'
#' @export
#'
#' @examples
#'
#' data('AD_data')
#' ad.clr <- AD_data$EgData$X.clr
#' ad.batch = AD_data$EgData$Y.bat
#' ad.trt = AD_data$EgData$Y.trt
#'
#' batch1.idx <- which(ad.batch == '09/04/2015')
#' trt.batch1 <- ad.trt[batch1.idx]
#' ad.batch1 <- ad.clr[batch1.idx, ]
#' control.idx = which(trt.batch1 == '0-0.5')
#'
#' ad.batch1.pn = percentileofscore(ad.batch1, control.idx)
#'
#'
percentileofscore = function(df, control.index){
  df.percentile = df
  df.percentile[1:nrow(df), 1:ncol(df)] = NA
  for(i in 1:ncol(df)){
    control = sort(df[control.index, i])
    for(j in 1:nrow(df)){
      percentile.strick = sum(control < df[j, i])/length(control)
      percentile.weak = (length(control) - sum(control > df[j, i]))/length(control)
      percentile = (percentile.strick + percentile.weak)/2
      df.percentile[j, i] = percentile

    }
  }
  return(df.percentile)
}

################################
#' Percentile Normalisation
#'
#' This function corrects for batch effects in case-control microbiome studies.
#' Briefly, the relative abundance of microbial variables (i.e. bacterial taxa) in case (i.e. disease) samples are
#' converted to percentiles of the equivalent variables in control (i.e. healthy) samples within a batch prior to
#' pooling data across batches. Pooled batches must have similar case and control cohort definitions.
#'
#' @importFrom Rdpack reprompt
#' @param data A data frame that contains the microbial variables and required to be corrected for batch effects.
#' Samples as rows and variables as columns.
#' @param batch A factor or a class vector for the batch grouping information (categorical outcome variable).
#' @param trt A factor or a class vector for the treatment grouping information (categorical outcome variable).
#' @param ctrl.grp Character, the name of control samples (i.e. healthy).
#'
#' @return A data frame that corrected for batch effects.
#'
#' @references
#' \insertRef{gibbons2018correcting}{PLSDAbatch}
#'
#' @export
#'
#' @examples
#'
#' data('AD_data')
#' ad.clr <- AD_data$EgData$X.clr
#' ad.batch = AD_data$EgData$Y.bat
#' ad.trt = AD_data$EgData$Y.trt
#' ad.PN <- percentile_norm(data = ad.clr, batch = ad.batch,
#'                          trt = ad.trt, ctrl.grp = '0-0.5')
#'
percentile_norm = function(data = data, batch = batch, trt = trt, ctrl.grp){
  batch = as.factor(batch)
  trt = as.factor(trt)
  trt.list = list()
  data.pn.df = data.frame()
  for(i in 1:nlevels(batch)){
    trt.each.b = trt[batch == levels(batch)[i]]
    trt.list[[i]] = trt.each.b
    data.each.b.pn = percentileofscore(data[batch == levels(batch)[i],],
                                       which(trt.each.b == ctrl.grp))
    data.pn.df = rbind(data.pn.df,data.each.b.pn)
  }
  names(trt.list) = levels(batch)
  data.pn.df.reorder = data.pn.df[rownames(data), ]
  return(data.pn.df.reorder)
}


