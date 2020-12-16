###################################
#' Linear Regression
#'
#' This function fits linear regression (linear model or linear mixed model) on each microbial variable
#' and includes treatment and batch effects as covariates. It generates p-values, adjusted p-values for
#' multiple comparisons, and evaluation metrics of model quality.
#'
#'
#' @importFrom lmerTest lmer
#' @importFrom performance rmse r2
#'
#' @param data A data frame that contains the response variables for the linear regression.
#' Samples as rows and variables as columns.
#' @param trt A factor or a class vector for the treatment grouping information (categorical outcome variable).
#' @param batch A factor or a class vector for the batch grouping information (categorical outcome variable).
#' @param type The type of model to be used for fitting, either 'linear model' or 'linear mixed model'.
#' @param p.adjust.method The method to be used for p-value adjustment, either "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr" or "none".
#'
#' @return \code{linear_regres} returns a list that contains the following components:
#' \item{type}{The type of model used for fitting.}
#' \item{model}{Each object fitted.}
#' \item{raw.p}{The p-values for each response variable.}
#' \item{adj.p}{The adjusted p-values for each response variable.}
#' \item{p.adjust.method}{The method used for p-value adjustment.}
#' \item{R2}{The proportion of variation in the response variable that is explained by the predictor variables.
#' A higher R2 indicates a better model. Results for 'linear model' only.}
#' \item{adj.R2}{Adjusted R2 for many predictor variables in the model. Results for 'linear model' only.}
#' \item{cond.R2}{The proportion of variation in the response variable that is explained by the "complete" model with all covariates.
#' Results for 'linear mixed model' only. Similar to \code{R2} in linear model.}
#' \item{marg.R2}{The proportion of variation in the response variable that is explained by the fixed effects part only.
#' Results for 'linear mixed model' only.}
#' \item{RMSE}{The average error performed by the model in predicting the outcome for an observation.
#' A lower RMSE indicates a better model.}
#' \item{RSE}{also known as the model \eqn{sigma}, is a variant of the RMSE adjusted for the number of predictors in the model.
#' A lower RSE indicates a better model.}
#' \item{AIC}{A penalisation value for including additional predictor variables to a model.
#' A lower AIC indicates a better model.}
#' \item{BIC}{is a variant of AIC with a stronger penalty for including additional variables to the model.}
#'
#' @note \code{R2, adj.R2, cond.R2, marg.R2, RMSE, RSE, AIC, BIC} all include the results of two models:
#' (i) the full input model; (ii) a model without batch effects. It can help to decide whether it is better to include
#' batch effects.
#'
#' @references
#' \insertRef{daniel2020performance}{PLSDAbatch}
#'
#' @export
#'
#' @examples
#' data('AD_data')
#' ad.clr <- AD_data$EgData$X.clr
#' ad.batch = AD_data$EgData$Y.bat
#' ad.trt = AD_data$EgData$Y.trt
#' ad.lm <- linear_regres(data = ad.clr, trt = ad.trt, batch = ad.batch, type = 'linear model')
#' ad.p.adj <- ad.lm$p_adjusted
#' head(ad.lm$AIC)
#'
linear_regres <- function(data,
                          trt,
                          batch,
                          type = 'linear model',
                          p.adjust.method = 'fdr'){

  if(type == 'linear model'){
    model <- list()
    p <- c()
    R2 = adj.R2 = RMSE = RSE = AIC = BIC = data.frame(trt.only = NA, trt.batch = NA)
    for(i in 1:ncol(data)){
      res.lm0 <- lm(data[,i] ~ trt)
      res.lm <- lm(data[,i] ~ trt + batch)
      summary.res0 <- summary(res.lm0)
      summary.res <- summary(res.lm)

      model[[i]] <- res.lm
      p[i] <- summary.res$coefficients[2,4]

      R2[i, ] <- c(summary.res0$r.squared, summary.res$r.squared)
      adj.R2[i, ] <- c(summary.res0$adj.r.squared, summary.res$adj.r.squared)
      RMSE[i, ] <- c(rmse(res.lm0), rmse(res.lm))
      RSE[i, ] <- c(summary.res0$sigma, summary.res$sigma)
      AIC[i, ] <- c(AIC(res.lm0), AIC(res.lm))
      BIC[i, ] <- c(BIC(res.lm0), BIC(res.lm))
    }
    p.adj <- p.adjust(p, method = p.adjust.method)
    cond.R2 = marg.R2 = NA
    rownames(R2) = rownames(adj.R2) = colnames(data)
  }

  if(type == 'linear mixed model'){
    model <- list()
    p <- c()
    cond.R2 = marg.R2 = RMSE = RSE = AIC = BIC = data.frame(trt.only = NA, trt.batch = NA)
    for(i in 1:ncol(data)){
      res.lm0 <- lm(data[,i] ~ trt)
      res.lmm <- lmer(data[,i] ~ trt + (1|batch))
      summary.res0 <- summary(res.lm0)
      summary.res <- summary(res.lmm)

      model[[i]] <- res.lmm
      p[i] <- summary.res$coefficients[2,4]

      cond.R2[i, ] <- c(summary.res0$r.squared, as.numeric(r2(res.lmm)[1]))
      marg.R2[i, ] <- c(NA, as.numeric(r2(res.lmm)[2]))
      RMSE[i, ] <- c(rmse(res.lm0), rmse(res.lmm))
      RSE[i, ] <- c(summary.res0$sigma, summary.res$sigma)
      AIC[i, ] <- c(AIC(res.lm0), AIC(res.lmm))
      BIC[i, ] <- c(BIC(res.lm0), BIC(res.lmm))
    }
    p.adj <- p.adjust(p, method = p.adjust.method)
    R2 = adj.R2 = NA
    rownames(cond.R2) = rownames(marg.R2) = colnames(data)
  }
  names(model) = names(p) = names(p.adj) = rownames(RMSE) = rownames(RSE) = rownames(AIC) = rownames(BIC) = colnames(data)

  result = list(type = type,
                model = model,
                raw.p = p,
                adj.p = p.adj,
                p.adjust.method = p.adjust.method,
                R2 = R2,
                adj.R2 = adj.R2,
                cond.R2 = cond.R2,
                marg.R2 = marg.R2,
                RMSE = RMSE,
                RSE = RSE,
                AIC = AIC,
                BIC = BIC)

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


