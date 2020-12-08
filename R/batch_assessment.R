######################################
#' Partitioned Variance Plot
#'
#' This function draws a partitioned variance plot explained by different sources.
#'
#' @param prop.df A data frame that contains the proportion of variance explained by different sources.
#' @param text.cex Numeric, the size of text on the plot.
#' @param x.angle Numeric, angle of x axis, in the range of \eqn{0} to \eqn{360}.
#' @param x.hjust Numeric, horizontal justification of x axis, in the range of \eqn{0} to \eqn{1}.
#' @param title Character, the plot title.
#'
#' @return none
#' @export
#'
#' @examples
#' ## First example
#' library(vegan) # for function varpart()
#' data('AD_data')
#' ad.clr <- AD_data$EgData$X.clr
#' ad.clr.scale <- scale(ad.clr, center = T, scale = T)
#' ad.batch = AD_data$EgData$Y.bat
#' ad.trt = AD_data$EgData$Y.trt
#'
#' ad.factors.df <- data.frame(trt = ad.trt, batch = ad.batch)
#' rda.res <- varpart(ad.clr.scale, ~ trt, ~ batch, data = ad.factors.df)
#'
#' ad.prop.df <- data.frame(Treatment = NA, Intersection = NA, Batch = NA, Residuals = NA)
#' ad.prop.df[1,] <- rda.res$part$indfract$Adj.R.squared
#'
#' ad.prop.df[ad.prop.df < 0] = 0
#' ad.prop.df <- as.data.frame(t(apply(ad.prop.df, 1, function(x){x/sum(x)})))
#'
#' partVar_plot(prop.df = ad.prop.df)
#'
#' ## Second example
#' ad.corrected.list <- AD_data$CorrectData
#' ad.corr_scale.list <- list()
#' for(i in 1:length(ad.corrected.list)){
#'   # scale the data on OTUs
#'   scale.res <- scale(ad.corrected.list[[i]], center = T, scale = T)
#'   ad.corr_scale.list[[i]] <- scale.res
#' }
#' names(ad.corr_scale.list) <- names(ad.corrected.list)
#'
#' ad.prop.df <- data.frame(Treatment = NA, Intersection = NA, Batch = NA, Residuals = NA)
#' for(i in 1:length(ad.corr_scale.list)){
#'   rda.res = varpart(ad.corr_scale.list[[i]], ~ trt, ~ batch,
#'                     data = ad.factors.df)
#'   ad.prop.df[i, ] <- rda.res$part$indfract$Adj.R.squared}
#'
#' rownames(ad.prop.df) = names(ad.corr_scale.list)
#' ad.prop.df[ad.prop.df < 0] = 0
#' ad.prop.df <- as.data.frame(t(apply(ad.prop.df, 1, function(x){x/sum(x)})))
#' partVar_plot(prop.df = ad.prop.df)
#'
#'
partVar_plot <- function(prop.df,
                         text.cex = 3,
                         x.angle = 60,
                         x.hjust = 1,
                         title = NULL){

  rda.ggplot <- data.frame(Prop = c(t(prop.df)),
                           Methods = rep(rownames(prop.df), each = ncol(prop.df)),
                           Type = rep(colnames(prop.df), nrow(prop.df)))

  rda.ggplot$Methods <- factor(rda.ggplot$Methods,
                               levels = rownames(prop.df))

  rda.ggplot$Type <- factor(rda.ggplot$Type,
                            levels = rev(colnames(prop.df)))


  rda.ggplot.position <- as.matrix(prop.df)
  rda.ggplot.position[which(rda.ggplot.position <= 0.03)] = 0.03

  rda.ggplot.position <- t(apply(rda.ggplot.position, 1, cumsum))
  rda.ggplot.position[,1] <- 0.06
  rda.ggplot.position[,ncol(prop.df)] <- 1


  rda.ggplot$ypos <- c(t(rda.ggplot.position))

  ggplot(rda.ggplot, aes(x = Methods, y = Prop, fill = Type)) +
    geom_bar(stat = 'identity') + ylab('Explained variance (%)') +
    scale_fill_manual(name = 'Variation sources', values = pb_color(11:14)) +
    theme_bw() + theme(axis.text.x = element_text(angle = x.angle, hjust = x.hjust)) +
    geom_text(aes(y = ypos, label = round(Prop, digits = 3)),
              vjust = 1.6, color = 'black', size = text.cex) + labs(title = title)

}




