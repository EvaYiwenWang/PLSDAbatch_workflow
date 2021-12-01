######################################
#' Partitioned Variance Plot
#'
#' This function draws a partitioned variance plot explained
#' by different sources.
#'
#' @import ggplot2
#'
#' @param prop.df A data frame that contains the proportion of variance
#' explained by different sources.
#' @param text.cex Numeric, the size of text on the plot.
#' @param x.angle Numeric, angle of x axis, in the range of
#' \eqn{0} to \eqn{360}.
#' @param x.hjust Numeric, horizontal justification of x axis,
#' in the range of \eqn{0} to \eqn{1}.
#' @param title Character, the plot title.
#' @param color.set A vector of characters, indicating the set of colors to use.
#' The colors are represented by hexadecimal color code.
#'
#' @return None.
#'
#' @author Yiwen Wang, Kim-Anh Lê Cao
#'
#' @seealso \code{\link{Scatter_Density}}, \code{\link{box_plot}},
#' \code{\link{density_plot}} and \code{\link{alignment_score}} as the other
#' methods for batch effect detection and batch effect removal assessment.
#'
#' @export
#'
#' @examples
#' ## First example
#' library(vegan) # for function varpart()
#' data('AD_data')
#' ad.clr <- AD_data$EgData$X.clr # centered log ratio transformed data
#' ad.batch = AD_data$EgData$Y.bat # batch information
#' ad.trt = AD_data$EgData$Y.trt # treatment information
#'
#' ad.factors.df <- data.frame(trt = ad.trt, batch = ad.batch)
#' rda.res <- varpart(ad.clr, ~ trt, ~ batch,
#'                    data = ad.factors.df, scale = TRUE)
#'
#' ad.prop.df <- data.frame(Treatment = NA, Intersection = NA,
#'                          Batch = NA, Residuals = NA)
#' ad.prop.df[1,] <- rda.res$part$indfract$Adj.R.squared
#'
#' ad.prop.df[ad.prop.df < 0] = 0
#' ad.prop.df <- as.data.frame(t(apply(ad.prop.df, 1, function(x){x/sum(x)})))
#'
#' partVar_plot(prop.df = ad.prop.df)
#'
#' ## Second example
#' # a list of data corrected from different methods
#' ad.corrected.list <- AD_data$CorrectData
#' ad.prop.df <- data.frame(Treatment = NA, Intersection = NA,
#'                          Batch = NA, Residuals = NA)
#' for(i in seq_len(length(ad.corrected.list))){
#'   rda.res = varpart(ad.corrected.list[[i]], ~ trt, ~ batch,
#'                     data = ad.factors.df, scale = TRUE)
#'   ad.prop.df[i, ] <- rda.res$part$indfract$Adj.R.squared}
#'
#' rownames(ad.prop.df) = names(ad.corrected.list)

#' ad.prop.df[ad.prop.df < 0] = 0
#' ad.prop.df <- as.data.frame(t(apply(ad.prop.df, 1,
#'                                     function(x){x/sum(x)})))
#'
#' partVar_plot(prop.df = ad.prop.df)
#'
#'
partVar_plot <- function(prop.df,
                         text.cex = 3,
                         x.angle = 60,
                         x.hjust = 1,
                         title = NULL,
                         color.set = NULL){

  Prop = Methods = Type = ypos = NULL
  rda.ggplot <- data.frame(Prop = c(t(prop.df)),
                           Methods = rep(rownames(prop.df),
                                         each = ncol(prop.df)),
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

  # color set
  if(is.null(color.set)){
    color.set = pb_color(11:20)
  }else{
    color.set = color.set
  }


  ggplot(rda.ggplot, aes(x = Methods, y = Prop, fill = Type)) +
    geom_bar(stat = 'identity') + ylab('Explained variance (%)') +
    scale_fill_manual(name = 'Variation sources', values = color.set) +
    theme_bw() + theme(axis.text.x = element_text(angle = x.angle,
                                                  hjust = x.hjust)) +
    geom_text(aes(y = ypos, label = round(Prop, digits = 3)),
              vjust = 1.6, color = 'black', size = text.cex) +
    labs(title = title)

}




