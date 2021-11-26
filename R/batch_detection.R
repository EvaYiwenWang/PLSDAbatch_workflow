##########################################################################
#' Principal Component Analysis (PCA) with Density Plots per Component
#'
#' This function draws a PCA sample plot with density plots per
#' principal component.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom gridExtra grid.arrange
#' @importFrom ggpubr get_legend
#'
#' @param object The object of class PCA.
#' @param batch A factor or a class vector for the batch grouping information
#' (categorical outcome variable).
#' @param trt A factor or a class vector for the treatment grouping information
#' (categorical outcome variable).
#' @param xlim A numeric vector of length 2, indicating the x coordinate ranges.
#' @param ylim A numeric vector of length 2, indicating the y coordinate ranges.
#' @param color.set A vector of character, indicating the set of colors to use.
#' The colors are represented by hexadecimal color code.
#' @param batch.legend.title Character, the legend title of batches.
#' @param trt.legend.title Character, the legend title of treatments.
#' @param density.lwd Numeric, the thickness of density lines.
#' @param title Character, the plot title.
#' @param title.cex Numeric, the size of plot title.
#' @param legend.cex Numeric, the size of legends.
#' @param legend.title.cex Numeric, the size of legend title.
#'
#' @return none
#'
#' @author Yiwen Wang, Kim-Anh Lê Cao
#'
#' @seealso \code{\link{box_plot}}, \code{\link{density_plot}},
#' \code{\link{alignment_score}} and \code{\link{partVar_plot}} as the other
#' methods for batch effect detection and batch effect removal assessment.
#'
#' @export
#'
#' @examples
#' # The first example
#' library(mixOmics) # for function pca()
#' data('AD_data')
#' ad.clr <- AD_data$EgData$X.clr # centered log ratio transformed data
#' ad.pca.before <- pca(ad.clr, ncomp = 3, scale = TRUE)
#' ad.batch = AD_data$EgData$Y.bat # batch information
#' ad.trt = AD_data$EgData$Y.trt # treatment information
#' Scatter_Density(object = ad.pca.before, batch = ad.batch, trt = ad.trt)
#'
#' # The second example
#' colorlist <- rainbow(10)
#' Scatter_Density(object = ad.pca.before, batch = ad.batch, trt = ad.trt,
#'                 color.set = colorlist)
#'
Scatter_Density <- function(object,
                            batch = NULL, trt = NULL,
                            xlim = NULL, ylim = NULL,
                            color.set = NULL,
                            batch.legend.title = 'Batch',
                            trt.legend.title = 'Treatment',
                            density.lwd = 0.2,
                            title = NULL, title.cex = 1.5,
                            legend.cex = 0.7, legend.title.cex = 0.75){
  data = as.data.frame(object[['variates']][['X']])
  expl.var = object[['explained_variance']]

  if(is.null(batch)){batch <- batch}else{batch <- as.factor(batch)}
  if(is.null(trt)){trt <- trt}else{trt <- as.factor(trt)}

  # color set
  if(is.null(color.set)){
    color.set = color.mixo(seq_len(10))
  }else{
    color.set = color.set
  }

  # main plot
  pMain <- ggplot(data = data, aes(x = data[ ,1], y = data[ ,2], colour = batch,
                                   shape = trt)) +
    geom_point() + xlab(paste0('PC1: ', round(as.numeric(expl.var[1])*100),
                               '% expl.var')) +
    ylab(paste0('PC2: ', round(as.numeric(expl.var[2])*100), '% expl.var')) +
    scale_shape_manual(values = c(1,19,2,17,4)) +
    scale_color_manual(values = color.set) + theme_bw() +
    labs(colour = batch.legend.title, shape = trt.legend.title) +
    scale_x_continuous(limits = xlim) + scale_y_continuous(limits = ylim) +
    theme(legend.position = 'right', legend.box = 'horizontal',
          legend.direction = 'vertical',
          legend.key.height = unit(0.2, 'cm'),
          legend.key.width = unit(0.1, 'cm'),
          legend.title = element_text(size = rel(legend.title.cex)),
          legend.spacing.x = unit(0.1, 'cm'),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.text = element_text(size = rel(legend.cex)))

  xlim.update <- layer_scales(pMain)$x$get_limits()
  ylim.update <- layer_scales(pMain)$y$get_limits()

  # top density plot
  pTop <- ggplot(data = data, aes(x = data[ ,1], fill = batch,
                                  linetype = trt)) +
    geom_density(size = density.lwd, alpha = 0.5) + ylab('Density') +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = rel(0.8)),
          plot.title = element_text(hjust = 0.5, size = rel(title.cex)),
          axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(), legend.position = 'none') +
    scale_fill_manual(values = color.set) +
    scale_x_continuous(limits = xlim.update) + labs(title = title)

  # right density plot
  pRight <- ggplot(data = data, aes(x = data[ ,2],
                                    fill = batch, linetype = trt)) +
    geom_density(size = density.lwd, alpha = 0.5) +  coord_flip() +
    ylab('Density') +
    theme(axis.title.x = element_text(size = rel(0.8)),
          axis.title.y = element_blank(), axis.line = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.background = element_blank(), legend.position = 'none') +
    scale_fill_manual(values = color.set) +
    scale_x_continuous(limits = ylim.update)

  if(is.null(batch) && is.null(trt)){legend <-
    grid.rect(gp = gpar(col="white"))}else{
    legend <- get_legend(pMain)
  }

  grid.arrange(pTop, legend, pMain + theme(legend.position = 'none'), pRight,
               ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))

}


######################################
#' Box Plot
#'
#' This function draws side-by-side box plots for each batch.
#'
#' @importFrom ggplot2 ggplot
#' @param df A data frame used to draw the box plots.
#' @param title Character, the plot title.
#' @param batch.legend.title Character, the legend title of batches.
#' @param ylab Character, y-axis title.
#' @param color.set A vector of character, indicating the set of colors to use.
#' The colors are represented by hexadecimal color code.
#' @param x.angle Numeric, angle of x axis, in the range of
#' \eqn{0} to \eqn{360}.
#' @param x.hjust Numeric, horizontal justification of x axis, in the range of
#' \eqn{0} to \eqn{1}.
#' @param x.vjust Numeric, vertical justification of x axis, in the range of
#' \eqn{0} to \eqn{1}.
#'
#' @return none
#'
#' @author Yiwen Wang, Kim-Anh Lê Cao
#'
#' @seealso \code{\link{Scatter_Density}}, \code{\link{density_plot}},
#' \code{\link{alignment_score}} and \code{\link{partVar_plot}} as the other
#' methods for batch effect detection and batch effect removal assessment.
#'
#' @export
#'
#' @examples
#' # The first example
#' data('AD_data')
#' ad.clr <- AD_data$EgData$X.clr # centered log ratio transformed data
#' ad.batch = AD_data$EgData$Y.bat # batch information
#' ad.df <- data.frame(value = ad.clr[,1], batch = ad.batch)
#' box_plot(df = ad.df, title = 'OTU 12', x.angle = 30)
#'
#' # The second example
#' colorlist <- rainbow(10)
#' box_plot(df = ad.df, title = 'OTU 12', color.set = colorlist, x.angle = 30)
#'
box_plot <- function(df, title = NULL,
                     batch.legend.title = 'Batch',
                     ylab = 'Value',
                     color.set = NULL,
                     x.angle = 0, x.hjust = 0.5, x.vjust = 0.5){
  value <- df[,1]
  batch <- df[,2]

  # color set
  if(is.null(color.set)){
    color.set = color.mixo(seq_len(10))
  }else{
    color.set = color.set
  }


  ggplot(data = df, aes(x = batch, y = value, fill = batch)) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot() + scale_fill_manual(values = color.set) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = x.angle, hjust = x.hjust,
                                     vjust = x.vjust),
          panel.grid = element_blank(),
          axis.title.x = element_blank(), axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5,size = rel(1.2))) +
    labs(fill = batch.legend.title, y = ylab, title = title)
}


################################
#' Density Plot
#'
#' This function draws an overlap of multiple density plots for each batch.
#'
#' @importFrom ggplot2 ggplot
#' @param df A data frame used to draw the density plots.
#' @param title Character, the plot title.
#' @param batch.legend.title Character, the legend title of batches.
#' @param xlab Character, x-axis title.
#' @param color.set A vector of character, indicating the set of colors to use.
#' The colors are represented by hexadecimal color code.
#' @param title.hjust Numeric, horizontal justification of the plot title,
#' in the range of \eqn{0} to \eqn{1}.
#'
#' @return none
#'
#' @author Yiwen Wang, Kim-Anh Lê Cao
#'
#' @seealso \code{\link{Scatter_Density}}, \code{\link{box_plot}},
#' \code{\link{alignment_score}} and \code{\link{partVar_plot}} as the other
#' methods for batch effect detection and batch effect removal assessment.
#'
#' @export
#'
#' @examples
#' # The first example
#' data('AD_data')
#' ad.clr <- AD_data$EgData$X.clr # centered log ratio transformed data
#' ad.batch = AD_data$EgData$Y.bat # batch information
#' ad.df <- data.frame(value = ad.clr[,1], batch = ad.batch)
#' density_plot(df = ad.df, title = 'OTU 12')
#'
#' # The second example
#' colorlist <- rainbow(10)
#' density_plot(df = ad.df, title = 'OTU 12', color.set = colorlist)
#'
density_plot <- function(df, title = NULL,
                         batch.legend.title = 'Batch',
                         xlab = 'Value',
                         color.set = NULL,
                         title.hjust = 0.5){
  value <- df[,1]
  batch <- df[,2]

  # color set
  if(is.null(color.set)){
    color.set = color.mixo(seq_len(10))
  }else{
    color.set = color.set
  }

  ggplot(data = df, aes(x = value, fill =  batch)) + geom_density(alpha = 0.5) +
    scale_fill_manual(values = color.set) +
    labs(title = title, x = xlab, y = 'Density', fill = batch.legend.title) +
    theme_bw() + theme(plot.title = element_text(hjust = title.hjust),
                       panel.grid = element_blank())
}



