#' Data Simulation Based on Components and Multivariate Strategy
#'
#' This function simulates data with treatment or/and batch effects
#' under different \code{batch x treatment design}
#' using component-based and multivariate strategy.
#'
#' @importFrom Rdpack reprompt
#' @importFrom stats rnorm runif
#'
#' @param bat.mean Numeric, the absolute mean of two batches respectively
#' in batch associated components.
#' @param sd.bat Numeric, the variability of batch effects.
#' @param trt.mean Numeric, the absolute mean of two treatment groups
#' respectively in treatment associated components.
#' @param sd.trt Numeric, the variability of treatment effects.
#' @param noise Numeric, the variability (standard deviation) of
#' background noise.
#' @param N Integer, the number of samples.
#' @param p_total Integer, the number of variables.
#' @param p_trt_relevant Integer, the number of variables affected with
#' treatment effects.
#' @param p_bat_relevant Integer, the number of variables affected with
#' batch effects.
#' @param percentage_samples Numeric, in the range of \eqn{0} to \eqn{1},
#' the percentage of samples in the first batch and first treatment group.
#' @param percentage_overlap_variables Numeric, in the range of
#' \eqn{0} to \eqn{1},
#' the percentage of variables with both batch and treatment effects given
#' fixed number of variables with each effect.
#' @param seeds Integer. Default value is randomly selected from
#' \eqn{0} to \eqn{999}.
#'
#'
#' @return \code{simData} returns a list that contains the following components:
#' \item{data}{The simulated data.}
#' \item{cleanData}{The ground-truth data without the batch effect.}
#' \item{Y.trt}{The simulated outcome variable indicating the treatment
#' grouping information.}
#' \item{Y.bat}{The simulated outcome variable indicating the batch
#' grouping information.}
#' \item{true.trt}{The variables assigned with the treatment effect.}
#' \item{true.batch}{The variables assigned with the batch effect.}
#' \item{tComp_relevant}{The component simulated with the treatment effect.}
#' \item{bComp_relevant}{The component simulated with the batch effect.}
#'
#' @author Yiwen Wang, Kim-Anh Lê Cao
#'
#'
#' @references
#' \insertRef{wang2020multivariate}{PLSDAbatch}
#'
#' @examples
#' ## First example
#' ## Balanced batch x treatment design
#' i = 1
#' simulation1 <- simData(bat.mean = 7, sd.bat = 8, trt.mean = 3, sd.trt = 2,
#'                        N = 40, p_total = 300, p_trt_relevant = 60,
#'                        p_bat_relevant = 150,
#'                        percentage_samples = 0.5,
#'                        percentage_overlap_variables = 0.5)
#'
#' ## Second example
#' ## Unbalanced batch x treatment design
#' i = 1
#' simulation2 <- simData(bat.mean = 7, sd.bat = 8, trt.mean = 3, sd.trt = 2,
#'                        N = 40, p_total = 300, p_trt_relevant = 60,
#'                        p_bat_relevant = 150,
#'                        percentage_samples = 0.2,
#'                        percentage_overlap_variables = 0.5)
#'
#' @export
simData = function(bat.mean,
                   sd.bat = 1,
                   trt.mean,
                   sd.trt = 1,
                   noise = 0.2,
                   N,
                   p_total,
                   p_trt_relevant,
                   p_bat_relevant,
                   percentage_samples = 0.5,
                   percentage_overlap_variables = 0.5,
                   seeds = round(runif(1, 0, 999))){

  set.seed(seeds)
  ##
  n = N/2

  ### batch component
  ## variables for Group 1
  bCompG1 <- rnorm(n, bat.mean, sd.bat)
  ## variables for Group 2
  bCompG2 <- rnorm(n, -bat.mean, sd.bat)
  bComp_relevant <- c(bCompG1, bCompG2)

  ### treatment component
  ## variables for Group 1
  tCompG1 <- rnorm(n, trt.mean, sd.trt)
  ## variables for Group 2
  tCompG2 <- rnorm(n, -trt.mean, sd.trt)
  tComp_relevant <- c(tCompG1[seq_len(n*percentage_samples)],
                      tCompG2[seq_len(n*(1-percentage_samples))],
                      tCompG1[((n*percentage_samples)+1):n],
                      tCompG2[((n*(1-percentage_samples))+1):n])

  ## data noise
  # background noise is independent
  Data.noise <- matrix(rnorm(N*p_total, mean = 0, sd = noise), nrow = N)


  ### treatment loadings
  ## p_trt_relevant variables with trt effects
  set.seed(777) # fix loadings for each repeat
  w.relevant.t <- sample(c(runif(p_trt_relevant/2,-0.3,-0.2),
                           runif(p_trt_relevant/2,0.2,0.3)))
  w.relevant.t <- w.relevant.t/sqrt(sum(w.relevant.t^2))
  Dat.relevant.t <- matrix(tComp_relevant, ncol = 1) %*%
    matrix(w.relevant.t, nrow = 1)


  ### batch loadings
  ## p_bat_relevant variables with batch effects
  w.relevant.b <- sample(c(runif(p_bat_relevant/2,-0.3,-0.2),
                           runif(p_bat_relevant/2,0.2,0.3)))
  w.relevant.b <- w.relevant.b/sqrt(sum(w.relevant.b^2))
  Dat.relevant.b <- matrix(bComp_relevant, ncol = 1) %*%
    matrix(w.relevant.b, nrow = 1)


  ## final dataset
  ## the index of affected variables with trt effects
  index.t <- seq_len(p_trt_relevant)

  ## the largest overlap of variables with both batch and trt effects
  overlap.num <- min(p_trt_relevant, p_bat_relevant)

  ## the index of affected variables with batch effects
  index.b <- c(sample(index.t, overlap.num*percentage_overlap_variables),
               sample((seq_len(p_total))[-index.t],
                      (p_bat_relevant-overlap.num*
                         percentage_overlap_variables)))

  finalData <- Data.noise
  finalData[,index.t] <- finalData[,index.t] + Dat.relevant.t
  finalData[,index.b] <- finalData[,index.b] + Dat.relevant.b

  ## the ground-truth data
  cleanData <- Data.noise
  cleanData[,index.t] <- cleanData[,index.t] + Dat.relevant.t


  rownames(finalData) <- rownames(cleanData) <- paste0('sample', seq_len(N))


  ## outcome variables
  Y.bat <- rep(c("bat1", "bat2"), each = n)
  Y.trt <- c(rep("trt1", n*percentage_samples), rep("trt2",
                                                    n*(1-percentage_samples)),
             rep("trt1", n*(1-percentage_samples)), rep("trt2",
                                                        n*percentage_samples))

  names(Y.trt) <- names(Y.bat) <- rownames(finalData)


  return(list(data = finalData,
              cleanData = cleanData,
              Y.trt = Y.trt,
              Y.bat = Y.bat,
              true.trt = index.t,
              true.batch = index.b,
              tComp_relevant = tComp_relevant,
              bComp_relevant = bComp_relevant))
}
