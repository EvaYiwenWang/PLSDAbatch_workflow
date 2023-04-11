################################################################################
#' Data Simulation from Gaussian Distribution
#'
#' This function simulates data with treatment and batch effects
#' under different \code{batch x treatment designs} from Gaussian distribution
#' using component-based and multivariate strategy.
#'
#' @importFrom Rdpack reprompt
#' @importFrom stats rnorm runif
#'
#' @param mean.batch Numeric, the absolute mean of two batches respectively
#' in batch associated components.
#' @param sd.batch Numeric, the variability (standard deviation) of batch effects.
#' @param mean.trt Numeric, the absolute mean of two treatment groups
#' respectively in treatment associated components.
#' @param sd.trt Numeric, the variability (standard deviation) of
#' treatment effects.
#' @param noise Numeric, the variability (standard deviation) of
#' background noise.
#' @param N Integer, the number of samples.
#' @param p_total Integer, the number of variables.
#' @param p_trt_relevant Integer, the number of variables affected by
#' treatment effects.
#' @param p_bat_relevant Integer, the number of variables affected by
#' batch effects.
#' @param percentage_overlap_samples Numeric, in the range of \eqn{0}
#' to \eqn{1}, the percentage of samples in the first batch and first
#' treatment group.
#' @param percentage_overlap_variables Numeric, in the range of
#' \eqn{0} to \eqn{1}, the percentage of variables with both batch and
#' treatment effects given fixed number of variables with each effect.
#' @param seeds Integer. Default value is randomly selected from
#' \eqn{0} to \eqn{999}.
#'
#'
#' @return \code{simData_Gaussian} returns a list that contains
#' the following components:
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
#' ## Not run:
#' ## First example
#' ## Balanced batch x treatment design
#' simulation1 <- simData_Gaussian(mean.batch = 7,
#'                                 sd.batch = 8,
#'                                 mean.trt = 3,
#'                                 sd.trt = 2,
#'                                 N = 40,
#'                                 p_total = 300,
#'                                 p_trt_relevant = 60,
#'                                 p_bat_relevant = 150,
#'                                 percentage_overlap_samples = 0.5,
#'                                 percentage_overlap_variables = 0.5,
#'                                 seeds = 1)
#'
#' ## Second example
#' ## Unbalanced batch x treatment design
#' simulation2 <- simData_Gaussian(mean.batch = 7,
#'                                 sd.batch = 8,
#'                                 mean.trt = 3,
#'                                 sd.trt = 2,
#'                                 N = 40,
#'                                 p_total = 300,
#'                                 p_trt_relevant = 60,
#'                                 p_bat_relevant = 150,
#'                                 percentage_overlap_samples = 0.2,
#'                                 percentage_overlap_variables = 0.5,
#'                                 seeds = 1)
#' ## End(Not run)
#' @export
simData_Gaussian = function(mean.batch,
                            sd.batch = 1,
                            mean.trt,
                            sd.trt = 1,
                            noise = 0.2,
                            N,
                            p_total,
                            p_trt_relevant,
                            p_bat_relevant,
                            percentage_overlap_samples = 0.5,
                            percentage_overlap_variables = 0.5,
                            seeds = round(runif(1, 0, 999))){

  set.seed(seeds)
  ##
  n = N/2

  ### batch component
  ## variables for Group 1
  bCompG1 <- rnorm(n, mean.batch, sd.batch)
  ## variables for Group 2
  bCompG2 <- rnorm(n, -mean.batch, sd.batch)
  bComp_relevant <- c(bCompG1, bCompG2)

  ### treatment component
  ## variables for Group 1
  tCompG1 <- rnorm(n, mean.trt, sd.trt)
  ## variables for Group 2
  tCompG2 <- rnorm(n, -mean.trt, sd.trt)
  tComp_relevant <- c(tCompG1[seq_len(n*percentage_overlap_samples)],
                      tCompG2[seq_len(n*(1-percentage_overlap_samples))],
                      tCompG1[((n*percentage_overlap_samples)+1):n],
                      tCompG2[((n*(1-percentage_overlap_samples))+1):n])

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
  Y.trt <- c(rep("trt1", n*percentage_overlap_samples),
             rep("trt2", n*(1-percentage_overlap_samples)),
             rep("trt1", n*(1-percentage_overlap_samples)),
             rep("trt2", n*percentage_overlap_samples))

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

################################################################################
#' Data Simulation from Multivariate Negative Binomial Distribution
#'
#' This function simulates data with treatment and batch effects
#' under different \code{batch x treatment designs} from multivariate negative
#' binomial distribution using quantile-quantile transformation between
#' multivariate normal and negative binomial distribution. This function can
#' generate two or three batch groups.
#'
#'
#' @importFrom Rdpack reprompt
#' @importFrom stats rnorm runif model.matrix pnorm qnbinom
#' @importFrom mvtnorm rmvnorm
#'
#' @param batch.group Integer, the number of batch groups. Only valid for 2 or 3
#' batch groups. Default value is 2.
#' @param mean.batch Numeric, the mean of batch effect coefficients for all
#' related variables.
#' @param sd.batch Numeric, the variability (standard deviation) of batch effect
#' coefficients for all related variables.
#' @param mean.trt Numeric, the mean of treatment effect coefficients for all
#' related variables
#' @param sd.trt Numeric, the variability (standard deviation) of treatment
#' effect coefficients for all related variables.
#' @param mean.bg Numeric, the mean of background noise.
#' @param sd.bg Numeric, the variability (standard deviation) of
#' background noise.
#' @param N Integer, the number of samples.
#' @param p_total Integer, the number of variables.
#' @param p_trt_relevant Integer, the number of variables affected by
#' treatment effects.
#' @param p_bat_relevant Integer, the number of variables affected by
#' batch effects.
#' @param percentage_overlap_samples Numeric, in the range of \eqn{0}
#' to \eqn{1}, the percentage of samples from one batch group in one
#' treatment group.
#' @param percentage_overlap_variables Numeric, in the range of \eqn{0}
#' to \eqn{1}, the percentage of variables with both batch and treatment
#' effects given fixed number of variables with each effect.
#' @param data.cor Correlation matrix. The correlation between each pair
#' of variables.
#' @param disp Numeric. The dispersion parameter of negative binomial
#' distribution.
#' @param prob_zero Numeric, in the range of \eqn{0} to \eqn{1}. The probability
#' of each value to be 0.
#' @param seeds Integer. Default value is randomly selected from
#' \eqn{0} to \eqn{999}.
#'
#'
#'
#' @return \code{simData_mnegbinom} returns a list that contains
#' the following components:
#' \item{data}{The simulated data.}
#' \item{cleanData}{The ground-truth data without the batch effect.}
#' \item{data_mu_matrix}{The \eqn{ln} value of the mean to model
#' negative binomial distribution for the simulated data.}
#' \item{clean_mu_matrix}{The \eqn{ln} value of the mean to model
#' negative binomial distribution for the ground-truth data.}
#' \item{Y.trt}{The simulated outcome variable indicating the treatment
#' grouping information.}
#' \item{Y.bat}{The simulated outcome variable indicating the batch
#' grouping information.}
#' \item{true.trt}{The variables assigned with the treatment effect.}
#' \item{true.batch}{The variables assigned with the batch effect.}
#' \item{zero_index}{The indexes of cells containing the resulting counts from
#' the Multivariate Negative Binomial Distribution that are inflated as zero.}
#'
#' @author Yiwen Wang, Kim-Anh Lê Cao
#'
#'
#' @references
#' \insertRef{wang2020multivariate}{PLSDAbatch}
#'
#' @examples
#' ## Not run:
#' ## First example
#' ## To simulate the correlation structure
#' data.cor.res = corStruct(p = 300, zero_prob = 0.7)
#'
#' ## Balanced batch x treatment design
#' simulation <- simData_mnegbinom(mean.batch = 7,
#'                                 sd.batch = 8,
#'                                 mean.trt = 3,
#'                                 sd.trt = 2,
#'                                 mean.bg = 0,
#'                                 sd.bg = 0.2,
#'                                 N = 40,
#'                                 p_total = 300,
#'                                 p_trt_relevant = 100,
#'                                 p_bat_relevant = 200,
#'                                 percentage_overlap_samples = 0.5,
#'                                 percentage_overlap_variables = 0.5,
#'                                 data.cor = data.cor.res$data.cor,
#'                                 disp = 10,
#'                                 prob_zero = 0,
#'                                 seeds = 1)
#'
#' ## Second example
#' ## Unbalanced batch x treatment design
#' simulation <- simData_mnegbinom(mean.batch = 7,
#'                                 sd.batch = 8,
#'                                 mean.trt = 3,
#'                                 sd.trt = 2,
#'                                 mean.bg = 0,
#'                                 sd.bg = 0.2,
#'                                 N = 40,
#'                                 p_total = 300,
#'                                 p_trt_relevant = 100,
#'                                 p_bat_relevant = 200,
#'                                 percentage_overlap_samples = 0.2,
#'                                 percentage_overlap_variables = 0.5,
#'                                 data.cor = data.cor.res$data.cor,
#'                                 disp = 10,
#'                                 prob_zero = 0,
#'                                 seeds = 1)
#' ## End(Not run)
#'
#' @export
simData_mnegbinom = function(batch.group = 2,
                             mean.batch,
                             sd.batch = 1,
                             mean.trt,
                             sd.trt = 1,
                             mean.bg = 0,
                             sd.bg = 0.5,
                             N,
                             p_total,
                             p_trt_relevant,
                             p_bat_relevant,
                             percentage_overlap_samples = 0.5,
                             percentage_overlap_variables = 0.5,
                             data.cor,
                             disp = 10,
                             prob_zero = 0,
                             seeds = round(runif(1, 0, 999))){

  set.seed(seeds)

  ## check if the number of batch groups is other than 2 or 3
  if(batch.group != 2 & batch.group != 3){
    stop('The number of batch groups has to be 2 or 3.')
  }

  if(batch.group == 2){
    # batch info
    Y.batch <- rep(c("batch1", "batch2"), each = N/2)

    # treatment info
    Y.trt <- c(rep("trt1", (N/2)*percentage_overlap_samples),
               rep("trt2", (N/2)*(1-percentage_overlap_samples)),
               rep("trt1", (N/2)*(1-percentage_overlap_samples)),
               rep("trt2", (N/2)*percentage_overlap_samples))
  }

  if(batch.group == 3){
    # batch info
    Y.batch <- rep(c("batch1", "batch2", "batch3"), each = N/3)

    # treatment info
    Y.trt <- c(rep("trt1", (N/3)*percentage_overlap_samples),
               rep("trt2", (N/3)*(1-percentage_overlap_samples)),
               rep("trt1", (N/3)*(1-percentage_overlap_samples)),
               rep("trt2", (N/3)*percentage_overlap_samples),
               rep("trt1", (N/3)*percentage_overlap_samples),
               rep("trt2", (N/3)*(1-percentage_overlap_samples)))
  }

  # design matrix
  X <- model.matrix(~ as.factor(Y.trt) + as.factor(Y.batch))

  if(batch.group == 2){
    X <- X[,2:3] # no intercept

    # coefficients matrix
    beta = matrix(0, nrow = p_total, ncol = 2)
  }

  if(batch.group == 3){
    X <- X[,2:4] # no intercept

    # coefficients matrix
    beta = matrix(0, nrow = p_total, ncol = 3)
  }


  # treatment effects related variable coefficients
  ## p_trt_relevant variables with trt effects
  index.t <- 1:p_trt_relevant
  beta[index.t, 1] = rnorm(p_trt_relevant, mean = mean.trt, sd = sd.trt)

  # batch effects related variable coefficients
  ## p_bat_relevant variables with batch effects
  min.num <- min(p_trt_relevant, p_bat_relevant)
  set.seed(777)
  index.b <- c(sample(index.t, min.num*percentage_overlap_variables),
               sample((1:p_total)[-index.t],
                      (p_bat_relevant-min.num*percentage_overlap_variables)))

  set.seed(seeds)
  if(batch.group == 2){
    beta[index.b, 2] = rnorm(p_bat_relevant, mean = mean.batch, sd = sd.batch)

    # clean data coefficients
    beta_clean = beta
    beta_clean[, 2] = 0
  }

  if(batch.group == 3){
    beta[index.b[1:(p_bat_relevant/2)], 2] = rnorm(p_bat_relevant/2,
                                                   mean = mean.batch,
                                                   sd = sd.batch)
    beta[index.b[(p_bat_relevant/2+1):p_bat_relevant],
         3] = rnorm(p_bat_relevant/2,
                    mean = mean.batch,
                    sd = sd.batch)

    # clean data coefficients
    beta_clean = beta
    beta_clean[, 2:3] = 0
  }

  # build mu matrix
  ## simu data
  mu.simu = X%*%t(beta)

  ## clean data without batch effects
  mu.clean = X%*%t(beta_clean)

  # add noise
  noise.matrix <- matrix(rnorm(N*p_total, mean = mean.bg, sd = sd.bg), nrow = N)
  mu.simu <- mu.simu + noise.matrix
  mu.clean <- mu.clean + noise.matrix

  # exponentiation
  mu.simu.exp = exp(mu.simu)
  mu.clean.exp = exp(mu.clean)


  gama.simu = rep(disp, N)
  normd <- rmvnorm(n = N, mean = rep(0, p_total), sigma = data.cor)
  normcdf <- pnorm(normd)

  # simu data
  Y <- mapply(normcdf,
              mu.simu.exp,
              gama.simu,
              FUN = function(p, u, gamma){
                qnbinom(p = p, mu = u, size = gamma)
              })


  # clean data
  Y_clean <- mapply(normcdf,
                    mu.clean.exp,
                    gama.simu,
                    FUN = function(p, u, gamma){
                      qnbinom(p = p, mu = u, size = gamma)
                    })

  # inflate zero
  ## simu data
  Y_f <- Y
  zero_index = runif(length(Y_f)) < prob_zero
  Y_f[zero_index] = 0

  ## clean data
  Y_clean_f <- Y_clean
  Y_clean_f[zero_index] = 0

  Y_f = matrix(Y_f, nrow = N)
  Y_clean_f = matrix(Y_clean_f, nrow = N)

  return(list(data = Y_f,
              cleanData = Y_clean_f,
              data_mu_matrix = mu.simu,
              clean_mu_matrix = mu.clean,
              Y.trt = Y.trt,
              Y.batch = Y.batch,
              true.trt = index.t,
              true.batch = index.b,
              zero_index = zero_index))
}


################################################################################
#' Simulation of Correlation Structure for Microbiome Data
#'
#' This function simulates the correlation between variables
#' within the microbiome data.
#'
#'
#' @importFrom Rdpack reprompt
#' @importFrom stats cov2cor runif
#'
#'
#' @param p Integer, the number of variables to simulate.
#' @param zero_prob Numeric, in the range of \eqn{0} to \eqn{1}. The percentage
#' of zeroes in the non-diagonal lower-triangular matrix used to
#' generate the precision matrix.
#'
#' @return \code{corStruct} returns a list that contains
#' the following components:
#' \item{data.cor}{The correlation matrix.}
#' \item{data.precision.m}{The precision matrix, which is the
#' inverse of covariance matrix.}
#' \item{data.cov}{The covariance matrix.}
#'
#'
#' @author Yiwen Wang, Kim-Anh Lê Cao
#'
#'
#' @references
#' \insertRef{mcgregor2020mdine}{PLSDAbatch}
#'
#' @examples
#' data.cor.res = corStruct(p = 300, zero_prob = 0.7)
#'
#' @export
corStruct = function(p,
                     zero_prob = 0){

  unif1 <- runif(p*(p-1)/2, min = -1.5, max = 1.5)
  unif2 <- runif(p, min = 1.5, max = 2.5)

  # Randomly set lower-triangular elements to zero with probability 'zero_prob'.
  unif1_wzero = unif1
  unif1_wzero[runif(length(unif1)) < zero_prob] = 0

  data.prec.half <- matrix(0, nrow = p, ncol = p)
  data.prec.half[lower.tri(data.prec.half, diag = F)] = unif1_wzero
  diag(data.prec.half) <- unif2

  data.prec <- data.prec.half %*% t(data.prec.half)

  data.cov <- solve(data.prec)

  data.cor <- cov2cor(data.cov)

  rownames(data.cor) = colnames(data.cor) = paste0('OTU', 1:p)

  return(list(data.cor = data.cor,
              data.precision.m = data.prec,
              data.cov = data.cov))
}
