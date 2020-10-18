
#require(mixOmics)
##############
#' Matrix Deflation
#'
#' This function removes the variance of given component \code{t} from the input matrix \code{X}.
#' \deqn{\hat{X} = X - t (t^{\top}t)^{-1}t^{\top}X}
#' It is a built-in funciton of \code{PLSDA_batch}.
#'
#' @param X A numeric matrix to be deflated. It assumes that samples are on the row,
#' while variables are on the column. \code{NA}s are not allowed.
#' @param t A component to be deflated out from the matrix.
#'
#' @return A deflated matrix with the same dimension as the input matrix.
#'
#' @examples
#' A built-in funciton of PLSDA_batch, not separately used.
#'
#'
#' @export
mtx_deflation <- function(X, t){
  X.res = X - t %*% (solve(crossprod(t))) %*% (t(t) %*% X)
  return(invisible(X.res))
}

###############
#' Partial Least Squares Discriminant Analysis
#'
#' This function estimates latent dimensions from the explanatory matrix \code{X}.
#' The latent dimensions are maximally associated with the outcome matrix \code{Y}.
#' It is a built-in funciton of \code{PLSDA_batch}.
#'
#' @param X A numeric matrix that is centered and scaled as an explanatory matrix. \code{NA}s are not allowed.
#' @param Y A dummy matrix that is centered and scaled as an outcome matrix.
#' @param ncomp Integer, the number of dimensions to include in the model.
#' @param keepX A numeric vector of length \code{ncomp}, the number of variables to keep in \eqn{X}-loadings.
#' By default all variables are kept in the model. A valid input of \code{keepX} extends \code{PLSDA} to a sparse version.
#' @param tol Numeric, convergence stopping value.
#' @param max.iter Integer, the maximum number of iterations.
#'
#' @return \code{PLSDA} returns a list that contains the following components:
#'
#' \item{original_data}{The original explanatory matrix \code{X} and outcome matrix \code{Y}.}
#' \item{defl_data}{The centered and scaled deflated matrices (\eqn{\hat{X}} and \eqn{\hat{Y}})
#' after removing the variance of latent components calculated with estimated latent dimensions.}
#' \item{latent_var}{The latent components calculated with estimated latent dimensions.}
#' \item{loadings}{The estimated latent dimensions.}
#' \item{iters}{Number of iterations of the algorthm for each component.}
#' \item{exp_var}{The amount of data variance explained per component (note that contrary to \code{PCA},
#' this amount may not decrease as the aim of the method is not to maximise the variance,
#' but the covariance between \code{X} and the dummy matrix \code{Y}).}
#'
#'
#' @examples
#' A built-in funciton of PLSDA_batch, not separately used.
#'
#' @export
PLSDA <- function(X, Y, ncomp, keepX = rep(ncol(X), ncomp), tol = 1e-06, max.iter = 500){
  # Y is dummy matrix
  # input X, Y are scaled
  # don't consider NA value

  # require(mixOmics)

  mat.t = matrix(nrow = nrow(X), ncol = ncomp)
  mat.u = matrix(nrow = nrow(X), ncol = ncomp)
  mat.a = matrix(nrow = ncol(X), ncol = ncomp)
  mat.b = matrix(nrow = ncol(Y), ncol = ncomp)

  c.iter = NULL
  X.temp = X
  Y.temp = Y
  for(h in 1:ncomp){
    nx = ncol(X) - keepX[h]

    # Initialisation
    M = crossprod(X.temp, Y.temp)
    svd.M = svd(M, nu = 1, nv = 1)
    a.old = svd.M$u
    b.old = svd.M$v

    t = X.temp %*% a.old
    u = Y.temp %*% b.old

    iter = 1

    # Iteration
    repeat{
      a.new = t(X.temp) %*%  u

		  if (nx != 0){
        abs_a = abs(a.new)
        if(any(rank(abs_a, ties.method = "max") <= nx)){
          a.new = ifelse(rank(abs_a, ties.method = "max") <= nx, 0,
            sign(a.new) * (abs_a - max(abs_a[rank(abs_a, ties.method = "max") <= nx])))
        }
      }

      a.new = a.new / drop(sqrt(crossprod(a.new)))

      t = X.temp %*% a.new

      b.new = t(Y.temp) %*% t

      b.new = b.new / drop(sqrt(crossprod(b.new)))

      u = Y.temp %*% b.new

      if (crossprod(a.new - a.old) < tol) {break}
      if (iter == max.iter){
        warning(paste("Maximum number of iterations reached for the component", h), call. = FALSE)
        break
      }

      a.old = a.new
      b.old = b.new
      iter = iter + 1
    }

    # deflation
    X.temp <- mtx_deflation(X.temp, t)
    Y.temp <- mtx_deflation(Y.temp, u)

    mat.t[,h] = t
    mat.u[,h] = u
    mat.a[,h] = a.new
    mat.b[,h] = b.new
    c.iter[h] = iter
  }

  rownames(mat.t) = rownames(mat.u) = rownames(X)
  rownames(mat.a) = colnames(X)
  rownames(mat.b) = colnames(Y)
  colnames(mat.t) = colnames(mat.u) = colnames(mat.a) = colnames(mat.b) = names(c.iter) = paste('comp', 1:ncomp)

  exp.var.X = mixOmics::explained_variance(X, mat.t, ncomp = ncomp)
  exp.var.Y = mixOmics::explained_variance(Y, mat.u, ncomp = ncomp)

  result = list(original_data = list(X = X, Y = Y),
                defl_data = list(X = X.temp, Y = Y.temp),
                latent_comp = list(t = mat.t, u = mat.u),
                loadings = list(a = mat.a, b = mat.b),
                iters = c.iter,
                exp_var = list(X = exp.var.X, Y = exp.var.Y))
  return(invisible(result))
}

########################
#' Partial Least Squares Discriminant Analysis for Batch Effect Correction
#'
#' This function removes batch variation from the input data given batch grouping information
#' and the number of associated components.
#'
#' @param X A numeric matrix as an explanatory matrix. \code{NA}s are not allowed.
#' @param Y.trt A factor or a class vector for the treatment grouping information (categorical outcome variable).
#' Without the input of \code{Y.trt}, treatment variation cannot be preserved before correcting for batch effects.
#' @param Y.bat A factor or a class vector for the batch grouping information (categorical outcome variable).
#' @param ncomp.trt Integer, the number of treatment associated dimensions to include in the model.
#' @param ncomp.bat Integer, the number of batch associated dimensions to include in the model.
#' @param keepX.trt A numeric vector of length \code{ncomp.trt}, the number of variables to keep in \eqn{X}-loadings.
#' By default all variables are kept in the model. A valid input of \code{keepX.trt} extends \code{PLSDA_batch} to a sparse version.
#' @param keepX.bat A numeric vector of length \code{ncomp.bat}, the number of variables to keep in \eqn{X}-loadings.
#' By default all variables are kept in the model. We usually use default setting.
#' @param max.iter Integer, the maximum number of iterations.
#' @param tol Numeric, convergence stopping value.
#' @param near.zero.var Logical, should be set to \code{TRUE} in particular for data with many zero values.
#' Setting this argument to \code{FALSE} (when appropriate) will speed up the computations. Default value is \code{TRUE}.
#' @param balance Logical, should be set to \code{TRUE}, if the \code{batch x treatment design} is balanced (or complete).
#' Setting this argument to \code{FALSE} extends \code{PLSDA_batch} to \code{weighted PLSDA_batch}. Default value is \code{TRUE}.
#'
#' @return \code{PLSDA_batch} returns a list that contains the following components:
#'
#' \item{X}{The original explanatory matrix \code{X}.}
#' \item{X.nobatch}{The batch corrected matrix with the same dimension as the input matrix.}
#' \item{X.notrt}{The matrix from which treatment variation is removed.}
#' \item{Y}{The original outcome variables \code{Y.trt} and \code{Y.bat}.}
#' \item{latent_var.trt}{The treatment associated latent components calculated with corresponding latent dimensions.}
#' \item{latent_var.bat}{The batch associated latent components calculated with corresponding latent dimensions.}
#' \item{loadings.trt}{The estimated treatment associated latent dimensions.}
#' \item{loadings.bat}{The estimated batch associated latent dimensions.}
#' \item{tol}{The tolerance used in the iterative algorithm, convergence stopping value.}
#' \item{max.iter}{The maximum number of iterations.}
#' \item{iter.trt}{Number of iterations of the algorthm for each treatment associated component.}
#' \item{iter.bat}{Number of iterations of the algorthm for each batch associated component.}
#' \item{explained_variance.trt}{The amount of data variance explained per treatment associated component.}
#' \item{explained_variance.bat}{The amount of data variance explained per batch associated component.}
#' \item{weight}{The sample weights, all \eqn{1} for a balanced \code{batch x treatment design}.}
#'
#'
#' @examples
#' ## First example
#' data('AD_data')
#' X = AD_data$X.clr #centered log ratio transformed data
#' Y.trt = AD_data$Y.trt
#' Y.bat = AD_data$Y.bat
#' ad_plsda_batch <- PLSDA_batch(X, Y.trt, Y.bat, ncomp.trt = 1, ncomp.bat = 5)
#' ad_X.corrected <- ad_plsda_batch$X.nobatch #batch corrected data
#'
#' ## Second example
#' ## sparse PLSDA-batch
#' ad_splsda_batch <- PLSDA_batch(X, Y.trt, Y.bat, ncomp.trt = 1, keepX.trt = 30, ncomp.bat = 5)
#'
#' ## Third example
#' ## weighted PLSDA-batch
#' ## The simulated data with an unbalanced batch x treatment design
#' i = 1
#' simulation <- simData(bat.mean = 7, sd.bat = 8, trt.mean = 3, sd.trt = 2,
#'                        N = 40, p_total = 300, p_trt_relevant = 60, p_bat_relevant = 150,
#'                        percentage_samples = 0.2, percentage_overlap_variables = 0.5)
#'
#' sdata_wplsda_batch <- PLSDA_batch(simulation$data, simulation$Y.trt, simulation$Y.bat,
#'                                   ncomp.trt = 1, ncomp.bat = 1, balance = F)
#'
#' ## Fourth example
#' ## sparse weighted PLSDA-batch
#' sdata_swplsda_batch <- PLSDA_batch(simulation$data, simulation$Y.trt, simulation$Y.bat,
#'                                    ncomp.trt = 1, keepX.trt = length(simulation$true.trt),
#'                                    ncomp.bat = 1, balance = F)
#'
#'
#' @export
PLSDA_batch <- function(X,
                        Y.trt = NULL,
                        Y.bat,
                        ncomp.trt = 2,
                        ncomp.bat = 2,
                        keepX.trt = rep(ncol(X), ncomp.trt),
                        keepX.bat = rep(ncol(X), ncomp.bat),
                        max.iter = 500,
                        tol = 1e-06,
                        near.zero.var = TRUE,
                        balance = TRUE){

  # balance = T: the design of batch-treatment is balanced

  #-- validation of arguments --#
  if (length(dim(X)) != 2){
    stop("'X' must be a numeric matrix.")}

  X = as.matrix(X)

  if(!is.numeric(X)){
    stop("'X' must be a numeric matrix.")}

  n = nrow(X)

  if(!is.null(dim(Y.bat))){
    stop("'Y.bat' should be a factor or a class vector.")}

  Y.bat = as.factor(Y.bat)
  Y.bat.mat = mixOmics::unmap(as.numeric(Y.bat))
  Y.bat.mat = as.matrix(Y.bat.mat)
  q.bat = ncol(Y.bat.mat)

  if(n != nrow(Y.bat.mat)){
    stop("unequal number of rows in 'X' and 'Y.bat'.")}
  if(is.null(ncomp.bat) || !is.numeric(ncomp.bat) || ncomp.bat <= 0){
    stop("invalid number of variates, 'ncomp.bat'.")}

  ncomp.bat = round(ncomp.bat)

  ###################
  if(near.zero.var == TRUE){
    nzv = mixOmics::nearZeroVar(X)

    if (length(nzv$Position > 0)){
      warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
      X = X[, -nzv$Position,drop=FALSE]

      if(ncol(X)==0){
        stop("No more predictors")}
    }
  }

  p = ncol(X)

  ##################

  if(ncomp.bat > p){
    warning("Reset maximum number of variates 'ncomp.bat' to ncol(X) = ", p, ".")
    ncomp.bat = p}


  if (length(keepX.bat) != ncomp.bat){
    stop("length of 'keepX.bat' must be equal to ", ncomp.bat, ".")}

  if (any(keepX.bat > p)){
    stop("each component of 'keepX' must be lower than or equal to ", p, ".")}

  #########
  #-- initialisation of matrices --#
  if(is.null(rownames(X))){
    if(is.null(names(Y.bat))){
      rownames(X) = rownames(Y.bat.mat) = names(Y.bat) = paste('S', 1:n)
    }
    rownames(X) = rownames(Y.bat.mat) = names(Y.bat)
  }
  rownames(Y.bat.mat) = rownames(X)

  if(is.null(colnames(X))){
    colnames(X) = paste0('X', 1:p)
  }

  colnames(Y.bat.mat) = levels(Y.bat)

  weight <- rep(1,nrow(X))
  ##########################
  if(!is.null(Y.trt)){
    # Testing the input Y
    if (!is.null(dim(Y.trt))){
      stop("'Y.trt' should be a factor or a class vector.")}

    Y.trt = as.factor(Y.trt)
    Y.trt.mat = mixOmics::unmap(as.numeric(Y.trt))
    Y.trt.mat = as.matrix(Y.trt.mat)
    q.trt = ncol(Y.trt.mat)

    if(n != nrow(Y.trt.mat)){
      stop("unequal number of rows in 'X' and 'Y.trt'.")}

    if(is.null(ncomp.trt) || !is.numeric(ncomp.trt) || ncomp.trt <= 0){
      stop("invalid number of variates, 'ncomp.trt'.")}

    ncomp.trt = round(ncomp.trt)

    ###########
    if(ncomp.trt > p){
      warning("Reset maximum number of variates 'ncomp.trt' to ncol(X) = ", p, ".")
      ncomp.trt = p
    }

    if (length(keepX.trt) != ncomp.trt){
      stop("length of 'keepX.trt' must be equal to ", ncomp.trt, ".")}

    if (any(keepX.trt > p)){
      stop("each component of 'keepX' must be lower than or equal to ", p, ".")}

    colnames(Y.trt.mat) = levels(Y.trt)
    rownames(Y.trt.mat) = names(Y.trt) = rownames(X)

    ######################
    if(balance == F){
      #-- weight --#
      group <- data.frame(trt = Y.trt, bat = Y.bat)
      weight.num <- table(Y.trt, Y.bat)

      for(i in 1:nlevels(Y.trt)){
        for(j in 1:nlevels(Y.bat)){
          weight[group$trt == levels(Y.trt)[i] & group$bat == levels(Y.bat)[j]] = weight.num[i,j]
        }
      }

      weight <- 1/ weight
      weight <- weight/min(weight)
      weight <- sqrt(weight)

      X.wt <- weight*X
      X.scale <- scale(X.wt, center = TRUE, scale = TRUE)
      X.mean <- attributes(X.scale)$`scaled:center`
      X.sd <- attributes(X.scale)$`scaled:scale`

      Y.bat.wt <-  weight*Y.bat.mat
      Y.bat.scale = scale(Y.bat.wt, center = TRUE, scale = TRUE)

      Y.trt.wt <- weight*Y.trt.mat
      Y.trt.scale = scale(Y.trt.wt, center = TRUE, scale = TRUE)
    }else{
      X.scale <- scale(X, center = TRUE, scale = TRUE)
      X.mean <- attributes(X.scale)$`scaled:center`
      X.sd <- attributes(X.scale)$`scaled:scale`

      Y.bat.scale = scale(Y.bat.mat, center = TRUE, scale = TRUE)
      Y.trt.scale = scale(Y.trt.mat, center = TRUE, scale = TRUE)
    }

    ##########
    # stage 1: fit the treatment effect first --#
    plsda_trt <- PLSDA(X = X.scale, Y = Y.trt.scale, ncomp = ncomp.trt,
      keepX = keepX.trt, tol = tol, max.iter = max.iter)
    X.notrt <- plsda_trt$defl_data$X

  }else{
    X.scale <- scale(X, center = TRUE, scale = TRUE)
    X.mean <- attributes(X.scale)$`scaled:center`
    X.sd <- attributes(X.scale)$`scaled:scale`

    Y.bat.scale = scale(Y.bat.mat, center = TRUE, scale = TRUE)

    plsda_trt = NULL
    X.notrt = X.scale
  }

  ##########
  # stage 2: fit the batch effect then --#
  plsda_bat <- PLSDA(X = X.notrt, Y = Y.bat.scale, ncomp = ncomp.bat,
    keepX = keepX.bat, tol = tol, max.iter = max.iter)

  bat_loadings <- plsda_bat$loadings$a

  ##########
  #stage 3: deflation of batch from the original matrix --#
  X.temp <- X.scale
  for(h in 1:ncomp.bat){
    a.bat = bat_loadings[,h]
    t.bat = X.temp %*% a.bat
    X.temp <- mtx_deflation(X.temp, t.bat)
  }

  X.nobat <- X.temp

  X.nobat.final <- t(t(X.nobat)*X.sd + X.mean)
  X.nobat.final <- X.nobat.final/weight

  if(!is.null(Y.trt)){
    X.notrt.final <- t(t(X.notrt)*X.sd + X.mean)
    X.notrt.final <- X.notrt.final/weight
  }else{
    X.notrt.final = NULL
  }


  cl = match.call()
  cl[[1]] = as.name('PLSDA_batch')

  result = list(call = cl,
                X = X,
                X.nobatch = X.nobat.final,
                X.notrt = X.notrt.final,
                Y = list(trt = Y.trt, bat= Y.bat),
                latent_var.trt = plsda_trt$latent_var,
                latent_var.bat = plsda_bat$latent_var,
                loadings.trt = plsda_trt$loadings,
                loadings.bat = plsda_bat$loadings,
                tol = tol,
                max.iter = max.iter,
                iter.trt = plsda_trt$iters,
                iter.bat = plsda_bat$iters,
                explained_variance.trt = plsda_trt$exp_var,
                explained_variance.bat = plsda_bat$exp_var,
                weight = weight)

  if(near.zero.var == TRUE) result$nzv = nzv

  return(invisible(result))
}


