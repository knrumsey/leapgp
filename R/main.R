#' Localized Ensemble of Approximate Gaussian Processes
#'
#' Function to train or initialize a leapGP model, as described in Rumsey et al. (2023).
#'
#' @param X a matrix of training locations (1 row for each training instance)
#' @param y a vector of training responses (\code{length(y)} should equal \code{nrow(X)})
#' @param M0 the number of prediction hubs desired. Defaults to \code{ceiling(sqrt(length(Y)))}.
#' @param rho (optional). The parameter controlling time-accuracy tradeoff. Can also be specified during prediction.
#' @param scale logical. Do we want the scale parameter to be returned for predictions? If TRUE,
#'            the matrix \eqn{K^{-1}} will be stored for each hub.
#' @param n local neighborhood size (for laGP)
#' @param start number of starting points for neighborhood (between 6 and n inclusive)
#' @param verbose  logical. Should status be printed? Deault is FALSE
#' @param justdoit logical. Force leapGP to run using specified parameters (may take a long time and/or cause R to crash).
#' @param ... optional arguments to be passed to \code{laGP()}
#' @return an object of class \code{leapGP} with fields \code{X}, \code{y}, and \code{hubs}.  Also returns scale parameter if \code{scale=TRUE}
#' @details The leapGP is extends the laGP framework of Gramacy & Apley (2015). The methods are equivalent for \code{rho=1},
#'          but leapGP trades memory for speed when \code{rho < 1}. The method is described in Rumsey et al. (2023) where they demonstrate
#'          that leapGP is faster than laGP for sequential predictions and is also generally more accurate for some settings of \code{rho}.
#' @references
#' Gramacy, R. B., & Apley, D. W. (2015). Local Gaussian process approximation for large computer experiments. Journal of Computational and Graphical Statistics, 24(2), 561-578.
#'
#' Rumsey, K. N., Huerta, G., & Derek Tucker, J. (2023). A localized ensemble of approximate Gaussian processes for fast sequential emulation. Stat, 12(1), e576.
#' @examples
#' # Generate data
#' f <- function(x){
#'    1.3356*(1.5*(1-x[1]) + exp(2*x[1] - 1)*sin(3*pi*(x[1] - 0.6)^2) +
#'    exp(3*(x[2]-0.5))*sin(4*pi*(x[2] - 0.9)^2))
#' }
#' X <- matrix(runif(200), ncol=2)
#' y <- apply(X, 1, f)
#'
#' # Generate data for prediction
#' Xtest <- matrix(runif(200), ncol=2)
#' ytest <- apply(Xtest, 1, f)
#'
#' # Train initial model
#' mod <- leapGP(X, y, M0 = 30)
#' # Make sequential predictions
#' pred <- rep(NA, 100)
#' for(i in 1:100){
#'   mod <- predict_leapGP(mod, matrix(Xtest[i,], nrow=1), rho=0.9)
#'   pred[i] <- mod$mean
#' }
#' @export
leapGP <- function(X, y,
                   M0 = ceiling(sqrt(length(y))),
                   rho = NA,
                   scale=FALSE,
                   n = ceiling(sqrt(length(y))),
                   start = NA,
                   verbose = FALSE,
                   justdoit = FALSE, ...){
  if(M0 == 0){
    out <- list(X=X, y=y, hubs=list())
    class(out) <- "leapGP"
    return(out)
  }

  cnt <- 1
  #Get parameters
  N <- length(y)
  if(is.na(n)){
    n <- ceiling(sqrt(N))
  }
  if(n > 300 & !justdoit){
    stop('n is very large - consider choosing a smaller value. Set justdoit=TRUE to override.')
  }

  if(is.na(M0)){
    M0 <- ceiling(sqrt(N))
  }
  if(is.na(start)){
    start <- floor(max(6, 0.1*n))
  }

  if(M0 < N){
    if(N <= 5000 | justdoit){
      #Use PAM to find medoids
      coord_ind <- cluster::pam(X, M0, pamonce=5)$id.med
    }else{
      warning('N is very large, a subset of 5000 points was used for inital hub placement. Set justdoit=TRUE to override.')
      rand_ind <- sample(N,5000)
      coord_ind <- rand_ind[cluster::pam(X[rand_ind,], M0, pamonce=5)$id.med]
    }
  }else{
    coord_ind <- 1:N
  }

  #Start building hubs
  hubs <- list()
  for(hh in 1:M0){
    h <- laGP::laGP(Xref=matrix(X[coord_ind[hh], ], ncol=ncol(X)),
                    X=X, Z=y,
                    start=start, end=n, ...)
    new_hub <- list(coord_id=coord_ind[hh], neigh=h$Xi, kappa=1/h$mle$d,
                    nugget=1/h$mle$g)
    #Calculate psi (called a in code)
    K <- diag(rep(1, n))
    for(i in 2:n){
      for(j in 1:(n-1)){
        ii <- new_hub$neigh[i]
        jj <- new_hub$neigh[j]
        K[i,j] <- K[j, i] <- exp(-new_hub$kappa*sum((X[ii,]-X[jj,])^2))
      }
    }
    Kinv <- solve(K+diag(rep(1e-9, n)))
    new_hub$a <- Kinv%*%matrix(y[new_hub$neigh], nrow=n)
    if(scale){
      new_hub$Kinv <- Kinv
      new_hub$phi <- matrix(y[new_hub$neigh], ncol=n)%*%(new_hub$a)
    }
    hubs[[hh]] <- new_hub
    if(verbose){
      if((hh %% ((M0 - M0%%10)/10)) == 0){
        cat(cnt*10, "%,\t", as.character(Sys.time()), "\n", sep="")
        cnt <- cnt + 1
      }
    }
  }

  out <- list(hubs=hubs, X=X, y=y)
  class(out) <- "leapGP"
  return(out)
}


#' Predict Method for leapGP
#'
#' Predict method for an object of class leapGP.
#' Returns a (possibly modified) leapGP object as well as a prediction (with uncertainty, if requested).
#'
#' @param object An object of class \code{leapGP}
#' @param newdata New data
#' @param rho parameter controlling time-accuracy tradeoff (default is \code{rho=0.95})
#' @param scale logical. Do we want the scale parameter to be returned for predictions? If TRUE,
#'            the matrix \eqn{K^{-1}} will be stored for each hub.
#' @param n local neighborhood size
#' @param start number of starting points for neighborhood (between 6 and n inclusive)
#' @param M_max the maximum number of hubs allowed (used to upper bound the run time)
#' @param ... optional arguments to be passed to \code{laGP()}
#' @return A list containing values \code{mean}, \code{hubs} \code{X} and \code{y}. If \code{scale=TRUE} the list also contains field \code{sd}.
#' @details The leapGP is extends the laGP framework of Gramacy & Apley (2015). The methods are equivalent for \code{rho=1},
#'          but leapGP trades memory for speed when \code{rho < 1}. The method is described in Rumsey et al. (2023) where they demonstrate
#'          that leapGP is faster than laGP for sequential predictions and is also generally more accurate for some settings of \code{rho}.
#' @references
#' Gramacy, R. B., & Apley, D. W. (2015). Local Gaussian process approximation for large computer experiments. Journal of Computational and Graphical Statistics, 24(2), 561-578.
#'
#' Rumsey, K. N., Huerta, G., & Derek Tucker, J. (2023). A localized ensemble of approximate Gaussian processes for fast sequential emulation. Stat, 12(1), e576.
#' @examples
#' # Generate data
#' f <- function(x){
#'    1.3356*(1.5*(1-x[1]) + exp(2*x[1] - 1)*sin(3*pi*(x[1] - 0.6)^2) +
#'    exp(3*(x[2]-0.5))*sin(4*pi*(x[2] - 0.9)^2))
#' }
#' X <- matrix(runif(200), ncol=2)
#' y <- apply(X, 1, f)
#'
#' # Generate data for prediction
#' Xtest <- matrix(runif(200), ncol=2)
#' ytest <- apply(Xtest, 1, f)
#'
#' # Train initial model
#' mod <- leapGP(X, y, M0 = 30)
#' # Make sequential predictions
#' pred <- rep(NA, 100)
#' for(i in 1:100){
#'   mod <- predict_leapGP(mod, matrix(Xtest[i,], nrow=1), rho=0.9)
#'   pred[i] <- mod$mean
#' }
#' @export
predict_leapGP <- function(object, newdata, rho=0.95,
                           scale=FALSE, n=ceiling(sqrt(length(y))), start=NA, M_max=Inf, ...){
  X <- object$X
  y <- object$y
  if(is.na(start)){
    start <- floor(max(6, 0.1*n))
  }
  if(is.null(dim(newdata))){
    newdata <- matrix(newdata, nrow=1)
  }
  adapt_hubs <- TRUE
  if(is.na(rho) | rho == 0){
    adapt_hubs <- FALSE
  }
  if(length(object$hubs) >= M_max){
    adapt_hubs <- FALSE
  }

  # If newdata is a matrix,
  # loop over the rows
  if(nrow(newdata) > 1){
    mean <- sd <- rep(NA, nrow(newdata))
    hubs <- object$hubs
    for(i in 1:nrow(newdata)){
      emulator <- predict_leapGP(object, newdata[i,], rho, scale, n, start, M_max, ...)
      mean[i] <- emulator$mean
      if(scale) sd[i] <- emulator$sd
      hubs <- emulator$hubs
    }
    out <- list(hubs=hubs, X=X, y=y, mean=mean)
    if(scale) out$sd <- sd
    class(out) <- "leapGP"
    return(out)
  }

  if(adapt_hubs){
    out <- leap_pred_adapt(object, newdata, rho, scale, n, start, M_max, ...)
  }else{
    if(is.null(object$hubs[[1]]$epsilon)){
      out <- leap_pred_no_adapt_convert(object, newdata, scale, ...)
    }else{
      out <- leap_pred_no_adapt(object, newdata, scale, ...)
    }
  }
  class(out) <- "leapGP"
  return(out)
}
