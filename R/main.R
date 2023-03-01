#' Localized Ensemble of Approximate Gaussian Processes
#'
#' This function is a modification of the laGP framework of Gramacy and Apley
#' designed for cases where parallel predictions are not possible (e.g. MCMC).
#' The leapGP offers a quadratic training algorithm which leads to fast predictions.
#'
#' @param X a matrix of training locations (1 row for each training instance)
#' @param y a vector of training responses (length(y) == nrow(X))
#' @param M0 the number of prediction hubs desired. Defaults to ceiling(sqrt(length(Y))).
#' @param rho (optional). The parameter controlling time-accuracy tradeoff. Can alternatively be specified during prediction.
#' @param scale logical. Do we want the scale parameter to be returned for predictions? If TRUE,
#'            the matrix K^{-1} will be stored for each hub.
#' @param n local neighborhood size (for laGP)
#' @param start number of starting points for neighborhood (between 6 and n inclusive)
#' @param verbose  logical. Should status be printed? Deault is FALSE
#' @param justdoit logical. Force leapGP to run using specified parameters (may take a long time and/or cause R to crash).
#' @param ... optional arguments to be passed to laGP()
#' @return a univariate prediction and an updated list of hubs. Also returns scale parameter if scale=TRUE
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


#' Predict metthod for leapGP class
#'
#' This function is a modification of the LA-GP framework of Gramacy and Apley
#' designed for cases where parallel predictions are not possible (i.e. MCMC).
#' The slapGP framework offers users a time-accuracy tradeoff based on the rho parameter.
#'
#' @param object An object of class `leapGP`
#' @param newdata New data
#' @param rho parameter controlling time-accuracy tradeoff (default = 0.95)
#' @param scale logical. Do we want the scale parameter to be returned for predictions? If TRUE,
#'            the matrix K^{-1} will be stored for each hub.
#' @param n local neighborhood size
#' @param start number of starting points for neighborhood (between 6 and n inclusive)
#' @param M_max the maximum number of hubs allowed (used to upper bound the run time)
#' @param ... optional arguments to be passed to laGP()
#' @return A list containing values `mean`, `hubs` `X` and `y`. If `scale=TRUE` the list also contains field `sd`.
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
  return(out)
}
