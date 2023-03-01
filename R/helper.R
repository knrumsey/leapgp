

leap_pred_adapt <- function(object, newdata, rho,
                            scale, n, start, M_max, ...){
  X <- object$X
  y <- object$y
  hubs <- object$hubs
  Xnew <- newdata

  N <- length(y)
  if(is.na(n)){
    n <- ceiling(sqrt(N))
  }
  if(is.na(start)){
    start <- floor(max(6, 0.1*n))
  }

  #Find nearest hub (if possible)
  if(length(hubs) > 0){
    #Convert coord ID to coord locs if needed
    if(is.null(hubs[[1]]$epsilon)){
      for(h in 1:length(hubs)){
        hubs[[h]]$epsilon <- sqrt(log(1/rho)/hubs[[h]]$kappa)
        hubs[[h]]$coord <- X[hubs[[h]]$coord_id, ]
      }
    }

    #Matrix of hub locations
    H <- matrix(unlist(lapply(hubs, function(h) h$coord)), nrow=length(hubs), byrow=TRUE)
    #Find nearest hub using kd-tree
    #browser()
    NH <- RANN::nn2(H, matrix(Xnew, nrow=1), k=1,
                    searchtype = ifelse(nrow(H) < 100, "standard", "priority"), eps = 0.025*floor(nrow(H/100)))
    curr_hub <- hubs[[NH$nn.idx]]
    #Is hub close enough to make a prediction?
    #alternatively, are we already at the maximum number of hubs
    if(NH$nn.dists < curr_hub$epsilon | nrow(H) >= M_max){
      #Make prediction and return hubs unmodified
      kvec <- rep(NA, n)
      for(i in 1:n){
        ii <- curr_hub$neigh[i]
        kvec[i] <- exp(-curr_hub$kappa*sum((X[ii,]-Xnew)^2))
      }
      #Prepare object for return
      out <- list(X=X, y=y)
      out$mean <- matrix(kvec, nrow=1)%*%curr_hub$a
      #out$pred_hub <- NH$nn.idx
      if(scale){
        out$sd <- sqrt(curr_hub$phi/n*(1 - matrix(kvec, nrow=1)%*%curr_hub$Kinv%*%matrix(kvec, ncol=1)))
      }
      out$hubs <- hubs
      return(out)
    }
  }
  #Make a new hub (and use it for prediction)
  h <- laGP::laGP(Xref=matrix(Xnew, ncol=ncol(X)), X=X, Z=y,
                  start=start, end=n, ...)
  #Prepare object for return
  out <- list(X=X, y=y,
              mean=h$mean)
  if(scale) out$sd <- sqrt(h$s2)

  #Make a new hub list
  new_hub <- list(coord=Xnew, neigh=h$Xi, kappa=1/h$mle$d,
                  nugget=1/h$mle$g)
  new_hub$epsilon <- sqrt(-log(rho)/new_hub$kappa)
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
  #Add new hub to hubs
  i <- length(hubs) + 1
  hubs[[i]] <- new_hub
  out$hubs <- hubs
  return(out)
}


leap_pred_no_adapt <- function(object, newdata, scale=FALSE, ...){
  X <- object$X
  y <- object$y
  hubs <- object$hubs
  Xnew <- newdata

  #Get parameters
  N <- length(y); n <- length(hubs[[1]]$neigh)
  if(ncol(X) != length(Xnew)){
    stop('wrong dimensions for Xnew')
  }
  #Assume single prediction requested
  H <- X[unlist(lapply(hubs, function(h) h$coord_id)),]
  #Find nearest hub using kd-tree
  NH <- RANN::nn2(H, matrix(Xnew, nrow=1), k=1)
  curr_hub <- hubs[[NH$nn.idx]]
  #Make prediction and return hubs unmodified
  kvec <- rep(NA, n)
  for(i in 1:n){
    ii <- curr_hub$neigh[i]
    kvec[i] <- exp(-curr_hub$kappa*sum((X[ii,]-Xnew)^2))
  }
  pred <- as.numeric(matrix(kvec, nrow=1)%*%curr_hub$a)
  out <- list(X=X, y=y, hubs=hubs, mean=pred)
  if(scale){
    scale <- curr_hub$phi/n*(1 - matrix(kvec, nrow=1)%*%curr_hub$Kinv%*%matrix(kvec, ncol=1))
    out$sd <- sqrt(scale)
  }
  return(out)
}

leap_pred_no_adapt_convert <- function(object, newdata, rho, scale=FALSE, ...){
  X <- object$X
  y <- object$y
  hubs <- object$hubs
  Xnew <- newdata

  #Convert coord ID to coord locs
  for(h in 1:length(hubs)){
    hubs[[h]]$epsilon <- sqrt(log(1/rho)/hubs[[h]]$kappa)
    hubs[[h]]$coord <- X[hubs[[h]]$coord_id, ]
  }
  object$hubs <- hubs

  #Get parameters
  X <- object$X; y <- object$y; hubs <- object$hubs
  N <- length(y); n <- length(hubs[[1]]$neigh)
  if(ncol(X) != length(Xnew)){
    stop('wrong dimensions for Xnew')
  }
  #Assume single prediction requested
  H <- X[unlist(lapply(hubs, function(h) h$coord_id)),]
  #Find nearest hub using kd-tree
  NH <- RANN::nn2(H, matrix(Xnew, nrow=1), k=1)
  curr_hub <- hubs[[NH$nn.idx]]
  #Make prediction and return hubs unmodified
  kvec <- rep(NA, n)
  for(i in 1:n){
    ii <- curr_hub$neigh[i]
    kvec[i] <- exp(-curr_hub$kappa*sum((X[ii,]-Xnew)^2))
  }
  pred <- as.numeric(matrix(kvec, nrow=1)%*%curr_hub$a)
  out <- list(X=X, y=y, hubs=hubs, mean=pred)
  if(scale){
    scale <- curr_hub$phi/n*(1 - matrix(kvec, nrow=1)%*%curr_hub$Kinv%*%matrix(kvec, ncol=1))
    out$sd <- sqrt(scale)
  }
  return(out)
}
