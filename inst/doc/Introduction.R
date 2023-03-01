## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(leapgp)
library(RColorBrewer)
library(lhs)
library(laGP)
library(tictoc)

## -----------------------------------------------------------------------------
f1 <- function(x){
  x1 <- x[1]
  x2 <- x[2]

  term1 <- (5/2)*x1 - (35/2)*x2
  term2 <- (5/2)*x1*x2 + 19*x2^2
  term3 <- -(15/2)*x1^3 - (5/2)*x1*x2^2
  term4 <- -(11/2)*x2^4 + (x1^3)*(x2^2)

  y <- (9 + term1 + term2 + term3 + term4)*11/20
  return(y)
}

f2 <- function(x){
  5*exp(-((8*x[1]-2)^2 + (8*x[2] - 2)^2))*(8*x[1] - 2)
}

twin_galaxies <- function(x) f1(x) + f2(x)

xx <- seq(0, 1, length.out=51)
zz <- expand.grid(xx, xx)
yy <- matrix(apply(zz, 1, twin_galaxies), nrow=51)
image(yy, 
      col=colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
      xlab="x1", ylab="x2")

## -----------------------------------------------------------------------------
X <- lhs::maximinLHS(500, 2)
y <- apply(X, 1, twin_galaxies)

Xtest <- matrix(runif(2*50), ncol=2)
ytest <- apply(Xtest, 1, twin_galaxies)

## -----------------------------------------------------------------------------
tic()
mod <- leapGP(X, y, M0=30)
toc()

tic()
pred1 <- rep(NA, 50)
for(i in 1:50){
  mod <- predict_leapGP(mod, Xtest[i,], rho=0.95)
  pred1[i] <- mod$mean
}
toc()
plot(pred1, ytest, main=paste("RMSE = ", round(sqrt(mean((pred1-ytest)^2)), 4)))


## -----------------------------------------------------------------------------
tic()
mod <- leapGP(X, y, M0=0)
toc()

## -----------------------------------------------------------------------------
tic()
pred2 <- rep(NA, 50)
for(i in 1:50){
  mod <- predict_leapGP(mod, Xtest[i,], rho=0.95)
  pred2[i] <- mod$mean
}
toc()
plot(pred2, ytest, main=paste("RMSE = ", round(sqrt(mean((pred2-ytest)^2)), 4)))

## -----------------------------------------------------------------------------
tic()
mod <- leapGP(X, y, M0=30)
toc()

tic()
pred3 <- rep(NA, 50)
for(i in 1:50){
  mod <- predict_leapGP(mod, Xtest[i,], rho=0)
  pred3[i] <- mod$mean
}
toc()
plot(pred3, ytest, main=paste("RMSE = ", round(sqrt(mean((pred3-ytest)^2)), 4)))

## ---- eval=FALSE--------------------------------------------------------------
#  
#  mod <- leapGP(X, y, M0=30, scale=TRUE)
#  
#  pred <-  sd <- rep(NA, 50)
#  for(i in 1:50){
#    mod <- predict_leapGP(mod, Xtest[i,], rho=0.95, scale=TRUE)
#    pred[i] <- mod$mean
#    sd[i] <- mod$sd
#  }

## ---- eval=FALSE--------------------------------------------------------------
#  mod  <- leapGP(X, y, M0=30, scale=TRUE)
#  mod  <- predict_leapGP(mod, Xtest, rho=0.95, scale=TRUE)
#  pred <- mod$mean
#  sd   <- mod$sd

