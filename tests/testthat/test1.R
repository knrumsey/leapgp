test_that("Simple bivariate function", {
  cat('A quick emulation problem')

  f <- function(x){
    1.3356*(1.5*(1-x[1]) + exp(2*x[1] - 1)*sin(3*pi*(x[1] - 0.6)^2) + exp(3*(x[2]-0.5))*sin(4*pi*(x[2] - 0.9)^2))
  }

  # Generate data
  N <- 500
  X <- matrix(runif(N * 2), nrow=N, ncol=2)
  y <- apply(X, 1, f)

  Xtest <- matrix(runif(50 * 2), ncol=2)
  ytest <- apply(Xtest, 1, f)

  # Train leapGP model
  mod <- leapGP(X, y, M0=10)

  pred <- rep(NA, 50)
  for(i in 1:50){
    mod <- predict_leapGP(mod, Xtest[i,], rho=0.95)
    pred[i] <- mod$mean
  }

  # Test prediction
  d1 <- cor(ytest, pred)
  expect_that(d1, is_more_than(0.9))
})
