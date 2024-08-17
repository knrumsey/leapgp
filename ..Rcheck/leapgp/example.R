library(leapgp)

# Generate data
f <- function(x) 1.3356*(1.5*(1-x[1]) + exp(2*x[1] - 1)*sin(3*pi*(x[1] - 0.6)^2) + exp(3*(x[2]-0.5))*sin(4*pi*(x[2] - 0.9)^2))
X <- matrix(runif(2000), ncol=2)
y <- apply(X, 1, f)

# Generate data for prediction
Xtest <- matrix(runif(2000), ncol=2)
ytest <- apply(Xtest, 1, f)

# ============================
#    CASE ONE: rho = 0
# ============================
# Train initial model
mod <- leapGP(X, y, M0 = 30)
# Make sequential predictions
pred <- rep(NA, 1000)
for(i in 1:1000){
  mod <- predict_leapGP(mod, matrix(Xtest[i,], nrow=1), rho=0)
  pred[i] <- mod$mean
}


# ============================
#    CASE TWO: M0 = 0
# ============================
# Initalize empty model
mod <- leapGP(X, y, M0 = 0)
# Make sequential predictions
pred <- rep(NA, 1000)
for(i in 1:1000){
  mod <- predict_leapGP(mod, matrix(Xtest[i,], nrow=1), rho=0.95)
  pred[i] <- mod$mean
}


# ============================
#    CASE THREE: Full leapGP
# ============================
# Train leapGP model
mod <- leapGP(X, y, M0 = 30)
# Make sequential predictions
pred <- rep(NA, 1000)
for(i in 1:1000){
  mod <- predict_leapGP(mod, matrix(Xtest[i,], nrow=1), rho=0.99)
  pred[i] <- mod$mean
}

# TEST CHUNK PREDICTION (NOT THE INTENDED USECASE)
mod <- leapGP(X, y, M0=30)
mod <- predict_leapGP(mod, Xtest, rho=0.95)
mod$mean

# TEST SCALE FEATURE
mod <- leapGP(X, y, M0=30, scale=TRUE)
mod <- predict_leapGP(mod, Xtest, rho=0.95, scale=TRUE)
mod$mean
mod$sd



