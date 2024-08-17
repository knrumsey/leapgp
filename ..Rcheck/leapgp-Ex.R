pkgname <- "leapgp"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "leapgp-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('leapgp')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("leapGP")
### * leapGP

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: leapGP
### Title: Localized Ensemble of Approximate Gaussian Processes
### Aliases: leapGP

### ** Examples

# Generate data
f <- function(x){
   1.3356*(1.5*(1-x[1]) + exp(2*x[1] - 1)*sin(3*pi*(x[1] - 0.6)^2) +
   exp(3*(x[2]-0.5))*sin(4*pi*(x[2] - 0.9)^2))
}
X <- matrix(runif(200), ncol=2)
y <- apply(X, 1, f)

# Generate data for prediction
Xtest <- matrix(runif(200), ncol=2)
ytest <- apply(Xtest, 1, f)

# Train initial model
mod <- leapGP(X, y, M0 = 30)
# Make sequential predictions
pred <- rep(NA, 100)
for(i in 1:100){
  mod <- predict_leapGP(mod, matrix(Xtest[i,], nrow=1), rho=0.9)
  pred[i] <- mod$mean
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("leapGP", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict_leapGP")
### * predict_leapGP

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict_leapGP
### Title: Predict Method for leapGP
### Aliases: predict_leapGP

### ** Examples

# Generate data
f <- function(x){
   1.3356*(1.5*(1-x[1]) + exp(2*x[1] - 1)*sin(3*pi*(x[1] - 0.6)^2) +
   exp(3*(x[2]-0.5))*sin(4*pi*(x[2] - 0.9)^2))
}
X <- matrix(runif(200), ncol=2)
y <- apply(X, 1, f)

# Generate data for prediction
Xtest <- matrix(runif(200), ncol=2)
ytest <- apply(Xtest, 1, f)

# Train initial model
mod <- leapGP(X, y, M0 = 30)
# Make sequential predictions
pred <- rep(NA, 100)
for(i in 1:100){
  mod <- predict_leapGP(mod, matrix(Xtest[i,], nrow=1), rho=0.9)
  pred[i] <- mod$mean
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict_leapGP", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
