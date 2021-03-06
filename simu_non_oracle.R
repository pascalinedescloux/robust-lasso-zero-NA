#### Simulation with automatic tuning


## y is generated by y = X beta + epsilon
## (beta is s-sparse)
## one observes y and an incomplete version of X
## 
## THE COMPLETE X IS GENERATED AS FOLLOWS:
## - X has size n x p
## - step1: the rows of X are iid from N(0, Sigma), Sigma toeplitz matrix
## - step2: the columns of X are then mean centered and standardized to have 
##      unit standard deviation (before introduction of NAs)
## 
## INTRODUCING NAs IN X:
## We specify pi.tot: the overall probability that an entry is missing.
## Then: an entry with value x is declared missing with probability 
##          1 / (1 + exp(-a * abs(x) - b)),
##          where a >= 0 is a param (a = 0: MCAR, a > 0: MNAR),
##          and b is determined so that P(missing) ~= pi.tot (targeted total proportion of NA)
## 
## 
## ESTIMATORS:
## (1) Lasso-Zero
## (2) Robust Lasso-Zero with 1 parameter (fix lambda = 1)
## (3) Lasso
## (4) ABSLOPE
## 
## Methods (1) - (3) require preliminary imputation of the missing entries.
## Imputation & standardization is done in the following way:
## - for each column j: compute the observed mean mu_j
## - impute each NAs by the observed mean of the corresponding column
## - mean center each column
## - devide each column by (full) standard deviation

## -- load functions, data and pkgs
source("functions_non_oracle.R")
source("fun_ABSLOPE.R")
load("miss.proba.logpar.Rda")
library(lass0)
library(qut)
library(doParallel)
library(MASS)
registerDoParallel(cores = detectCores()) # parallel computing

## -- simulation parameters:
a <- 0 # 0(MCAR) or 5(MNAR)
s <- 3 # sparsity index of beta^0
pi.tot <- 0.05 # expected proportion of missing values in X
# finding parameter b:
a.ix <- which(dimnames(miss.proba.logpar)$a.value == a)
miss.proba <- miss.proba.logpar[a.ix, ]
b.ix <- which.min(abs(miss.proba - pi.tot))
b <- as.numeric(dimnames(miss.proba.logpar)$b.value[b.ix])
R <- 100 # number of replications
n <- 100 # sample size
p <- 200 # number of covariates
sigma <- 0.5 # noise level
rho <- 0 # correlation parameter
Sigma <- toeplitz(rho^(0:(p-1))) # covariance matrix of covariates
alpha <- 0.05 # level for QUT
MCrep <- 200 # number of replications for estimating QUT
GEVapprox <- TRUE # estimating QUT from a GEV fit
simu.parameters <- list(R = R, n = n, p = p, s = s,
                        rho = rho, sigma = sigma, a = a, b = b,
                        pi.tot = pi.tot, alpha = alpha, MCrep = MCrep)
parallel <- TRUE # whether to use parallel foreach loops inside the main loop

## -- function computing various criteria
comp.crit <- function(betahat, beta) {
    s <- sum(beta != 0)
    exactrec <- all((betahat == 0) == (beta == 0))
    TPP <- sum(beta != 0 & betahat != 0) / s
    FDP <- sum(betahat != 0 & beta == 0) / max(c(1, sum(betahat != 0)))
    exactsign <- all(sign(betahat) == sign(beta))
    TPPs <- (sum(betahat > 0 & beta > 0) + sum(betahat < 0 & beta < 0)) / s
    FDPs <- (sum(betahat != 0) - sum(beta > 0 & betahat > 0) - sum(beta < 0 & betahat < 0)) / max(c(1, sum(betahat != 0)))
    c(exactrec, TPP, FDP, exactsign, TPPs, FDPs)
}

## -- estimators:
est.names <- c( "Rlass0", "lass0",  "lasso", "ABslope")

## -- array for results:
res <- array(NA, dim = c(length(est.names), 6, R),
             dimnames = list(estimator = est.names, 
                             criterion = c("exactrec", "TPP", "FDP", "exactsign", "TPPs", "FDPs"),
                             rep = 1:R))

## -- simu:
for (r in 1:R) {
    print(r)
    set.seed(r)
    # random support and signs
    S <- sample(p, s)
    beta <- rep(0, p)
    beta[S] <- 1 * sign(rnorm(s))
    # generating complete X
    X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    # mean centering and standardization of X:
    X <- t(t(X) - colMeans(X))
    X <- t(t(X) / apply(X, 2, sd))
    # generating response vector:
    y <- X[, S, drop = FALSE] %*% beta[S] + rnorm(n, 0, sigma)
    y <- y - mean(y)
    # missing entries in X:
    U <- matrix(runif(n*p), n, p)
    isNA <- (U < 1 / (1 + exp(-a * abs(X) - b)))
    XNA <- X
    XNA[isNA] <- NA
    
    # impute & standardize XNA:
    rowsIx <- which(apply(is.na(XNA), 1, sum) > 0)
    XNAmeans <- apply(XNA, 2, mean, na.rm = TRUE)
    XNA <- t(t(XNA) - XNAmeans)
    Ximp <- XNA
    Ximp[is.na(Ximp)] <- 0 # at this point: matrix is imputed and mean-centered
    sds <- apply(Ximp, 2, sd)
    Ximp <- t(t(Ximp) / sds)
    
    # sub-identity matrix for robust lasso and lass0:
    subI <- diag(n)
    subI <- subI[, rowsIx, drop = FALSE]
    subI <- subI - 1/n # X is mean-centered -> so is subI
    
    # ratio lambda_beta / lambda_w used when robust lasso/lass0 are tuned with 1 parameter
    delta <- sqrt(n) # this way, all columns in [Ximp, G, delta*I] have the same observed std.
    delta.seq <- unique(c(seq(0.0001 * delta, delta, length = 20), 
                          seq(delta, 10 * delta, length = 20)))
    beta.names <- paste0("beta.", est.names)
    
    # estimation of sigma:
    sigma.est <- sigmarcv(y, base::cbind(Ximp, delta * subI), cv = TRUE)$sigmahat
    
    if("lass0" %in% est.names) {
        beta.lass0 <- lass0(Ximp, y, alpha = alpha, standardizeG = TRUE, 
                                parallel = parallel, intercept = TRUE, ols = FALSE, 
                                MCrep = MCrep, GEVapprox = GEVapprox)$coefficients
    }
    if ("Rlass0" %in% est.names) {
        beta.Rlass0.1par.QUT <- lass0(base::cbind(Ximp, delta * subI), y, alpha = alpha, 
                                                 standardizeG = TRUE, 
                                                 parallel = parallel, intercept = TRUE, ols = FALSE, 
                                                 MCrep = MCrep, GEVapprox = GEVapprox, var.subset = 1:p)$coefficients
        beta.Rlass0 <- beta.Rlass0.1par.QUT[1:p]
    }
    if ("lasso" %in% est.names) {
        cv.lasso.res <- cv.glmnet(Ximp, y, family = "gaussian", standardize = FALSE,
                                  intercept = TRUE, parallel = parallel)
        beta.lasso <- coef(cv.lasso.res, s = "lambda.min")[-1]
        
    }
    if("ABslope" %in% est.names) {
        ABslope.result <- ABSLOPE(XNA, y, lambda = qnorm(1-alpha*(1:p)/(2*p)), 
                         a = .01*n, b = .1*n)
        beta.ABslope <- ABslope.result$beta.new
    }
    
    estimators <- do.call(base::cbind, mget(beta.names))
    colnames(estimators) <- est.names
    res[, , r] <- t(apply(estimators, 2, comp.crit, beta = beta))
}

results <- list(res = res,
                simu.parameters = simu.parameters)

## -- save results:
assign(paste0("automatic.s", s, ".pi", pi.tot, ".a", a, ".rho", rho), results)
save(list=paste0("automatic.s", s, ".pi", pi.tot, ".a", a, ".rho", rho), 
     file=paste0("automatic.s", s, ".pi", pi.tot, ".a", a, ".rho", rho, ".Rda"))
