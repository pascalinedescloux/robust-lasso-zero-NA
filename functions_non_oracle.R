### Functions computing estimators in the non oracle (automatic) tuning case

## -- load pkgs
library(knockoff)
library(lass0)
library(glmnet)
library(abind)
library(doParallel)


## -- Robust Lasso-Zero (2 parameters: lambda & tau):
robustLass0.QUTGIC <- function(X, y, sigma, alpha, rowsIx = 1:nrow(X), delta.seq,
                               MCrep, GEVapprox = FALSE, parallel = TRUE,
                               intercept = TRUE, standardizeG, ...) {
    
    registerDoParallel(cores = detectCores())
    
    n <- nrow(X)
    p <- ncol(X)
    # sub-identity matrix:
    subI <- diag(n)
    subI <- subI[, rowsIx, drop = FALSE]
    # mean-centering:
    if (intercept) {
        meansX <- colMeans(X)
        X <- t(t(X) - meansX)
        meany <- mean(y)
        y <- y - meany
        subI <- subI - 1/n
    }
    taus.qut <- numeric(length(delta.seq))
    betas <- array(NA, dim = c(p, length(delta.seq)))
    GICs <- numeric(length(delta.seq))
    for(i in 1:length(delta.seq)) {
        print(paste0("for loop over delta.seq, par n° ", i, " out of ", length(delta.seq)))
        lass0.res <- lass0(base::cbind(X, delta.seq[i] * subI), y, alpha = alpha, 
                           standardizeX = FALSE, standardizeG = TRUE, MCrep = MCrep,
                           GEVapprox = GEVapprox, parallel = parallel, var.subset = 1:p,
                           ols = TRUE)
        betas[, i] <- lass0.res$coefficients[1:p]
        GICs[i] <- sum((y - base::cbind(X, delta.seq[i] * subI) %*% lass0.res$coefficients)^2) + 
            sigma^2 * log(log(n)) * log(p+n) * sum(lass0.res$coefficients != 0)
    }
    ix <- which.min(GICs)
    betahat <- betas[, ix]
    
    
    out <- NULL
    out$betahat <- betahat
    out$deltaGIC <- delta.seq[ix]
    out
}


## -- Lasso-Zero
lass0.byGIC<- function(X, y, sigma, standardizeG, intercept = TRUE, parallel = TRUE) {
    
    registerDoParallel(cores = detectCores())
    
    n <- nrow(X)
    p <- ncol(X)
    if (intercept) {
        meany <- mean(y)
        y <- y - meany
        meansX <- colMeans(X)
        X <- t(t(X) - meansX)
    }
    
    betamed <- lass0(X, y, tau = 0, standardizeX = FALSE, intercept = intercept,
                     standardizeG = standardizeG, parallel = parallel, ols = FALSE)$coefficients
    taus <- sort(abs(betamed), decreasing = TRUE)
    all.GIC <- sapply(1:length(taus), function(i) {
        tau <- taus[i]
        betahat <- betamed
        betahat[abs(betahat) <= tau] <- 0
        sum((y - X %*% betahat)^2) + sigma^2 * log(log(n)) * log(p) * sum(betahat != 0)
    })
    ix <- which.min(all.GIC)
    tau <- taus[ix]
    betahat <- betamed
    betahat[abs(betahat) <= tau] <- 0
    
    
    out <- list()
    out$betahat <- betahat
    if (intercept) {
        out$intercept <- meany - betahat %*% meansX
    } else {
        out$intercept <- NULL
    }
    return(out)
}


## -- Lasso:
lasso.byGIC <- function(X, y, sigma, intercept=FALSE, standardize=FALSE){
    n <- nrow(X)
    p <- ncol(X)
    if(intercept){
        meansX <- colMeans(X)
        X <- t(t(X) - meansX)
        meany <- mean(y)
        y <- y - meany
    }
    if(standardize){
        sdsX <- apply(X, 2, sd)
        X <- t(t(X) / sdsX)
    }else{sdsX <- rep(1, p)}
    
    res <- glmnet(X, y, standardize=FALSE, intercept=FALSE)
    lambda.seq <- res$lambda
    GIC <- sapply(1:length(lambda.seq), function(i){
        sum((y - X %*% res$beta[, i])^2) + sigma^2 * log(log(n)) * log(p) * sum(res$beta[, i] != 0)
    })
    ix <- which.min(GIC)
    betahat <- res$beta[, ix]
    betahat <- betahat / sdsX
    if(intercept) intercept <- meany - betahat %*% meansX
    
    out <- list()
    out$betahat <- betahat
    if(intercept) out$intercept <- intercept
    return(out)
}