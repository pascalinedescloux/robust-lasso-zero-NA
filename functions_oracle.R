### Functions computing estimators in the s-oracle tuning case

## -- load pkgs
library(lass0)
library(lars)
library(glmnet)
library(abind)
library(doParallel)

## -- Robust Lasso-Zero (2 parameters: lambda & tau)
robustLass0.Soracle <- function(X, y, rowsIx = 1:nrow(X), delta.seq, S, 
                                intercept = TRUE, standardizeG, deltaref = NULL, ...) {
    
    # Reparametrization: delta = sqrt(n) / lambda
    # deltaref: if non NULL, then among all delta values (in delta.seq) giving the highest TPP,
    # one chooses the one as close to deltaref as possible
    
    registerDoParallel(cores = detectCores())
    
    n <- nrow(X)
    s <- length(S)
    # tpp function
    comp.tpp <- function(beta, S) sum(which(beta != 0) %in% S) / length(S)
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
    all.betas <- array(NA, dim = c(ncol(X), length(delta.seq)))
    for (i in 1:length(delta.seq)) {
        message(paste0("robustlass0 delta value ", i, " out of ", length(delta.seq)))
        delta <- delta.seq[i]
        betahat <- lass0(base::cbind(X, delta * subI), y, tau = 0, 
                         standardizeX = FALSE, intercept = intercept,
                         standardizeG = standardizeG, parallel = TRUE, ols = FALSE)$coefficients
        betahat <- betahat[1:ncol(X)] # remove w coefficients
        thresh <- sort(abs(betahat), decreasing = TRUE)[s]
        betahat[abs(betahat) < thresh] <- 0
        all.betas[, i] <- betahat
    }
    
    
    all.tpp <- apply(all.betas, 2, comp.tpp, S = S)
    all.beta.ix <- which(all.tpp == max(all.tpp))
    if (!is.null(deltaref)) {
        deltadiff <- abs(delta.seq - deltaref)[all.beta.ix]
        beta.ix <- all.beta.ix[which.min(deltadiff)]
    } else {
        beta.ix <- which.max(all.tpp)
    }
    
    delta <- delta.seq[beta.ix]
    betahat <- all.betas[, beta.ix]
    
    out <- NULL
    out$betahat <- betahat
    out$delta <- delta
    
    out
}


## -- Lasso-Zero:
lass0.soracle<- function(X, y, s, standardizeG, intercept = TRUE, parallel = TRUE) {
    
    registerDoParallel(cores = detectCores())
    
    n <- nrow(X)
    if (intercept) {
        meany <- mean(y)
        y <- y - meany
        meansX <- colMeans(X)
        X <- t(t(X) - meansX)
    }
    
    betahat <- lass0(X, y, tau = 0, standardizeX = FALSE, intercept = intercept,
                     standardizeG = standardizeG, parallel = parallel, ols = FALSE)$coefficients
    thresh <- sort(abs(betahat), decreasing = TRUE)[s]
    betahat[abs(betahat) < thresh] <- 0
    
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
lasso.soracle <- function(X, y, s, intercept = TRUE, standardize = FALSE){
    
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
    
    lars.res <- lars(X, y, normalize = FALSE, intercept = FALSE)
    Msize <- apply(lars.res$beta!=0, 1, sum)
    betahat <- lars.res$beta[min(which(Msize >= s)), ]
    betahat <- betahat / sdsX
    
    out <- list()
    out$betahat <- betahat
    if (intercept) {
        out$intercept <- meany - betahat %*% meansX
    } else {
        out$intercept <- NULL
    }
    return(out)
}
