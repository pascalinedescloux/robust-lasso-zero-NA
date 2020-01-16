### Missing data mechanism based on logistic function:
### 
### P(X_ij is missing | X_ij = x) = 1 / (1 + e^{-a |x| - b}),
### 
### with a >= 0, b real.
### 
### We considered a = 0 (data is MCAR) and a = 5 (data is MNAR)
### 
### This scripts builds a table containing estimations of P(X is missing) with 
### a = 0 and a = 5 respectively, and with a large grid of b values; with 
### X ~ N(0,1).


## -- a and b values
a.values <- c(0,5)
b.values <- seq(-15, 0, length=10000)

## -- array for results
miss.proba <- array(dim = c(length(a.values), length(b.values)))
dimnames(miss.proba) <- list(a.value = a.values, b.value = b.values)

## -- looping over all a and b
set.seed(2019)
R <- 1.e5 # drawn sample size for fixed (a,b)
for(i in 1:length(a.values)) {
    a <- a.values[i]
    print(i)
    for(j in 1:length(b.values)) {
        print(j)
        b <- b.values[j]
        x <- rnorm(R)
        u <- runif(R)
        miss.proba[i, j] <- sum(as.numeric(u < 1 / (1 + exp(-a*abs(x) - b)))) / R
    }
}

## -- save data
miss.proba.logpar <- miss.proba
save(miss.proba.logpar, file="miss.proba.logpar.Rda")

        