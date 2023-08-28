

simulate_multivariate_poisson <- function(n, mu, Sigma) {
    Cor <- cov2cor(Sigma)
    SDs <- sqrt(diag(Sigma))
    d   <- length(SDs)
    normal  <- rmvnorm(n, rep(0, d), Cor)
    unif   <- pnorm(normal)
    data <- t(qpois(t(unif), mu, ...))
    return(data)
}

n = 100
mu  = 
