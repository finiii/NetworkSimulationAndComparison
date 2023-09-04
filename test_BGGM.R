library(BGGM)

#hier kann ich mich dann an der vignette orientieren

# need these packages
library(BGGM)
library(ggplot2)
library(assortnet)
library(networktools)
library(MASS)


# group 1
Yg1 <- MASS::mvrnorm(n = 926,
                     mu = rep(0, 16),
                    Sigma = ptsd_cor3,
                    empirical = TRUE)

# group 2
Yg2 <- MASS::mvrnorm(n = 956,
                     mu = rep(0, 16),
                     Sigma = ptsd_cor4,
                     empirical = TRUE)



# test partial correlation matrix distance

f <- function(Yg1, Yg2){
  # nodes
  p <- ncol(Yg1)

  # index of off-diagonal
  indices <- upper.tri( diag(p))

  # fit models
  fit1 <-  BGGM::estimate(Yg1, analytic = TRUE)
  fit2 <-  BGGM::estimate(Yg2, analytic = TRUE)

  pcor1 <- BGGM::pcor_mat(fit1)
  pcor2 <- BGGM::pcor_mat(fit2)

  # CDM for partial correlations
  # note: numerator is the trace; denominator is the Frobenius norm
  1 - (sum(diag(pcor1 %*% pcor2)) / (norm(pcor1, type = "f") * norm(pcor2, type = "f")))
}

obs <- f(Yg1, Yg2)

# observed
obs

ppc <- BGGM::ggm_compare_ppc(Yg1, Yg2,
                       FUN = f,
                       custom_obs = obs,
                       iter = 1000)

ppc




# Hamming distance

f <- function(Yg1, Yg2){
  # nodes
  p <- ncol(Yg1)

  # index of off-diagonal
  indices <- upper.tri( diag(p))

  # fit models
  fit1 <-  BGGM::estimate(Yg1, analytic = TRUE)
  fit2 <-  BGGM::estimate(Yg2, analytic = TRUE)

  # select graphs
  sel1 <- BGGM::select(fit1)
  sel2 <- BGGM::select(fit2)

  # hamming distance
  sum((sel1$adj[indices] - sel2$adj[indices]) ^ 2)
}

obs <- f(Yg1, Yg2)

# observed
obs

ppc <- BGGM::ggm_compare_ppc(Yg1, Yg2,
                       FUN = f,
                       custom_obs = obs,
                       iter = 1000)

ppc


pdf(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/BGGM.pdf")
plot(ppc)
dev.off()

