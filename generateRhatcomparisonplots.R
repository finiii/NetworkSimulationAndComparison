#give the function needed to calculate the rhat values
#I took this function from the coda package and changed it a bit
"gelman.preplot" <-
  function (x, max.bins = 5000, confidence = 0.95, bin.width = 1, autoburnin = TRUE) 
{

  nbin <- min(floor((niter(x) - 50)/thin(x)),max.bins)
  binw <- min(floor((niter(x) - 50)/nbin), bin.width)
  last.iter <- c(seq(from = start(x) + 50 * thin(x), by = binw * 
                     thin(x), length = nbin), end(x))
  shrink <- array(dim = c(nbin + 1, nvar(x), 2))
  dimnames(shrink) <- list(last.iter, varnames(x),
                           c("median", paste(50 * (confidence + 1), "%",
                                             sep = ""))
                           )
  #the first entries correponding to the burn in should be zero! autobutnin = TRUE cuts half of the samples
  for (i in 1:(nbin)) {

    shrink[i, , ] <- gelman.diag(window(x, end = last.iter[i]),
                                 multivariate = TRUE,
                                 autoburnin = TRUE)$psrf
  }
  return(list(shrink = shrink, last.iter = last.iter))
}





load("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/pPGM_simulation_dim10_lower_bound_-0.1_upper_bound_0_n_samples_6_iterations_5000.RData")
result10 <- result
load("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/pPGM_simulation_dim100_lower_bound_-0.1_upper_bound_0_n_samples_6_iterations_20000.RData")
result100 <- result

data10 <- result10$X_1_simulations
data100 <- result100$X_1_simulations

###prepare for R hat 
mcmc.list.10 <- mcmc.list(as.mcmc(data10[1,]), as.mcmc(data10[2,]), as.mcmc(data10[3,]) ,as.mcmc(data10[4,]), as.mcmc(data10[5,]), as.mcmc(data10[6,]))
mcmc.list.100 <-  mcmc.list(as.mcmc(data100[1,]), as.mcmc(data100[2,]), as.mcmc(data100[3,]) ,as.mcmc(data100[4,]), as.mcmc(data100[5,]), as.mcmc(data100[6,]))


R_hat_10 <- gelman.preplot(mcmc.list.10, max.bins = 2000)$shrink[,1,1]
R_hat_100 <- gelman.preplot(mcmc.list.100, max.bins = 2000)$shrink[,1,1]


pdf("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/Rhat_comparison_10_100.pdf")
#ich muss die normalen plots machen!!! - aber wie dann übereinanderlegen???
plot(R_hat_10, type = "l", col = "blue", ylim = c(0,2), xlab = "iteration", ylab = "R hat", main = "R hat comparison")
dev.off()