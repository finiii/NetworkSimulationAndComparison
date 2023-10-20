rm(list = ls())


#load the required packages
library("Matrix")
library("matrixcalc")
library("MASS")
library("igraph")
library("ggplot2")
library("dplyr")



# set the dimension/length of X (p int he thesis)
dim = 10

#set the parameter values
eta_0 = rep(1, dim)
#eta_0 = rep(1, 10) *runif(10, 0, 1)
#eta_1 = eta_0

#  set the number of chains
n_samples = 4

# Create a empty correlation matrix with negative values
theta <- matrix(0, nrow = dim, ncol = dim)

#swet the upper and the lower bound for the uniform distribution the edge weights are drawn from
lower_bound = -1
upper_bound = 0

# Fill the upper triangular part with random negative values
for (i in 1:(dim - 1)) {
  for (j in (i + 1):dim) {
    theta[i, j] <- runif(1, lower_bound, upper_bound)  # Random value
    theta[ j,i]=theta[i, j]
  }
}

# Copy the upper triangular values to the lower triangular part
theta[lower.tri(theta)] <- t(theta)[lower.tri(theta)]

# Set the diagonal elements to 0
diag(theta) <- 0 

#generate a matrix from a random graph and weight it

#probability of an edge
#edge_probability = 0.5

#samples a random graph with d nodes
#every possible edge is present with probability prob
#adj_matrix.ig <- erdos.renyi.game(dim,edge_probability)
#get adjacency matrix
#adj_matrix <- as.matrix(get.adjacency(adj_matrix.ig))

#simulates a weighted_adj matrix, this is a weighted version of the adjecency matrix
#weighted_adj = matrix(runif(dim^2, lower_bound, upper_bound), dim, dim) * adj_matrix

#diag is set to 0
#diag(weighted_adj) = 0
#theta = weighted_adj
#theta = as.matrix(forceSymmetric(theta))



##########################################################################

# Step 1: Compute eta_y
# since we have a fixed eta, we do not need to compute it 
eta_y = eta_0

# Step 2: generate X_0 from the independent model with eta_y

# write an empty matrix for the samples
X_0 = matrix(NA, nrow = dim, ncol = n_samples)

#simulate the starting values for the samples
for (i in 1:dim){
  eta_i = eta_y[i]
  theta_ii =  0 #theta[i,i], we have theta_ii = 0 from the def of the model
  X_0[i,] = rpois(n_samples, exp(eta_i))
  
}

# iterations
iterations = 5000  #N
burn_in = 200 # burn in length
log_likelihood_values = numeric(iterations)
eta_minus_i_old = matrix(eta_0, nrow = dim, ncol = n_samples)
eta_minus_i= matrix(NA, nrow = dim, ncol = n_samples)
log_likelihood_variance = numeric(iterations)
X_old = X_0
X_new = X_0
mean_X_new = matrix(NA, nrow = dim, ncol = iterations)
mean_per_chain = matrix(NA, nrow = n_samples, ncol = iterations)
report = numeric(iterations)
for (i in 1:dim){
  name <- paste0("X_", i, "_simulations")
  assign(name, matrix(X_0[i,], nrow = n_samples, ncol = 1))
}

#here the gibbs sampling starts
for(loop in 1:iterations){
  
  #step 3
  #a) compute the conditional parameters
  for(j in 1:n_samples){
    for (i in 1:dim){
      #i-th row of theta after removing the i-th column from theta
      theta_minus_i = theta[i,-i]
      eta_minus_i[i,j] = eta_minus_i_old[i,j] + theta_minus_i %*% X_new[-i,j]
      
      X_new[i,j] = rpois(1, exp(eta_minus_i[i,j]))
    }


  }
  

  
  # Step 4: caluclate the mean 
    mean_X_new[,loop] = rowMeans(X_new)

    mean_per_chain[,loop] = colMeans(X_new)

    if (loop>1){
    report[loop] = mean(X_new)}


    #save the simulated values
    for(i in 1:dim){
      name <- paste0("X_", i, "_simulations")
      assign(name, cbind(get(name), X_new[i,]))
    }
  }


#save the results
result = list(X_new = X_new, mean_X_new = mean_X_new)
save(result, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/pPGM_simulation_dim", dim, "_lower_bound_", lower_bound, "_upper_bound_", upper_bound, "_n_samples_", n_samples, "_iterations_", iterations,".RData"))

  #delete the samples from the burn in phase

  X_new = X_new[, -c(1:burn_in)]
  mean_X_new = mean_X_new[, -c(1:burn_in)]
  mean_per_chain = mean_per_chain[, -c(1:burn_in)]
  for(i in 1:dim){
      name <- paste0("X_", i, "_simulations")
      simulated_data <- get(name)
      assign(name, simulated_data[, -c(1:burn_in)])
  }

#hier plotte ich alle running means der variablen über die chains 
#also x.B. mean_X_new[1,] ist der running mean der ersten variable über die 4 chains
picture_file  = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/pPGM_simulation_dim", dim, "_lower_bound_", lower_bound, "_upper_bound_", upper_bound, "_n_samples_", n_samples, "_iterations_", iterations,".pdf")

pdf(file = picture_file)
for(i in 1:dim){
  plot(mean_X_new[i,], type = "l")
}
plot(report, type = "l")
dev.off()


#density plots
pdf(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/density_plots.pdf")
for(i in 1:dim){
  name <- paste0("X_", i, "_simulations")
  simulated_data <- as_tibble(get(name))
  #but which iterations do we actually take??
  simulated_data <- mutate(simulated_data, chain = 1:n_samples)
  simulated_data <- tidyr::pivot_longer(simulated_data, cols = -chain, names_to = "variable", values_to = "value")
  
  plot <- ggplot2::ggplot(simulated_data, aes(x= value, fill = as.factor(chain))) +
    geom_density(alpha = 0.5) + ggtitle(paste0("Density plot of X_", i, " over the chains"))
  

  print(plot)
}
dev.off()

#plot the acf functions
pdf(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/acf_plots.pdf")
    for(i in 1:dim){
      name <- paste0("X_", i, "_simulations")
      simulated_data <- get(name)
      par(mfrow=c(2,2))
      for (j in 1:4)  {
        acf(simulated_data[j, ], main = paste0("X_", i, " chain ",j))
      }
      par(mfrow=c(1,1)) #restore default
    }
dev.off()


####gelman-rubin statistics
#calculate the test statistic (mean)
#test it here for X_1
phi <- t(apply(X_1_simulations, 1, cumsum))
for (i in 1:nrow(phi)){
phi[i,] <- phi[i,] / (1:ncol(phi))
}

Gelman.Rubin <- function(phi) {
# phi[i,j] is the statistic phi(X[i,1:j])
# for i-th chain and up to j-th iteration
phi <- as.matrix(phi)
n <- ncol(phi) # this is the number of iterations
m <- nrow(phi) # this is the number of chains
B <- n * var(rowMeans(phi)) #between variance estimate
W <- 0 #Initialize W, the within variance estimate
# Calculate row-wise means
mean_phi_i <- rowMeans(phi)
# Calculate the sum of squared differences
for (i in 1:m) {
  for (j in 1:n) {
    W <- W + (phi[i, j] - mean_phi_i[i])^2
  }
}

# Divide by (m * (n - 1))
W <- W / (m * (n - 1))
v_hat <- W*(n-1)/n + (B/n) #upper variance estimate
R_hat <- sqrt(v_hat / W) #Gelman-Rubin statistic
return(R_hat)
}
print(Gelman.Rubin(phi))


#plot the sequence of R-hat statistics
pdf(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/gelman_rubin_statistic.pdf")
for(i in 1:dim){
      name <- paste0("X_", i, "_simulations")
      simulated_data <- get(name)
      psi <- t(apply(simulated_data, 1, cumsum))
      for (i in 1:nrow(psi)){
      psi[i,] <- psi[i,] / (1:ncol(psi))
      }
      rhat <- rep(0, (iterations-burn_in))
      for (j in 1:(iterations-burn_in)){
      rhat[j] <- Gelman.Rubin(psi[,1:j])
      }
      plot(rhat, type="l", xlab="", ylab="R")
      abline(h=1.1, lty=2)}
dev.off()





#Gelman Rubin on the mean

mean_per_chain_mod <- t(apply(mean_per_chain, 1, cumsum))
for (i in 1:nrow(mean_per_chain_mod)){
mean_per_chain_mod[i,] <- mean_per_chain_mod[i,] / (1:ncol(mean_per_chain_mod))}
Gelman.Rubin(mean_per_chain_mod)

pdf(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/gelman_rubin_mean.pdf")
#plot mean_X_new for the four chains
par(mfrow=c(2,2))
for (i in 1:4)
plot(mean_per_chain[i, ], type="l",
xlab=i, ylab=bquote(mean_X_new))
par(mfrow=c(1,1)) #restore default

#plot the sequence of R-hat statistics
rhat <- rep(0, iterations-burn_in)
for (j in 1:(iterations-burn_in)){

rhat[j] <- Gelman.Rubin(mean_per_chain_mod[,1:j])
}
plot(rhat, type="l", xlab="", ylab="R hat")
abline(h=1.1, lty=2)
dev.off()


###R hat with coda

library(coda)
#this needs to be done for every variable
for (i in 1:dim){
      name <- paste0("X_", i, "_simulations")
      simulated_data <- get(name)
      mcmc.list <- mcmc.list(as.mcmc(simulated_data[1,]), as.mcmc(simulated_data[2,]), as.mcmc(simulated_data[3,]) ,as.mcmc(simulated_data[4,]))
      new_name <- paste0("mcmc.list.X", i)
      assign(new_name, mcmc.list)
}
mcmc.list.X1 <- mcmc.list(as.mcmc(X_1_simulations[1,]), as.mcmc(X_1_simulations[2,]), as.mcmc(X_1_simulations[3,]) ,as.mcmc(X_1_simulations[4,]))
gelman.diag(mcmc.list.X1)



pdf(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/gelman_rubin_coda.pdf")
for(i in 1:dim){
  mcmc.list <- paste0("mcmc.list.X", i)
  mcmc.list <- get(mcmc.list)
  gelman.plot(mcmc.list)
}
dev.off()