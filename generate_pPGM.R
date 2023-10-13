rm(list = ls())

library("Matrix")
library("matrixcalc")
library("MASS")
library("igraph")
library("ggplot2")
library("dplyr")
#set.seed(1234)
# set the dimension/length of X (p int he thesis)
dim = 30


#set the parameter values
eta_0 = rep(1, 30)
#eta_0 = rep(1, 10) *runif(10, 0, 1)
#eta_1 = eta_0




#  set the number of samples
n_samples = 4

# Create a sample correlation matrix with negative values
# Generate a random correlation matrix with negative values
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

for (i in 1:dim){
  eta_i = eta_y[i]
  theta_ii =  0 #theta[i,i], we have theta_ii = 0 from the def of the model
  X_0[i,] = rpois(n_samples, exp(eta_i))
  
}

# iterations
iterations = 1000
log_likelihood_values = numeric(iterations)
eta_minus_i_old = matrix(eta_0, nrow = dim, ncol = n_samples)
eta_minus_i= matrix(NA, nrow = dim, ncol = n_samples)
log_likelihood_variance = numeric(iterations)
X_old = X_0
X_new = X_0
mean_X_new = matrix(NA, nrow = dim, ncol = iterations)
report = numeric(iterations)
for (i in 1:dim){
  name <- paste0("X_", i, "_simulations")
  assign(name, matrix(X_0[i,], nrow = n_samples, ncol = 1))
}
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
  #update X_old and eta_minus_i_minus_1
  #eta_minus_i_old = eta_minus_i
  #X_old = X_new
  
  #c) update
  
  # Step 4: caluclate the mean 
  mean_X_new[,loop] = rowMeans(X_new)
    if (loop>1){
    report[loop] = mean(X_new)

    #save the simulated values
    for(i in 1:dim){
      name <- paste0("X_", i, "_simulations")
      assign(name, cbind(get(name), X_new[i,]))
    }
  }
  #log_likelihood_samples = numeric(n_samples)
  
  #for(j in 1:n_samples){
  #  log_likelihood_samples[j] = t(eta_0) %*% X_new[,j]+ 0.5 * t(X_new[,j]) %*% theta %*% X_new[,j]+ t(rep(1, dim)) %*% -log(factorial(X_new[,j])) 
  #}
  #log_likelihood_values[loop] = sum(log_likelihood_samples)
  #print(log_likelihood_values[loop])
  #log_likelihood_variance[loop] = var(log_likelihood_samples)

}


#pdf(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/likelihood.pdf")
#plot(log_likelihood_values, type = "l")  
#plot(log_likelihood_variance, type = "l") 
#plot(density(X_new[10, 250:1000]))
#dev.off()
picture_file  = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/pPGM_simulation_dim", dim, "_lower_bound_", lower_bound, "_upper_bound_", upper_bound, "_n_samples_", n_samples, "_iterations_", iterations,".pdf")

pdf(file = picture_file)
#das evtl als ggplot schreiben für alle variablen
for(i in 1:dim){
  plot(mean_X_new[i,], type = "l")
}
plot(report, type = "l")
dev.off()


exp = matrix(NA, nrow = dim, ncol = n_samples)
sample_exp = numeric(dim)
final_exp = numeric(dim)
for (i in 1:dim){
  for(j in 1:n_samples){
    exp[i,j] = exp(eta_0[i] + theta[i,] %*% X_new[,j])
  }
  final_exp[i] = mean(exp[i,])
  sample_exp[i] = mean(X_new[i,])
}


save(X_new, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/pPGM_simulation_dim", dim, "_lower_bound_", lower_bound, "_upper_bound_", upper_bound, "_n_samples_", n_samples, "_iterations_", iterations,".RData"))


###gelman-rubin statistics

pdf(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/chain_test.pdf")
plot(X_1_simulations[1,], type = "l")
dev.off()

#calculate the test statistic (mean)
psi <- t(apply(X_1_simulations, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))

Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}

print(Gelman.Rubin(psi))




pdf(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/chain_test.pdf")
#plot psi for the four chains
par(mfrow=c(2,2))
for (i in 1:4)
plot(psi[i, 1:iterations], type="l",
xlab=i, ylab=bquote(psi))
par(mfrow=c(1,1)) #restore default

#plot the sequence of R-hat statistics
rhat <- rep(0, iterations)
for (j in 1:iterations)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[1:iterations], type="l", xlab="", ylab="R")
abline(h=1.1, lty=2)
dev.off()





#Gelman Rubin on the mean
mean_X_new

Gelman.Rubin(mean_X_new)

pdf(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/chain_test_mean.pdf")
#plot mean_X_new for the four chains
par(mfrow=c(2,2))
for (i in 1:4)
plot(mean_X_new[i, 1:iterations], type="l",
xlab=i, ylab=bquote(mean_X_new))
par(mfrow=c(1,1)) #restore default

#plot the sequence of R-hat statistics
rhat <- rep(0, iterations)
for (j in 1:iterations)
rhat[j] <- Gelman.Rubin(mean_X_new[,1:j])
plot(rhat[1:iterations], type="l", xlab="", ylab="R")
abline(h=1.1, lty=2)
dev.off()


