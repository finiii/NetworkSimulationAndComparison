rm(list = ls())

library("Matrix")
library("matrixcalc")
library("MASS")
library("igraph")
#set.seed(1234)
#set the parameter values
eta_0 = rep(3, 10)
#eta_0 = rep(1, 10) *runif(10, 0, 1)
#eta_1 = eta_0

# set the dimension/length of X
dim = 10

#  set the number of samples
n_samples = 5000

# Create a sample correlation matrix with negative values
# Generate a random correlation matrix with negative values
theta <- matrix(0, nrow = dim, ncol = dim)

# Fill the upper triangular part with random negative values
for (i in 1:(dim - 1)) {
  for (j in (i + 1):dim) {
    theta[i, j] <- runif(1, -20, -10)  # Random value
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


#swet the upper and the lower bound for the uniform distribution the edge weights are drawn from
#lower_bound = -0.5
#upper_bound = 0

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
  X_old = X_new
  
  #c) update
  
  # Step 4: caluclate the likelihood
  log_likelihood_samples = numeric(n_samples)
  
  for(j in 1:n_samples){
    log_likelihood_samples[j] = t(eta_0) %*% X_new[,j]+ 0.5 * t(X_new[,j]) %*% theta %*% X_new[,j]+ t(rep(1, dim)) %*% -log(factorial(X_new[,j])) 
  }
  log_likelihood_values[loop] = sum(log_likelihood_samples)
  print(log_likelihood_values[loop])
  log_likelihood_variance[loop] = var(log_likelihood_samples)

}


pdf(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/likelihood.pdf")
plot(log_likelihood_values, type = "l")  
plot(log_likelihood_variance, type = "l") 
plot(density(X_new[10, 250:1000]))
dev.off()


mean(X_new[1,])
mean(X_new[2,])
mean(X_new[3,])
mean(X_new[4,])
mean(X_new[5,])

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

final_exp - sample_exp

 #exp(eta_0 + theta %*% X_new[,4])