rm(list = ls())


#load the required packages
library("Matrix")
library("matrixcalc")
library("MASS")
library("igraph")


# set the dimension/length of X (p in the thesis)
dim = 50

#set the parameter values
eta_0 = rep(5, dim)
#eta_0 = rep(1, 10) *runif(10, 0, 1)
#eta_1 = eta_0

#  set the number of chains
n_samples = 6

# Create a empty correlation matrix with negative values
theta <- matrix(0, nrow = dim, ncol = dim)

#swet the upper and the lower bound for the uniform distribution the edge weights are drawn from
lower_bound = -0.01
upper_bound = -0.005


#probability of an edge
edge_probability = 0.8

#samples a random graph with d nodes
#every possible edge is present with probability prob
#
adj_matrix.ig <- erdos.renyi.game(dim,edge_probability)

adj_matrix <- as.matrix(get.adjacency(adj_matrix.ig))

#simulates a weighted_adj matrix, this is a weighted version of the adjecency matrix
weighted_adj = matrix(runif(dim^2, lower_bound, upper_bound), dim, dim) * adj_matrix

#diag is set to 0
diag(weighted_adj) = 0
theta = weighted_adj
theta = as.matrix(forceSymmetric(theta))

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
  X_0[i,] = rpois(n_samples, exp(eta_i))
}

# iterations
iterations = 500000  #N




X_old = X_0
X_new = X_0
simulation_name <- paste0("pPGM_simulation_dim", dim, "_lower_bound_", lower_bound, "_upper_bound_", upper_bound, "_n_samples_", n_samples, "_iterations_", iterations, "_edge_probability_", edge_probability)





input_function = function(eta_0, theta, X_0, iterations, n_samples, dim, X_new){
#here the gibbs sampling starts
for (i in 1:dim){
  name <- paste0("X_", i, "_simulations")
  assign(name, matrix(X_0[i,], nrow = n_samples, ncol = 1))
}

eta_minus_i_old = matrix(eta_0, nrow = dim, ncol = n_samples)
eta_minus_i= matrix(NA, nrow = dim, ncol = n_samples)
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
    #here the eta update happens
    #eta_minus_i_old = eta_minus_i
  }
  

  
  # Step 4: caluclate the correlation matrix


    #save the simulated values
    for(i in 1:dim){
      name <- paste0("X_", i, "_simulations")
      assign(name, cbind(get(name), X_new[i,]))
    }
  }

result = list(X_new = X_new, theta = theta)
#add all X_i_simulations to results
for (i in 1:dim){
  name <- paste0("X_", i, "_simulations")
  result[[name]] <- get(name)
}

return(result)
}


library(foreach)
library(doParallel)
library(doRNG)
library(PRROC)
library("Matrix")
library("matrixcalc")
library("MASS")
library("igraph")
library(pulsar)
library(parallel)
no_cores <- 30#55
cl <- makeCluster(no_cores, outfile = "")
clusterEvalQ(cl, {library("Matrix")
library("matrixcalc")
library("MASS")
library("igraph")
  library(PRROC)
  library(pulsar)
  library(parallel)})
registerDoParallel(cl)

fct <- function(i){
  out <- input_function(eta_0, theta, X_0, iterations, n_samples, dim,X_new)
  return(out)
}

reps = 1

time_parallel <- system.time(
  temp <- foreach(i=1:reps) %dopar% {fct(i)})
stopCluster(cl)
time_parallel

save(temp, time_parallel, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_poisson_final/", simulation_name, ".RData"))



