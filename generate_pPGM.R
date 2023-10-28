rm(list = ls())


#load the required packages
library("Matrix")
library("matrixcalc")
library("MASS")
library("igraph")
library("ggplot2")
library("dplyr")
#for the diagnostics
library(coda)



# set the dimension/length of X (p int he thesis)
dim = 50

#set the parameter values
eta_0 = rep(1, dim)
#eta_0 = rep(1, 10) *runif(10, 0, 1)
#eta_1 = eta_0

#  set the number of chains
n_samples = 6

# Create a empty correlation matrix with negative values
theta <- matrix(0, nrow = dim, ncol = dim)

#swet the upper and the lower bound for the uniform distribution the edge weights are drawn from
lower_bound = -20
upper_bound = -10

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
iterations = 100000  #N
burn_in = iterations/2 # burn in length
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

simulation_name <- paste0("pPGM_simulation_dim", dim, "_lower_bound_", lower_bound, "_upper_bound_", upper_bound, "_n_samples_", n_samples, "_iterations_", iterations)
#save the results
result = list(X_new = X_new, mean_X_new = mean_X_new)
save(result, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/", simulation_name,".RData"))

#I took this function from the coda package and changed it a bit
#I think something is still off here since it does not coincide with the plots.....
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


###prepare for R hat with coda
for (i in 1:dim){
      name <- paste0("X_", i, "_simulations")
      simulated_data <- get(name)
      mcmc.list <- mcmc.list(as.mcmc(simulated_data[1,]), as.mcmc(simulated_data[2,]), as.mcmc(simulated_data[3,]) ,as.mcmc(simulated_data[4,]))
      new_name <- paste0("mcmc.list.X", i)
      assign(new_name, mcmc.list)
}


#R hat per iteration and dimension
for (i in 1:dim){
  name <- paste0("X_", i, "_Rhat")
  assign(name, rep(NA, iterations))
}


for (i in 1:dim){
  name <- paste0("X_", i, "_Rhat")
    mcmc.list <- paste0("mcmc.list.X", i)
    mcmc.list <- get(mcmc.list)
    assign(name, as.vector(gelman.preplot(mcmc.list, max.bins = 2000)$shrink[,1,1]))
    print(name)
}

#build the Rhat convergence check
convergence_check_Rhat <- function(limit = 1.1){
  non_converged_chains <- list()
  for (i in 1:dim){
    name <- paste0("X_", i, "_Rhat")
    Rhat <- get(name)
    #I discard the first few Rhat values since they do not use as many values already
    #I choose to look at the last thrid of R hat values that should be enough
      if (any(Rhat[-c(1:(2*iterations/3))] > limit,na.rm = TRUE)){
        sum_bigger <- sum(Rhat > limit, na.rm = TRUE)
       message <- paste0("X_", i, " did not converge ", sum_bigger," Rhat values are not smaller than ", limit)
        non_converged_chains <- c(non_converged_chains, message)
      }
    }
  if (length(non_converged_chains) == 0) {
    return("All chains converged.")
  } else {
    return(paste("Chains that did not converge:", non_converged_chains))
  }
  }

convergence_check_Rhat(limit = 1.1)

#build the acf convergence check
convergence_check_acf <- function(burn_in, limit = 0.05){
  non_converged_chains <- list()
  for (i in 1:dim){
    name <- paste0("X_", i, "_simulations")
    simulated_data <- get(name)
    for (j in 1:n_samples){
      acf_values <- acf(simulated_data[j, -c(1:burn_in)], plot = FALSE)$acf
      if (any(acf_values[length(acf_values)] > limit)){
        max_lag <- length(acf_values)
              message <- paste0("X_", i, " chain ", j, " did not converge, the last acf value, lag ", max_lag, " is ", acf_values[length(acf_values)], " and should be smaller than ", limit)
        non_converged_chains <- c(non_converged_chains, message)
      }
    }
  }
  
  if (length(non_converged_chains) == 0) {
    return("All chains converged.")
  } else {
    return(paste("Chains that did not converge:", non_converged_chains))
  }
}
convergence_check_acf(burn_in)

#diagnostic plots
pdf(file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/diagnostic_plots/", simulation_name,".pdf"))



#density plots
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



#plot the acf functions
    for(i in 1:dim){
      name <- paste0("X_", i, "_simulations")
      simulated_data <- get(name)
      for (j in 1:n_samples)  {
        acf(simulated_data[j, -c(1:burn_in)], main = paste0("X_", i, " chain ",j))
      }
    }



#plot the r hat diagnostics
for(i in 1:dim){
  mcmc.list <- paste0("mcmc.list.X", i)
  mcmc.list <- get(mcmc.list)
  gelman.plot(mcmc.list, autoburnin = T) #, transform = TRUE
}

dev.off()

txt_file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/diagnostic_plots/", simulation_name,".txt")
sink(txt_file)
#here I print the automatic diagnostics
print("Convergence check with the Rubin-Gelman diagnostic")
print(convergence_check_Rhat(limit = 1.1))
print("Convergence check with the acf diagnostic")
print(convergence_check_acf(burn_in))
sink()

