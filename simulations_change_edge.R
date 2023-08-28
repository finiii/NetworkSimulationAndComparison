#this code is leant on what is done in the paper Claudia van Borkulo et al. “Comparing Network Structures on Three Aspects: A Permutation Test”. In:
#Psychological Methods (Dec. 2021)

#ich probiere das umzuschreiben dass positiv und negativ gewichtet wird


library("Matrix")
library("matrixcalc")
library("MASS")
library("igraph")
library("SpiecEasi")

load("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_different_weights/simulation_10_edges_0.3_prob_0.2_0.3_0.2_0.3_bounds_250_sample_size_100_repetitions.RData")

### same parameters for both graphs

    #get adjacency matrix
    adjacency_matrix <- simulation$adj_matrix

number_edges = 10 
edge_probability = 0.3

#swet the upper and the lower bound for the uniform distribution the edge weights are drawn from

# graph 1
upper_bound_1 = 0.9
lower_bound_1 = 0.8


#set the sample size, same for both graphs
sample_size = 250

#number of repetitions, same for both graphs
n_reps = 100

#größten und kleinsten korr in tabelle
#von 100 auf 200 perm gehen
#verbal resultate beschreiben






      #samples a random graph with d nodes
      #every possible edge is present with probability prob, this is the same for both graphs
 
      #simulates a weighted_adj matrix, this is a weighted version of the adjecency matrix
      weighted_adj_1 = matrix(runif(number_edges^2, lower_bound_1, upper_bound_1), number_edges, number_edges) * adjacency_matrix
       #diag is set to 0
      diag(weighted_adj_1) = 0
      max_index <- which.max(weighted_adj_1)
    weighted_adj_2 = weighted_adj_1
    # Set the largest element to zero
    weighted_adj_2[max_index] <- 0 #kann es vlt auch auf 0.5 mal das setzen oder so 
      

     
      omega_1 = weighted_adj_1
      omega_2 = weighted_adj_2
      #diag is set to the smallest eigenvalue and 0.1 is added
      diag(omega_1) = abs(min(Re(eigen(omega_1)$values))) + 0.1
        diag(omega_2) = abs(min(Re(eigen(omega_2)$values))) + 0.1
      #omega is forced to be symmetric
      omega_1 <- as.matrix(forceSymmetric(omega_1))
      omega_2 <- as.matrix(forceSymmetric(omega_2)) 
      #if omega is positive definite, break
      if(!is.positive.definite(omega_1)) adj_matrix=weighted_adj_1=NA 
      if(!is.positive.definite(omega_2)) adj_matrix=weighted_adj_2=NA 
    

  
  #converts covariance matrix to correlation matrix
  # prec2cov checkt ob omega invertierbar ist und falls nicht gibt es eine Fehlermeldung, falls ja wird die Inverse berechnet
  sigma_1 = cov2cor(SpiecEasi::prec2cov(omega_1))
  sigma_2 = cov2cor(SpiecEasi::prec2cov(omega_2))
  #generates random sample from multivariate normal distribution with length n and mena zero (rep(0,d)) and covariance matrix sigma
  gen_data1 <- list()
  gen_data2 <- list()
 for (i in 1: n_reps){
 x1 = mvrnorm(sample_size, rep(0, number_edges), sigma_1)
 x2 = mvrnorm(sample_size, rep(0, number_edges), sigma_2)
 gen_data1[[i]] <- x1
 gen_data2[[i]] <- x2
 
 }
 #sigmahat = cor(x1)
 results <- list(data_1=gen_data1, data_2 = gen_data2, sigma_1=sigma_1, omega_1=omega_1, adj_matrix=adjacency_matrix, weighted.adj_1 = weighted_adj_1, sigma_2=sigma_2, omega_2=omega_2, weighted.adj_2 = weighted_adj_2)



simulation = results
save(simulation, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulation_change_edge/simulation_changed_biggest_edge_", number_edges, "_edges_", edge_probability, "_prob_", lower_bound_1, "_", upper_bound_1, "_bounds_", sample_size, "_sample_size_", n_reps, "_repetitions.RData"))


