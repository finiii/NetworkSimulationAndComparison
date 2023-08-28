#this code is leant on what is done in the paper Claudia van Borkulo et al. “Comparing Network Structures on Three Aspects: A Permutation Test”. In:
#Psychological Methods (Dec. 2021)

#ich probiere das umzuschreiben dass positiv und negativ gewichtet wird


library("Matrix")
library("matrixcalc")
library("MASS")
library("igraph")
library("SpiecEasi")

### same parameters for both graphs

#number of edges
number_edges = 10

#probability of an edge
edge_probability = 0.3

#es muss immer der gleiche graph sein, ansonsten kann man es ja nicht verlgeichen!!
     adj_matrix.ig <- erdos.renyi.game(number_edges,edge_probability)
    #get adjacency matrix
    adjacency_matrix <- as.matrix(get.adjacency(adj_matrix.ig))


#swet the upper and the lower bound for the uniform distribution the edge weights are drawn from

# graph 1
upper_bound_1 = 0.3
lower_bound_1 = 0.2

# graph 2
upper_bound_2 = 0.8
lower_bound_2 = 0.7

#set the sample size, same for both graphs
sample_size = 250

#number of repetitions, same for both graphs
n_reps = 100

#größten und kleinsten korr in tabelle
#von 100 auf 200 perm gehen
#verbal resultate beschreiben



#this part is from NCT
simulate_graph_data <- function(d, n, adj_matrix=adjacency_matrix, weighted_adj=NA, prob,
                       upper_bound_1=.9, lower_bound_1=.3,upper_bound_2=.9, lower_bound_2=.3, n_reps = 100){




      #samples a random graph with d nodes
      #every possible edge is present with probability prob, this is the same for both graphs
 
      #simulates a weighted_adj matrix, this is a weighted version of the adjecency matrix
      weighted_adj_1 = matrix(runif(d^2, lower_bound_1, upper_bound_1), d, d) * adj_matrix

      #now we have a second weighted matrix with weights from another distribution
      weighted_adj_2 = matrix(runif(d^2, lower_bound_2, upper_bound_2), d, d) * adj_matrix

      #diag is set to 0
      diag(weighted_adj_1) = 0
      diag(weighted_adj_2) = 0
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
 x1 = mvrnorm(n, rep(0, d), sigma_1)
 x2 = mvrnorm(n, rep(0, d), sigma_2)
 gen_data1[[i]] <- x1
 gen_data2[[i]] <- x2
 }

 #sigmahat = cor(x1)
 results <- list(data_1=gen_data1, data_2 = gen_data2, sigma_1=sigma_1, omega_1=omega_1, adj_matrix=adj_matrix, weighted.adj_1 = weighted_adj_1, sigma_2=sigma_2, omega_2=omega_2, weighted.adj_2 = weighted_adj_2)
 return(results)
}

simulation = simulate_graph_data(d = number_edges, prob = edge_probability, n = sample_size, upper_bound_1=upper_bound_1, lower_bound_1=lower_bound_1,upper_bound_2=upper_bound_2, lower_bound_2=lower_bound_2, adj_matrix = adjacency_matrix)

save(simulation, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_different_weights/simulation_", number_edges, "_edges_", edge_probability, "_prob_", lower_bound_1, "_", upper_bound_1,"_", lower_bound_2, "_", upper_bound_2, "_bounds_", sample_size, "_sample_size_", n_reps, "_repetitions.RData"))


