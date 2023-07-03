#this code is leant on what is done in the paper Claudia van Borkulo et al. “Comparing Network Structures on Three Aspects: A Permutation Test”. In:
#Psychological Methods (Dec. 2021)

#ich probiere das umzuschreiben dass positiv und negativ gewichtet wird


library("Matrix")
library("matrixcalc")
library("MASS")
library("igraph")


#should I set a seed??

#number of edges
number_edges = 20


#probability of an edge
edge_probability = 0.5


#swet the upper and the lower bound for the uniform distribution the edge weights are drawn from
upper_bound = 3
lower_bound = 2

#set the sample size
sample_size = 250

#number of repetitions
n_reps = 100

#größten und kleinsten korr in tabelle
#von 100 auf 200 perm gehen
#verbal resultate beschreiben


#this part is from NCT
simulate_graph_data <- function(d, n, adj_matrix=NA, weighted_adj=NA, prob, u=0.1,
                       upper.bound=.9, lower.bound=.3, n_reps = 100){

    #simulates a graph
    if(any(is.na(adj_matrix)) & any(is.na(weighted_adj))){
      #samples a random graph with d nodes
      #every possible edge is present with probability prob
      adj_matrix.ig <- erdos.renyi.game(d,prob)
      #get adjacency matrix
      adj_matrix <- as.matrix(get.adjacency(adj_matrix.ig))
      #simulates a weighted_adj matrix, this is a weighted version of the adjecency matrix
      weighted_adj = matrix(runif(d^2, lower.bound, upper.bound), d, d) * adj_matrix



# Generate a binary matrix indicating whether each element should be positive or negative
sign_matrix <- matrix(sample(c(-1, 1), d^2, replace = TRUE), d, d)

# Multiply the matrices element-wise, accounting for positive or negative signs
weighted_adj <- weighted_adj * sign_matrix 
      #diag is set to 0
      diag(weighted_adj) = 0
      omega = weighted_adj
      #diag is set to the smallest eigenvalue and 0.1 + u is added
      diag(omega) = abs(min(Re(eigen(omega)$values))) + 1
      #omega is forced to be symmetric
      omega <- as.matrix(forceSymmetric(omega))
      #if omega is positive definite, break
      if(!is.positive.definite(omega)) adj_matrix=weighted_adj=NA 
    }

  
  #converts covariance matrix to correlation matrix
  # prec2cov checkt ob omega invertierbar ist und falls nicht gibt es eine Fehlermeldung, falls ja wird die Inverse berechnet
  sigma = cov2cor(SpiecEasi::prec2cov(omega))
  #generates random sample from multivariate normal distribution with length n and mena zero (rep(0,d)) and covariance matrix sigma
  gen_data1 <- list()
  gen_data2 <- list()
 for (i in 1: n_reps){
 x1 = mvrnorm(n, rep(0, d), sigma)
 x2 = mvrnorm(n, rep(0, d), sigma)
 gen_data1[[i]] <- x1
 gen_data2[[i]] <- x2
 }

 #sigmahat = cor(x1)
 results <- list(data_1=gen_data1, data_2 = gen_data2, sigma=sigma, omega=omega, adj_matrix=adj_matrix, weighted.adj = weighted_adj)
 return(results)
}

simulation = simulate_graph_data(d = number_edges, prob = edge_probability, n = sample_size, upper.bound = upper_bound, lower.bound = lower_bound)
round(simulation$sigma,3)
round(simulation$omega, 3)
simulation$weighted.adj #das mit den negativen Werten scheint noch nicht zu klappen

#pdf(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/graphs.pdf")
#plot(graph1,edge.label = E(graph1)$weight)
#plot(graph2,edge.label = E(graph2)$weight)
#dev.off()

save(simulation, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations/simulation_", number_edges, "_edges_", edge_probability, "_prob_", lower_bound, "_", upper_bound, "_bounds_", sample_size, "_sample_size_", n_reps, "_repetitions.RData"))

simulation$adj_matrix
pracma::cond(simulation$sigma)
min(simulation$sigma)
round(simulation$omega,3)
