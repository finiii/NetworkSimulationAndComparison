#this code is leant on what is done in the paper Claudia van Borkulo et al. “Comparing Network Structures on Three Aspects: A Permutation Test”. In:
#Psychological Methods (Dec. 2021)

library(NetworkComparisonTest)
library("Matrix")
library("matrixcalc")
library("MASS")
library("igraph")


#this part is from NCT
simulate_graph_data <- function(d, n, adj_matrix=NA, weighted_adj=NA, prob, v=0.1, u=0.1,
                       upper.bound=.9, lower.bound=.3){
  # simulates graph
  count=0
  repeat{
    #simulates a graph
    if(any(is.na(adj_matrix)) & any(is.na(weighted_adj))){
      #samples a random graph with d nodes
      #every possible edge is present with probability prob
      adj_matrix.ig <- erdos.renyi.game(d,prob)
      #get adjacency matrix
      adj_matrix <- as.matrix(get.adjacency(adj_matrix.ig))
      #simulates a weighted_adj matrix, this is a weighted version of the adjecency matrix
      weighted_adj = matrix(runif(d^2, lower.bound, upper.bound), d, d) * adj_matrix
      #diag is set to 0 
      diag(weighted_adj) = 0
      #everything is multiplied by v (what is this v?)
      omega = weighted_adj * v
      #diag is set to the smallest eigenvalue and 0.1 + u is added
      diag(omega) = abs(min(Re(eigen(omega)$values))) + 0.1 + u
      #omega is forced to be symmetric
      omega <- as.matrix(forceSymmetric(omega))
      #if omega is positive definite, break
      if(!is.positive.definite(omega)) adj_matrix=weighted_adj=NA else break
    }
    #if adj_matrix and weighted_adj are specified
    if(!any(is.na(adj_matrix))& !any(is.na(weighted_adj))){
      omega = weighted_adj * v
      diag(omega) = abs(min(Re(eigen(omega)$values))) + 0.1 + u
      omega <- as.matrix(forceSymmetric(omega))
      if(is.positive.definite(omega)) break
    }
    #if adj_matrix is specified and weighted_adj is not
    if(!any(is.na(adj_matrix))& any(is.na(weighted_adj))){
      weighted_adj = matrix(runif(d^2, lower.bound, upper.bound), d, d) * adj_matrix
      omega = weighted_adj * v
      diag(omega) = abs(min(Re(eigen(omega)$values))) + 0.1 + u
      omega <- as.matrix(forceSymmetric(omega))
      if(!is.positive.definite(omega)) weighted_adj=NA else break
    }
    count=count+1
    print(count)
  }
  #converts covariance matrix to correlation matrix
  sigma = cov2cor(solve(omega))
  #calculates the inverse of sigma, warum wird das gemacht?
  omega = solve(sigma)
  #generates random sample from multivariate normal distribution with length n and mena zero (rep(0,d)) and covariance matrix sigma
  x = mvrnorm(n, rep(0, d), sigma)
  sigmahat = cor(x)
  results <- list(data=x, sigma=sigma, weighted_adj=weighted_adj, adj_matrix=adj_matrix)
  #returns the generated data, the covariance matrix sigma, the adjacency matrix adj_matrix and the weighted_adj matrix (which is equal to the weighted adjacency matrix)
  return(results)
}

set.seed(050623)

simulation1 <- simulate_graph_data(d = 10, n = 500, prob = 0.2)


graph1 <- graph.adjacency(simulation1$weighted_adj, weighted = TRUE, mode = "undirected")
graph1 <- graph_from_adjacency_matrix(
  simulation1$weighted_adj,
  mode = "undirected",
  weighted = TRUE
)

#wie simuliert man den zweiten graphen? 
#erste möglichkeit: andere prob
simulation2 <- simulate_graph_data(d = 10, n = 500, prob = 0.4)

graph2 <- graph.adjacency(simulation2$weighted_adj, weighted = TRUE, mode = "undirected")
graph2 <- graph_from_adjacency_matrix(
  simulation2$weighted_adj,
  mode = "undirected",
  weighted = TRUE
)

#man könnte auch die adjazenzmatrix von dem ersten graphen nehmen und die verändern


pdf(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/graphs.pdf")
plot(graph1,edge.label = E(graph1)$weight)
plot(graph2,edge.label = E(graph2)$weight)
dev.off()

save(simulation1, file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations/simulation1.RData")
save(simulation2, file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations/simulation2.RData")