rm(list = ls())


library("Matrix")
library("matrixcalc")
library("MASS")
library("igraph")
library("SpiecEasi")


#load the file

simulation_path <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_different_weights/simulation_10_edges_0.3_prob_0.2_0.3_0.2_0.3_bounds_250_sample_size_100_repetitions.RData"
load(simulation_path)

n_reps = 100

#set the sample size, same for both graphs
sample_size = 1000

sigma_1 <- simulation$sigma_1
sigma_2 <- simulation$sigma_2

#generate
  gen_data1 <- list()
  gen_data2 <- list()
 for (i in 1: n_reps){
 x1 = mvrnorm(sample_size, rep(0, 10), sigma_1)
 x2 = mvrnorm(sample_size, rep(0, 10), sigma_2)
 gen_data1[[i]] <- x1
 gen_data2[[i]] <- x2
 }
  simulation$data_1= gen_data1
  simulation$data_2 = gen_data2


 #save the result 

 save(simulation, file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_different_weights/simulation_10_edges_0.3_prob_0.2_0.3_0.2_0.3_bounds_1000_sample_size_100_repetitions.RData")
