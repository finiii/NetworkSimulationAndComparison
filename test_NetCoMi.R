library(NetCoMi)


#load my own data to tst if it works with them
simulation_path <-"/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulation_change_edge/simulation_changed_biggest_edge_10_edges_0.3_prob_0.3_0.4_bounds_250_sample_size_100_repetitions.RData"
load(simulation_path)


Yg1 <- simulation$data_1[[1]]
Yg2 <- simulation$data_2[[1]]


# try to fit the network

net_season <- netConstruct(data = Yg1, 
                           data2 = Yg2)