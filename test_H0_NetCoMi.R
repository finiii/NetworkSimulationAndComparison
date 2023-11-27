library(NetCoMi)


result1 <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/pPGM_simulation_dim10_lower_bound_0_upper_bound_0_n_samples_6_iterations_1e+05_edge_probability_0.8.RData"
load(result1)
data_null <- lapply(names(result), function(name) {
  if (grepl("^X_[0-9]+_simulations$", name)) {
    as.vector(result[[name]])
  }
})
data_null <- as.data.frame(do.call(cbind, data_null))
rownames(data_null) <- 1:nrow(data_null)

result2 <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/pPGM_simulation_dim10_lower_bound_0_upper_bound_0_n_samples_6_iterations_1e+05_edge_probability_0.5.RData"
load(result2)
data_null_2 <- lapply(names(result), function(name) {
  if (grepl("^X_[0-9]+_simulations$", name)) {
    as.vector(result[[name]])

  }
})
data_null_2 <- as.data.frame(do.call(cbind, data_null_2))
rownames(data_null_2) <- 1:nrow(data_null_2)

data1 <- data_null[1001:2000,]
data2 <- data_null_2[1001:2000,]


input_function <- function(data1, data2)
{

net <- netConstruct(data = data1, 
                           data2 = data2,
                           filtTax = "none",
                           filtSamp = "none",
             dataType = "counts")
netAna <- NetCoMi::netAnalyze(net, graphlet = F, gcmHeat = F)
test <- netCompare(netAna, permTest = TRUE,nPerm = 100L,
                            cores = 56L)
                            
}
save(test, file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/test_NetCoMi.RData")
