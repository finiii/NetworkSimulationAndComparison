library(BGGM)

#load my own data to test if it works with them
result <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/pPGM_simulation_dim10_lower_bound_-0.01_upper_bound_-0.005_n_samples_6_iterations_10000.RData"
load(result)
data1 <- lapply(names(result), function(name) {
  if (grepl("^X_[0-9]+_simulations$", name)) {
    as.vector(result[[name]][,seq(ncol(result[[name]])/2, ncol(result[[name]]), by = 35) ])
  }
})
data1 <- as.data.frame(do.call(cbind, data1))
rownames(data1) <- 1:nrow(data1)

#I need data dor the other group
result <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/pPGM_simulation_dim10_lower_bound_-0.1_upper_bound_-0.05_n_samples_6_iterations_50000.RData"
load(result)
data2 <- lapply(names(result), function(name) {
  if (grepl("^X_[0-9]+_simulations$", name)) {
    as.vector(result[[name]][,seq(ncol(result[[name]])/2, ncol(result[[name]]), by = 35) ])
  }
})
data2 <- as.data.frame(do.call(cbind, data2))



#bring them on the same length
data2 <- data2[(nrow(data2)-nrow(data1)):nrow(data2),]
rownames(data2) <- 1:nrow(data2)

#die daten transformieren

data1 <- clr(data1)
data2 <- clr(data2)

ppc <- BGGM::ggm_compare_ppc(data1, data2,
                       custom_obs = obs,
                       iter = 1000)

ppc