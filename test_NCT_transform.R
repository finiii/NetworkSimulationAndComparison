library(NetworkComparisonTest)
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



cor(data1)
cor(data2)
#bring them on the same length
data2 <- data2[(nrow(data2)-nrow(data1)):nrow(data2),]
rownames(data2) <- 1:nrow(data2)

#die daten transformieren

data1 <- clr(data1)
data2 <- clr(data2)

test <- NCT(data1, data2, it = 1000)

##############################################
#test the null hypothesis
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

#cut the data into slices with 1000 simulations each
#in this case we do not need to account for acf bc we use the simulations under the independent model
n_reps <- 100
NCT_pval <- numeric(n_reps)
BGGM_pval <- numeric(n_reps)
for(iteration in 1:n_reps)
{
  data1 <- data_null[(iteration*1000)-1000:iteration*1000,]
  data2 <- data_null_2[(iteration*1000)-1000:iteration*1000,]

  #do the NCT
  test_NCT <- NCT(data1, data2, it = 1000)
  NCT_pval[iteration] <- test_NCT$nwinv.pval

  ### do the BGGM
  test_BGGM <- BGGM::ggm_compare_ppc(data1, data2,
                       iter = 1000)

BGGM_pval[iteration] <- test_BGGM[["ppp_jsd"]]
}

result_NCT_BGGM <- data.frame(NCT_pval, BGGM_pval)
result_NCT_BGGM_name <- paste0(gsub("\\.RData$", "", basename(result1)),"_", gsub("\\.RData$", "", basename(result2)), "_NCT_BGGM.RData")
save(result_NCT_BGGM, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/diagnostic_plots/results_test_NCT_BGGM/",result_NCT_BGGM_name))