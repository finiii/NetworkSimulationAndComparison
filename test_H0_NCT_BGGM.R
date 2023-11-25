##############################################
#test the null hypothesis

library(NetworkComparisonTest)
library(BGGM)
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

data1 <- data_null[1:1000,]
data2 <- data_null_2[1:1000,]

input_function <- function(data1, data2)
{
#do the NCT
  test_NCT <- NCT(data1, data2, it = 1000)
  NCT_pval <- test_NCT$nwinv.pval

  ### do the BGGM
  test_BGGM <- BGGM::ggm_compare_ppc(data1, data2,
                       iter = 1000)

    BGGM_pval <- test_BGGM[["ppp_jsd"]]
  return(list(NCT_pval, BGGM_pval))
}

library(foreach)
library(doParallel)
library(doRNG)
library(PRROC)
library(pulsar)
library(NetworkComparisonTest)
library(BGGM)
library(parallel)
no_cores <- 55
cl <- makeCluster(no_cores, outfile = "")
clusterEvalQ(cl, {library(COZINE)
  library(PRROC)
  library(pulsar)
  library(parallel)})
registerDoParallel(cl)

reps = 100

fct <- function(i){
  out <- input_function(  data1 = data_null[(i*1000)-1000:i*1000,],
  data2 = data_null_2[(i*1000)-1000:i*1000,])
  return(out)
}



time_parallel <- system.time(
  temp <- foreach(i=1:reps) %dopar% {fct(i)})
stopCluster(cl)
time_parallel

result_NCT_BGGM <- data.frame(NCT_pval, BGGM_pval)
result_NCT_BGGM_name <- paste0(gsub("\\.RData$", "", basename(result1)),"_", gsub("\\.RData$", "", basename(result2)), "_NCT_BGGM.RData")
save(temp, time_parallel, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/diagnostic_plots/results_test_NCT_BGGM/", result_NCT_BGGM_name))
