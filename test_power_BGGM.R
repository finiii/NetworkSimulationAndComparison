##############################################
#test the null hypothesis

library(BGGM)
result1 <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_poisson_final/pPGM_simulation_dim20_lower_bound_-0.01_upper_bound_-0.005_n_samples_6_iterations_1e+05_edge_probability_0.5.RData"
load(result1)
data_null <- lapply(names(result), function(name) {
  if (grepl("^X_[0-9]+_simulations$", name)) {
    as.vector(result[[name]][,seq(ncol(result[[name]])/2, ncol(result[[name]]), by = 10) ])
  }
})
data_null <- as.data.frame(do.call(cbind, data_null))
rownames(data_null) <- 1:nrow(data_null)

rm(result)

result2 <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_poisson_final/pPGM_simulation_dim20_lower_bound_-0.01_upper_bound_-0.005_n_samples_6_iterations_1e+05_edge_probability_0.8.RData"
load(result2)
data_null_2 <- lapply(names(result), function(name) {
 if (grepl("^X_[0-9]+_simulations$", name)) {
    as.vector(result[[name]][,seq(ncol(result[[name]])/2, ncol(result[[name]]), by = 10) ])

  }
})
data_null_2 <- as.data.frame(do.call(cbind, data_null_2))
rownames(data_null_2) <- 1:nrow(data_null_2)

rm(result)


data_null <- lapply(names(result), function(name) {
  if (grepl("^X_[0-9]+_simulations$", name)) {
    #as.vector(result[[name]])
    as.vector(result[[name]][,seq(ncol(result[[name]])/2, ncol(result[[name]]), by = 10) ])
  }
})
data_null <- as.data.frame(do.call(cbind, data_null))

input_function <- function(data1, data2)
{

  ### do the BGGM
  test_BGGM <- BGGM::ggm_compare_ppc(data1, data2,
                       iter = 1000)

    BGGM_pval <- test_BGGM[["ppp_jsd"]]
  return(BGGM_pval)
}

library(foreach)
library(doParallel)
library(doRNG)
library(PRROC)
library(pulsar)
library(BGGM)
library(parallel)
no_cores <- 30#55
cl <- makeCluster(no_cores, outfile = "")
clusterEvalQ(cl, {
  library(BGGM)
  library(pulsar)
  library(parallel)})
registerDoParallel(cl)

reps = 100

fct <- function(i){
 # out <- input_function(  data1 = data_null[(i*1000)-1000:i*1000,],
  #data2 = data_null_2[(i*1000)-1000:i*1000,])
    out <- input_function(  data1 = data_null[sample(nrow(data_null), 1000),],
  data2 = data_null_2[sample(nrow(data_null), 1000),])
  return(out)
}



time_parallel <- system.time(
  temp <- foreach(i=1:reps) %dopar% {fct(i)})
stopCluster(cl)
time_parallel

result_BGGM_name <- paste0(gsub("\\.RData$", "", basename(result1)),"_", gsub("\\.RData$", "", basename(result2)), "_NCT_BGGM.RData")

save(temp, time_parallel, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/results_test_H0_BGGM/", result_BGGM_name))
