##############################################
#test the null hypothesis

library(BGGM)

result1 <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_poisson_final/pPGM_simulation_dim10_lower_bound_-0.06_upper_bound_-0.05_n_samples_6_iterations_1e+05_edge_probability_0.3.RData"
load(result1)
data_null <- lapply(names(result), function(name) {
  if (grepl("^X_[0-9]+_simulations$", name)) {
    as.vector(result[[name]][,seq(ncol(result[[name]])/2, ncol(result[[name]]), by = 10) ])
  }
})
data_null <- as.data.frame(do.call(cbind, data_null))
rownames(data_null) <- 1:nrow(data_null)

rm(result)

result2 <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_poisson_final/pPGM_simulation_dim10_lower_bound_-0.08_upper_bound_-0.07_n_samples_6_iterations_1e+05_edge_probability_0.3.RData"
load(result2)
data_null_2 <- lapply(names(result), function(name) {
 if (grepl("^X_[0-9]+_simulations$", name)) {
    as.vector(result[[name]][,seq(ncol(result[[name]])/2, ncol(result[[name]]), by = 10) ])

  }
})
data_null_2 <- as.data.frame(do.call(cbind, data_null_2))
rownames(data_null_2) <- 1:nrow(data_null_2)

rm(result)


input_function <- function(data1, data2)
{

  ### do the BGGM
  test_BGGM <- BGGM::ggm_compare_ppc(data1, data2,
                       iter = 1000)

    BGGM_pval <- test_BGGM[["ppp_jsd"]]
  return(BGGM_pval)
}


reps = 100

fct <- function(i){
 # out <- input_function(  data1 = data_null[(i*1000)-1000:i*1000,],
  #data2 = data_null_2[(i*1000)-1000:i*1000,])
    seed = 2023 + i
    set.seed(seed)
    out <- input_function(  data1 = data_null[sample(nrow(data_null), 1000),],
  data2 = data_null_2[sample(nrow(data_null), 1000),])
  return(out)
}

results <- numeric(reps)
foreach(i=1:reps) {
    results[i] <- fct(i)}


result_BGGM_name <- paste0(gsub("\\.RData$", "", basename(result1)),"_", gsub("\\.RData$", "", basename(result2)), "_NCT_BGGM.RData")

save(results, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/results_test_H0_BGGM/", result_BGGM_name))
