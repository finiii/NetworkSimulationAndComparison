##############################################
#test the null hypothesis

library(NetworkComparisonTest)

#result1 <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/pPGM_simulation_dim10_lower_bound_0_upper_bound_0_n_samples_6_iterations_1e+05_edge_probability_0.8.RData"
#load(result1)
#data_null <- lapply(names(result), function(name) {
#  if (grepl("^X_[0-9]+_simulations$", name)) {
#    as.vector(result[[name]])
#  }
#})
#data_null <- as.data.frame(do.call(cbind, data_null))
#rownames(data_null) <- 1:nrow(data_null)

###result2 <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/pPGM_simulation_dim10_lower_bound_0_upper_bound_0_n_samples_6_iterations_1e+05_edge_probability_0.5.RData"
#load(result2)
#data_null_2 <- lapply(names(result), function(name) {
# if (grepl("^X_[0-9]+_simulations$", name)) {
#   as.vector(result[[name]])

# }
#})
#data_null_2 <- as.data.frame(do.call(cbind, data_null_2))
#rownames(data_null_2) <- 1:nrow(data_null_2)

result1 <- #enter path to the data "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_poisson_final/pPGM_simulation_dim10_lower_bound_-0.1_upper_bound_-0.05_n_samples_6_iterations_1e+05_edge_probability_0.3.RData"
load(result1)

data_null <- lapply(names(result), function(name) {
  if (grepl("^X_[0-9]+_simulations$", name)) {
    #as.vector(result[[name]])
    as.vector(result[[name]][,seq(ncol(result[[name]])/2, ncol(result[[name]]), by = 10) ])
  }
})
data_null <- as.data.frame(do.call(cbind, data_null))


input_function <- function(data1, data2)
{
  #do the NCT
  test_NCT <- NCT(data1, data2, it = 1000)
  NCT_pval <- test_NCT$nwinv.pval
  
  return(NCT_pval)
}


library(NetworkComparisonTest)


reps = 100

fct <- function(i){
  # out <- input_function(  data1 = data_null[(i*1000)-1000:i*1000,],
  #data2 = data_null_2[(i*1000)-1000:i*1000,])
  out <- input_function(  data1 = data_null[sample(nrow(data_null), 1000),],
                          data2 = data_null[sample(nrow(data_null), 1000),])
  return(out)
}


results <- numeric(reps)
foreach(i=1:reps) {
  results[i] <- fct(i)}



#result_NCT_name <- paste0(gsub("\\.RData$", "", basename(result1)),"_", gsub("\\.RData$", "", basename(result2)), "_NCT_BGGM.RData")

result_NCT_name <- paste0(gsub("\\.RData$", "", basename(result1)), "_NCT.RData")
save(results, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/results_test_H0_NCT/", result_NCT_name))
