library(NetCoMi)


result1 <- "NetworkSimulationAndComparison/simulations_poisson_final/pPGM_simulation_dim10_lower_bound_-0.1_upper_bound_-0.05_n_samples_6_iterations_1e+05_edge_probability_0.3.RData"
load(result1)


data_null <- lapply(names(result), function(name) {
  if (grepl("^X_[0-9]+_simulations$", name)) {
    as.vector(result[[name]][,seq(ncol(result[[name]])/2, ncol(result[[name]]), by = 10) ])
  }
})
data_null <- as.data.frame(do.call(cbind, data_null))
rm(result)



input_function <- function(data1, data2)
{

net <- netConstruct(data = data1, 
                           data2 = data2,
                           filtTax = "none",
                           filtSamp = "none",
             dataType = "counts")
netAna <- NetCoMi::netAnalyze(net, graphlet = F, gcmHeat = F)
test <- netCompare(netAna, permTest = TRUE,nPerm = 100L, cores = 1)
rm(net)
rm(netAna)
return(data.frame(test$pvalDiffGlobal))  
rm(test)             
}


library(foreach)
library(doParallel)
library(doRNG)
library(pulsar)
library(NetCoMi)
library(parallel)
no_cores <-  55
cl <- makeCluster(no_cores, outfile = "")
clusterEvalQ(cl, {library(NetCoMi)
  library(PRROC)
  library(pulsar)
  library(parallel)})
registerDoParallel(cl)

reps = 100#100

fct <- function(i){
  #Idee: ich sample einfach 2000 Zeilen jeweils aus den Daten
  out <- input_function(  data1 = data_null[sample(nrow(data_null), 2000),],
  data2 = data_null[sample(nrow(data_null), 2000),])
  return(out)
}



time_parallel <- system.time(
  temp <- foreach(i=1:reps) %dopar% {fct(i)})
stopCluster(cl)
time_parallel

result_NetCoMi_name <- paste0(gsub("\\.RData$", "", basename(result1)), "_NetCoMi.RData") #"_", gsub("\\.RData$", "", basename(result2)),
save(temp, time_parallel, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/results_test_H0_NetCoMi/", result_NetCoMi_name))

