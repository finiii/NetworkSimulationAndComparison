rm(list = ls())

library(NetworkComparisonTest)
library(BGGM)

# same graph
#get list of all files in subfolder simulations
simulations_path <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations"

# Use list.files() to get a list of all files in the folder
simulations <- list.files(simulations_path)


# prepare the result file
result = data.frame(NCT = numeric(0), BGGM = numeric(0), simulation_file = character(0))

for (simulation_name in simulations)
{
    load(paste0(simulations_path, "/", simulation_name))
    # get the data
    gen_data1 <- simulation$data_1
    gen_data2 <- simulation$data_2

    n_reps <- length(gen_data2)

    # prepare the solution vectors
    NCT.pval <- numeric(n_reps)
    BGGM.pval <- numeric(n_reps)

#### Set parameters 
#set the gamma parameter for the NCT
gamma <- 0.5
#set the number of iterations for the NCT and for the BGGM
iterations = 100


for (i in 1:n_reps) { #
data1 <- gen_data1[[i]]
data2 <- gen_data2[[i]]

print(paste("currently doing repetition", i))

### do the NCT
test_NCT <- NCT(data1, data2, gamma = gamma, it = iterations)

#extract the p-values from the test

#The p value resulting from the permutation test concerning the maximum difference in edge weights.
NCT.pval[i] <- test_NCT$nwinv.pval

### do the BGGM
test_BGGM <- BGGM::ggm_compare_ppc(data1, data2,
                       iter = iterations)

BGGM.pval[i] <- test_BGGM[["ppp_jsd"]]
}


#save the results
    result_this_run = data.frame(NCT = NCT.pval, BGGM = BGGM.pval, simulation_file = simulation_name)
    result = rbind(result, result_this_run)
    print(result)

}

save(result, file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/results/result_NCT_BGGM_iter_100.RData")
#load("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/results/result_NCT_BGGM_iter_100.RData")
library(dplyr)
result %>% mutate(NCT_is_lower_than_alpha = (NCT < 0.05), BGGM_is_lower_than_alpha = (BGGM < 0.05)) %>% group_by(simulation_file) %>% summarise_at(vars(NCT_is_lower_than_alpha, BGGM_is_lower_than_alpha), sum)


unique(result$simulation_file)
