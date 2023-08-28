rm(list = ls())

library(NetworkComparisonTest)

simulation_path <-"/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulation_change_edge/simulation_changed_biggest_edge_10_edges_0.3_prob_0.3_0.4_bounds_250_sample_size_100_repetitions.RData"
load(simulation_path)

pracma::cond(simulation$sigma)

round(simulation$sigma_1,3)
round(simulation$sigma_2,3)
simulation$adj_matrix

gen_data1 <- simulation$data_1
gen_data2 <- simulation$data_2

n_reps <- length(gen_data2)

nwinv.pval <- numeric(n_reps)

#set the gamma
gamma <- 0.5

#set the number of iterations
iterations = 100

for (i in 1:n_reps) { #(
data1 <- gen_data1[[i]]
data2 <- gen_data2[[i]]

print(paste("currently doing repetition", i))
test <- NCT(data1, data2, gamma = gamma, it = iterations)

#extract the p-values from the test

#The p value resulting from the permutation test concerning the maximum difference in edge weights.
nwinv.pval[i] <- test$nwinv.pval

#p-values (corrected for multiple testing or not according to ’p.adjust.methods’)
#per edge from the permutation test concerning differences in edges weights
#einv.pvals[i] <- test$einv.pvals

#p-values(corrected for multiple testing or not according to ’p.adjust.methods’)
#per node from the permutation test concerning differences in centralities
#diffcen.pval[i] <- test$diffcen.pval
}



result_NCT <- data.frame(nwinv.pval)# einv.pvals, diffcen.pval)

#write somthing to save the results
save(result_NCT, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/results/result_NCT_", basename(simulation_path), "gamma", gamma, "_iterations_", iterations, ".RData"))

result_NCT


#number of rejections of the null hypothesis
sum(result_NCT$nwinv.pval < 0.05)


simulation$weighted.adj_1-simulation$weighted.adj_2
