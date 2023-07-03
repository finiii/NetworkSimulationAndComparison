library(NetworkComparisonTest)

simulation_path <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_SpiecEasi/gen_data_SpiecEasi_weighted_reps_100__samples_250_d_20_e_15_cn_1000_.RData"
load(simulation_path)

gen_data1 <- simulation$data_1
gen_data2 <- simulation$data_2

n_reps <- length(gen_data2)

nwinv.pval <- numeric(n_reps)

#set the gamma
gamma <- 0

for (i in 1:n_reps) { #(
data1 <- gen_data1[[i]]
data2 <- gen_data2[[i]]

print(paste("currently doing repetition", i))
test <- NCT(data1, data2, gamma = gamma)

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
save(result_NCT, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/results/result_NCT_", basename(simulation_path), "gamma", gamma, ".RData"))

result_NCT


#number of rejections of the null hypothesis
sum(result_NCT$nwinv.pval < 0.05)
