library(NetworkComparisonTest)


load(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations/simulation1.RData")
load(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations/simulation2.RData")

data1 <- simulation1$data
data2 <- simulation2$data

test1 <- NCT(data1, data2, test.edges = T, test.centrality = TRUE)

#The p value resulting from the permutation test concerning the maximum difference in edge weights.
test1$nwinv.pval
#p-values (corrected for multiple testing or not according to ’p.adjust.methods’)
#per edge from the permutation test concerning differences in edges weights
test1$einv.pvals
#p-values(corrected for multiple testing or not according to ’p.adjust.methods’)
#per node from the permutation test concerning differences in centralities
test1$diffcen.pval

test2 <- NCT(data1, data1, test.edges = T, test.centrality = TRUE)

#The p value resulting from the permutation test concerning the maximum difference in edge weights.
test2$nwinv.pval
#p-values (corrected for multiple testing or not according to ’p.adjust.methods’)
#per edge from the permutation test concerning differences in edges weights
test2$einv.pvals
#p-values(corrected for multiple testing or not according to ’p.adjust.methods’)
#per node from the permutation test concerning differences in centralities
test2$diffcen.pval