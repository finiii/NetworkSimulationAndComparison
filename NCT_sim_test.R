library(NetworkComparisonTest)

library(foreach)
library(parallel)
library(doParallel)
library(psych)
library("MASS")
#library(GPArotation)
library("corpcor")
library("qgraph")
#library("polycor")
#library(pdist)
library(matrixcalc)
library(ggm)
library(sm)
library(lavaan)
library(Matrix)
library(igraph)

n_reps = 100
######################## END OF CLAUDIAS NCT CODE
############################
graph.data <- function(d, n, gm=NA, theta=NA, prob, v=0.1, u=0.1, upper.bound=.9, lower.bound=.3, n_reps = n_reps){
 # simulates graph
 count=0
 repeat{
 if(any(is.na(gm)) & any(is.na(theta))){
 gm.ig <- erdos.renyi.game(d,prob)
 gm <- as.matrix(get.adjacency(gm.ig))
 theta = matrix(runif(d^2, lower.bound, upper.bound), d, d) * gm
 diag(theta) = 0
 omega = theta * v
 diag(omega) = abs(min(Re(eigen(omega)$values))) + 0.1 + u
 omega <- as.matrix(forceSymmetric(omega))
 if(!is.positive.definite(omega)) gm=theta=NA else break
 }
 if(!any(is.na(gm))& !any(is.na(theta))){
 omega = theta * v
 diag(omega) = abs(min(Re(eigen(omega)$values))) + 0.1 + u
 omega <- as.matrix(forceSymmetric(omega))
 if(is.positive.definite(omega)) break
 }
 if(!any(is.na(gm))& any(is.na(theta))){
 theta = matrix(runif(d^2, lower.bound, upper.bound), d, d) * gm
 omega = theta * v
 diag(omega) = abs(min(Re(eigen(omega)$values))) + 0.1 + u
 omega <- as.matrix(forceSymmetric(omega))
 if(!is.positive.definite(omega)) theta=NA else break
 }
 count=count+1
 print(count)
 }
 sigma = cov2cor(solve(omega))
 omega = solve(sigma)
  gen_data1 <- list()
  gen_data2 <- list()
 for (i in 1: n_reps){
 x1 = mvrnorm(n, rep(0, d), sigma)
 x2 = mvrnorm(n, rep(0, d), sigma)
 gen_data1[[i]] <- x1
 gen_data2[[i]] <- x2
 }

 #sigmahat = cor(x1)
 results <- list(gen_data1=gen_data1, gen_data2 = gen_data2, sigma=sigma, theta=theta, gm=gm)
 return(results)
}


graph1 <- graph.data( d= 20, n = 250, prob = 0.3, n_reps = n_reps)
round(graph1$sigma, 3)

nwinv.pval <- numeric(n_reps)

for (i in 1:n_reps) { #(
data1 <- graph1$gen_data1[[i]]
data2 <- graph1$gen_data2[[i]]

print(paste("currently doing repetition", i))
test <- NCT(data1, data2, gamma = 0.5, weighted = TRUE, it = 1000)

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

sum(result_NCT$nwinv.pval > 0.05)  #H_0 rejected if p-value <= 0.05 
sum(result_NCT$nwinv.pval == 1.00)

#write somthing to save the results
save(result_NCT, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/results/result_NCT_", "NCT_sim_d_10_prob_0.1_n_250_gamma_0.5_it_1000", ".RData"))


