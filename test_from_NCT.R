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


######################## END OF CLAUDIAS NCT CODE
############################
graph.data <- function(d, n, gm=NA, theta=NA, prob, v=0.1, u=0.1, upper.bound=.9, lower.bound=.3){
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
 x1 = mvrnorm(n, rep(0, d), sigma)
    x2 = mvrnorm(n, rep(0, d), sigma)
 #sigmahat = cor(x1)
 results <- list(data1=x1, data2 = x2, sigma=sigma, theta=theta, gm=gm)
 return(results)
}

graph1 <- graph.data( d= 20, n = 5, prob = 0.2)
NCT(graph1$data1, graph1$data2, test.edges = T, gamma = 0.5)

############ programming the permutation test in parallel test
############################
nCores <- 16
loops<-10 #1000
#conditions
#hypotheses: H0 =H0 H1 =maximum edge split halves, H2=maximum edge is set to zero
eq<-"equal" #group1 is denser than group2, equal= equal groupsizes, unequal1 = group1 (the denser group) is larger, unequal2= group2 the less dense group
#is larger
dens<-0.1
l<-4#samplesizes
threeGammas <- list(gamma0=list(), gamma25=list(),gamma50=list())
 for(gam in 1:3){
 if(gam==1){gamma<-0}else if(gam==2) {gamma<-0.25} else {gamma<-0.5}
threeVariableNumbers < -list(NV10=list(), NV20=list(), NV30=list())
for(nv in 1:3){
 nV<-nv*10
 threehypotheses <- list(H0=list(), H1=list(), H2=list()) #H0 is equal, H1
#is edge split half H2 is edge zero
 for(h in 0:2){
 hypothesis <- h
 allsamplesizes <- list(nobs250=list(), nobs500=list(), nobs750=list(),
nobs1000=list())
 for(n in 1:l){
 if(eq=="equal"){nobs <- n*250; nobs2<-n*250} else if(eq=="unequal1")
{nobs<-n*250; nobs2<-1.5*nobs} else {nobs2<-n*250; nobs<-1.5*nobs2}
 minnobs <- n*250
 cl <- nCores-1 # no of cores
 cl <- makeCluster(cl, outfile="") #registering
 registerDoParallel(cl)
 #saveoutput
 corgr1<-list(I(replicate(loops,matrix(rep(NA, nV*nV), nV))))
 corgr2<-list(I(replicate(loops,matrix(rep(NA, nV*nV), nV))))
 NW1<-list(I(replicate(loops,matrix(rep(NA, nV*nV), nV))))
 NW2<-list(I(replicate(loops,matrix(rep(NA, nV*nV), nV))))
 pvalues<-rep(NA,loops)
 #loop
 output<-foreach(gr=1:loops, .combine='rbind', .multicombine=TRUE,
.init=list(list(), list()), .packages=c("psych", "MASS", "GPArotation",
"corpcor", "qgraph", "polycor", "pdist", "matrixcalc", "ggm", "sm", "lavaan",
"igraph", "Matrix", "matrixcalc", "MASS")) %dopar% {
 graph<-graph.data(nV, nobs, prob=dens, upper.bound=-1,
lower.bound=-5)
 cor1<-graph$sigma
 pcor1<-round(cor2pcor(cor1),4)
 data1<-graph$data #mvrnorm(nobs,mu=rep(0,nV), Sigma=cor1)
 lowerpcor<-pcor1[lower.tri(pcor1)]
 edges<-which(lowerpcor!=0)
 nedges<-length(edges)
 if(hypothesis==0){
 pcor2<-pcor1
 data2<-mvrnorm(nobs2,mu=rep(0,nV),Sigma=cor1)
 }
 if(hypothesis==1){
 row<-as.numeric(which(pcor1 == max(lowerpcor), arr.ind = T)[1,1])
 col<-as.numeric(which(pcor1 == max(lowerpcor), arr.ind = T)[1,2])
 newvalue<-max(lowerpcor)/2
 pcor2<-pcor1
 pcor2[row,col]<-newvalue
 pcor2[col,row]<-newvalue
 cor2<-pcor2cor(pcor2)
 data2<-mvrnorm(nobs2, mu=rep(0,nV), Sigma=cor2)
 }
 if(hypothesis==2){
 row<-as.numeric(which(pcor1 == max(lowerpcor), arr.ind = T)[1,1])
 col<-as.numeric(which(pcor1 == max(lowerpcor), arr.ind = T)[1,2])
 newvalue<-0
 pcor2<-pcor1
 pcor2[row,col]<-newvalue
 pcor2[col,row]<-newvalue
 cor2<-pcor2cor(pcor2)
 data2<-mvrnorm(nobs2, mu=rep(0, nV), Sigma=cor2)
 }
 pvalues[gr]<-NCT(data1, data2, gamma=gamma,
progressbar=F)$nwinv.pval
 corgr1[[gr]]<-cor(data1)
 corgr2[[gr]]<-cor(data2)
 NW1[[gr]]<-pcor1
 NW2[[gr]]<-pcor2
 rejections<-list("pvalues"=pvalues[gr], "corgr1"=I(corgr1[[gr]]),
"corgr2"=I(corgr2[[gr]]), "NW1"=I(NW1[[gr]]), "NW2"=I(NW2[[gr]]))
 return(rejections)
 }
 stopCluster(cl)
 #str(output)
 output<-output[-1,]
 rownames(output)<-NULL

 #separate vectors from matrices (output1, elements are values,
#output2 elements are matrices)
 output1<-output[,1]
 output2<-output[,2:5]

 #when one column in output1:
 dataA<-unlist(output1)
 #when multiple columns in output1:
 #dataA<-unlist(output1[,1])
 #for(i in 2:ncol(output1)){
 # dataA<-cbind(dataA,unlist(output1[,i]))
 #}


 dataB<-I(output2[,1])
 for(i in 2:ncol(output2)){
 dataB<-data.frame(dataB,I(output2[,i]))
 }

 combinedoutput<-data.frame(dataA,dataB)
 colnames(combinedoutput)<-colnames(output)

 allsamplesizes[[n]]<-combinedoutput


 #sum(pvalues<0.05)/loops
 }
 threehypotheses[[h+1]]<-allsamplesizes
 }
 threeVariableNumbers[[nv]]<-threehypotheses
 }
threeGammas[[gam]]<-threeVariableNumbers
}
saveRDS(threeGammas,file=paste0("output_invariance_dens1_equal.rds")))