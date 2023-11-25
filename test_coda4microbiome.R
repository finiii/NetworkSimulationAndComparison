library(coda4microbiome)
library(dplyr)

#load my own data to test if it works with them
result <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/pPGM_simulation_dim10_lower_bound_-0.01_upper_bound_-0.005_n_samples_6_iterations_10000.RData"
load(result)
data1 <- lapply(names(result), function(name) {
  if (grepl("^X_[0-9]+_simulations$", name)) {
    as.vector(result[[name]][,seq(ncol(result[[name]])/2, ncol(result[[name]]), by = 35) ])
  }
})
data1 <- as.data.frame(do.call(cbind, data1))
rownames(data1) <- 1:nrow(data1)

#I need data dor the other group
result <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/pPGM_simulation_dim10_lower_bound_-0.1_upper_bound_-0.05_n_samples_6_iterations_50000.RData"
load(result)
data2 <- lapply(names(result), function(name) {
  if (grepl("^X_[0-9]+_simulations$", name)) {
    as.vector(result[[name]][,seq(ncol(result[[name]])/2, ncol(result[[name]]), by = 35) ])
  }
})
data2 <- as.data.frame(do.call(cbind, data2))


#bring them on the same length
data2 <- data2[(nrow(data2)-nrow(data1)):nrow(data2),]
rownames(data2) <- 1:nrow(data2)


#for this we need a bit of a different format

#it seems like the data needs to be transformed to proportions but this also does not make too much sense right??
data1 <- data1 %>% mutate(y = 0)
data2 <- data2  %>% mutate(y = 1)
data <- rbind(data1, data2) %>% mutate(y = as.factor(y))

#now we can use the coda4microbiome package


#change the permutation test slightly such that I can get a p-value
coda_glmnet_null<-function(x,y,niter=100,covar=NULL,lambda="lambda.1se", alpha=0.9,sig=0.05){

  alpha0<-alpha
  lambda0<-lambda
  covar0<-covar


  y.binary<-ifelse(dim(table(y))==2, TRUE, FALSE)
  y1<-y
  lrmatrix<-logratios_matrix(x)

  lrX<-lrmatrix[[1]]
  idlrX<-lrmatrix[[2]]
  nameslrX<-lrmatrix[[3]]

  #my change: also calculate the accuracy with the original data
  lr<-coda_glmnet0(x=x,lrX=lrX,idlrX=idlrX,nameslrX=nameslrX,y=y1,lambda=lambda0,covar=covar0, alpha=alpha0)
  if (y.binary==TRUE){
      res<-lr$`mean cv-AUC`
    } else {
      res<-lr$`mean cv-MSE`
    }
  accuracy0 <- res
  accuracy<-rep(0,niter)
  for(i in (1:niter)){
    y1<-sample(y1)
    lr<-coda_glmnet0(x=x,lrX=lrX,idlrX=idlrX,nameslrX=nameslrX,y=y1,lambda=lambda0,covar=covar0, alpha=alpha0)
    if (y.binary==TRUE){
      res<-lr$`mean cv-AUC`
    } else {
      res<-lr$`mean cv-MSE`
    }
    accuracy[i]<-res
    print(paste("iter",i))
  }
  results <- list(
    "accuracy"=accuracy,
    "confidence interval"=quantile(accuracy, c((sig/2),(1-(sig/2))))
    "p.value" = sum(accuracy >= accuracy0)/(niter)
  )
  return(results)
}

test_coda <- coda_glmnet_null(x = data[,1:(ncol(data)-1)], y = data$y)
test_coda$p.value


