
##library(coda4microbiome)
library(glmnet)
library(pROC)
library(ggpubr)

result1 <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_poisson_final/pPGM_simulation_dim20_lower_bound_-0.01_upper_bound_-0.005_n_samples_6_iterations_1e+05_edge_probability_0.5.RData"
load(result1)
data_null <- lapply(names(result), function(name) {
  if (grepl("^X_[0-9]+_simulations$", name)) {
    as.vector(result[[name]][,seq(ncol(result[[name]])/2, ncol(result[[name]]), by = 10) ])
  }
})
data_null <- as.data.frame(do.call(cbind, data_null))
rownames(data_null) <- 1:nrow(data_null)

rm(result)

result2 <- "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_poisson_final/pPGM_simulation_dim20_lower_bound_-0.01_upper_bound_-0.005_n_samples_6_iterations_1e+05_edge_probability_0.8.RData"
load(result2)
data_null_2 <- lapply(names(result), function(name) {
 if (grepl("^X_[0-9]+_simulations$", name)) {
    as.vector(result[[name]][,seq(ncol(result[[name]])/2, ncol(result[[name]]), by = 10) ])

  }
})
data_null_2 <- as.data.frame(do.call(cbind, data_null_2))
rownames(data_null_2) <- 1:nrow(data_null_2)

rm(result)



impute_zeros<-function(x){
  if (min(as.numeric(unlist(x)))< 0) {
    stop("Negative abundance values (check your data)")
  } else {
    if (sum(x==0)>0){
      #xmin = min(x[x > 0]);
      xmin = min(as.numeric(unlist(x))[as.numeric(unlist(x))>0])
      if (xmin >= 1){
        x = x + 1;
      }else{
        x = x+xmin/2;
      }
    }
  }
  return(x)
}


logratios_matrix<-function(x){

  x<-impute_zeros(x)

  if(is.null(colnames(x))) colnames(x)<-(1:ncol(x))

  k<-ncol(x)
  m<-nrow(x)

  logx <- log(x)
  lrcolnames<-NULL
  lrX <- matrix(0,m,k*(k-1)/2)
  idlrX <- matrix(0,k*(k-1)/2,2)
  nameslrX <- matrix(0,k*(k-1)/2,2)
  colnamesx <- colnames(x)
  lloc <-0
  for(i in (1:(k-1))){
    for(j in ((i+1):k)) {
      lloc=lloc+1
      idlrX[lloc,]<-c(i,j)
      nameslrX[lloc,] <- c(colnamesx[i],colnamesx[j])
      lrX[,lloc] <- logx[,i]-logx[,j]
      lrcolnames<-c(lrcolnames,paste(paste("lr",i,sep=""),j,sep="."))
    }
  }
  colnames(lrX)<-lrcolnames
  results <- list(
    "matrix of log-ratios" = lrX,
    "pairs of variables in the logratio" = idlrX,
    "names of the variables in the logratio" = nameslrX
  )

  return(results)
}


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

  accuracy<-rep(0,niter)
  accuracy0 <- coda_glmnet0(x=x,lrX=lrX,idlrX=idlrX,nameslrX=nameslrX,y=y1,lambda=lambda0,covar=covar0, alpha=alpha0)$`mean cv-AUC`
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
    "accuracy0"=accuracy0,
    "p-value"=sum(accuracy>=accuracy0)/niter,
    "confidence interval"=quantile(accuracy, c((sig/2),(1-(sig/2))))
  )
  return(results)
}

coda_glmnet0<-function(x,lrX,idlrX,nameslrX,y, covar=NULL, lambda="lambda.1se",alpha=0.9){
  #suppressWarnings()

  # library(glmnet)
  # library(pROC)
  # library(ggpubr)


  if (sum(x==0)>0){
    x<-impute_zeros(x)
  }


  kselect<-ncol(x)

  idlrXsub<-idlrX
  lrXsub<-lrX

  y.binary<-ifelse(dim(table(y))==2, TRUE, FALSE)

  if (y.binary==TRUE){
    if(is.null(covar)){
      lassocv<-glmnet::cv.glmnet(lrXsub,y, family = "binomial" , alpha=alpha, type.measure = "auc", keep=TRUE, parallel = TRUE)
    } else {
      df0<-data.frame(y,covar)
      model0<-glm(y~., family = "binomial", data=df0)
      x0<-predict(model0)
      lassocv<-glmnet::cv.glmnet(lrXsub,y, family = "binomial" , offset=x0,alpha=alpha, type.measure = "auc", keep=TRUE, parallel = TRUE)
    }
  } else {
    if(is.null(covar)){
      lassocv<-glmnet::cv.glmnet(lrXsub,y , alpha=alpha, type.measure = "deviance", keep=TRUE, parallel = TRUE)
    } else {
      df0<-data.frame(y,covar)
      model0<-lm(y~., data=df0)
      x0<-predict(model0)
      lassocv<-glmnet::cv.glmnet(lrXsub,y , offset=x0, alpha=alpha, type.measure = "deviance", keep=TRUE, parallel = TRUE)
    }
  }

  lambdavalue<-lambda
  if (is.character(lambda)){
    if (lambda=="lambda.1se") lambdavalue <-lassocv$lambda.1se
    if (lambda=="lambda.min") lambdavalue <-lassocv$lambda.min
  }
  idrow<-max(which(lassocv$glmnet.fit$lambda>=lambdavalue))  # idrow= row in glmnet.fit object corresponding to the specified lambda

  coeflr<-as.vector(coef(lassocv, s = lambda))[-1]
  lrselect<-which(coeflr!=0)

  idlrXsub[lrselect,]

  coeflogcontrast<-rep(0,ncol(x))
  for (i in (1:length(coeflr))){
    coeflogcontrast[idlrXsub[i,1]]<-coeflogcontrast[idlrXsub[i,1]]+coeflr[i]
    coeflogcontrast[idlrXsub[i,2]]<-coeflogcontrast[idlrXsub[i,2]]-coeflr[i]
  }

  coeflogcontrast<-2*coeflogcontrast/sum(abs(coeflogcontrast))
  varlogcontrast<-which(abs(coeflogcontrast)>0)
  coeflogcontrast<-coeflogcontrast[varlogcontrast]

  (names.select<-colnames(x)[varlogcontrast])

  (positive<-ifelse(coeflogcontrast>0,1,0))

  positive<-factor(positive, levels = c(0,1), labels = c("negative","positive"))

  #df<-data.frame(taxa.name=names.select, taxa.num=varlogcontrast, coefficient=round(coeflogcontrast,digits = 2), positive)

  logcontrast=as.matrix(log(x)[,varlogcontrast])%*%coeflogcontrast
  # logcontrast<-logcontrast-mean(logcontrast)

  if (is.null(covar)){
    predictions<-logcontrast
  } else {
    if (y.binary==TRUE){
      df1<-data.frame(y,logcontrast, covar)
      m1<-glm(y~., family = "binomial", data=df1)
      predictions<-predict(m1)

    } else {
      df1<-data.frame(y,logcontrast, covar)
      m1<-lm(y~., data=df1)
      predictions<-predict(m1)
    }

    # predictions<-predictions-mean(predictions)

  }

  if (y.binary==TRUE){
    AUC_signature<-pROC::auc(pROC::roc(y, as.numeric(predictions),quiet = TRUE))[[1]]
    if (length(varlogcontrast)==0) AUC_signature<- 0.5
    mcvAUC<-lassocv$cvm[idrow]
    sdcvAUC<-lassocv$cvsd[idrow]

  } else {
    mcvMSE<-lassocv$cvm[idrow]
    sdcvMSE<-lassocv$cvsd[idrow]
    Rsq<-as.numeric(cor(predictions,y)^2)
    if (length(varlogcontrast)==0) Rsq <- 0
  }

  if (y.binary==TRUE){
    results <- list(
      "taxa.num" = varlogcontrast,
      "taxa.name" = names.select,
      "log-contrast coefficients" = coeflogcontrast,
      "predictions"=predictions,
      "apparent AUC"= AUC_signature,
      "mean cv-AUC"= mcvAUC,
      "sd cv-AUC"= sdcvAUC)
  } else {
    results <- list(
      "taxa.num" = varlogcontrast,
      "taxa.name" = names.select,
      "log-contrast coefficients" = coeflogcontrast,
      "predictions"=predictions,
      "apparent Rsq" = Rsq,
      "mean cv-MSE"= mcvMSE,
      "sd cv-MSE"= sdcvMSE)
  }
  return(results)
}

input_function <- function(data1, data2)
{
# Set 'y' column to 0 in data1
data1$y <- 0

# Set 'y' column to 1 in data2
data2$y <- 1

# Combine data1 and data2 into a new data frame 'data'
data <- rbind(data1, data2)

# Convert 'y' column to a factor
data$y <- as.factor(data$y)

#do the coda4microbiome test
#test_coda <- coda_glmnet_null_altered(x = data[,1:(ncol(data)-1)], y = data$y, niter = 1000)
p_value_coda <- coda_glmnet_null(x = data[,1:(ncol(data)-1)], y = data$y, niter = 1000)$`p-value`
#return results
return(list(p_value_coda = p_value_coda))
}

library(foreach)
library(doParallel)
library(doRNG)
library(PRROC)
library(pulsar)
library(glmnet)
library(pROC)
library(ggpubr)
#library(coda4microbiome)
library(parallel)
no_cores <- 4
cl <- makeCluster(no_cores, outfile = "log.txt")
clusterEvalQ(cl, {#library(coda4microbiome)
library(glmnet)
library(pROC)
library(ggpubr)
 library(pulsar)
 library(parallel)
 })


registerDoParallel(cl)

reps = 100

fct <- function(i){
  out <- input_function(  data1 = data_null[sample(nrow(data_null), 1000),],
  data2 = data_null_2[sample(nrow(data_null_2), 1000),])
  return(out)
}


time_parallel <- system.time(
  temp <- foreach(i=1:reps) %dopar% {fct(i)})
stopCluster(cl)
time_parallel


#result_coda_name <- paste0(gsub("\\.RData$", "", basename(result1)), "_coda4microbiome.RData")
#save(temp, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/results_test_H0_coda4microbiome/", result_coda_name))

#time_parallel <- system.time(
#  temp <- foreach(i=1:reps) %dopar% {fct(i)})
#stopCluster(cl)
#time_parallel

result_coda_name <- paste0(gsub("\\.RData$", "", basename(result1)),"_", gsub("\\.RData$", "", basename(result2)), "_coda4microbiome.RData")
save(temp, time_parallel, file = paste0("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/results_test_H0_coda4microbiome/", result_coda_name))


