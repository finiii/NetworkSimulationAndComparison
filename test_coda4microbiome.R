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
coda_glmnet_null_altered<-function(x,y,niter=100,covar=NULL,lambda="lambda.1se", alpha=0.9,sig=0.05){

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
    "confidence interval"=quantile(accuracy, c((sig/2),(1-(sig/2)))),
    "p.value" = sum(accuracy >= accuracy0)/(niter),
    "accuracy0" = accuracy0
  )
  return(results)
}

#need to add this internal function from their github repo
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
      lassocv<-glmnet::cv.glmnet(lrXsub,y, family = "binomial" , alpha=alpha, type.measure = "auc", keep=TRUE)
    } else {
      df0<-data.frame(y,covar)
      model0<-glm(y~., family = "binomial", data=df0)
      x0<-predict(model0)
      lassocv<-glmnet::cv.glmnet(lrXsub,y, family = "binomial" , offset=x0,alpha=alpha, type.measure = "auc", keep=TRUE)
    }
  } else {
    if(is.null(covar)){
      lassocv<-glmnet::cv.glmnet(lrXsub,y , alpha=alpha, type.measure = "deviance", keep=TRUE)
    } else {
      df0<-data.frame(y,covar)
      model0<-lm(y~., data=df0)
      x0<-predict(model0)
      lassocv<-glmnet::cv.glmnet(lrXsub,y , offset=x0, alpha=alpha, type.measure = "deviance", keep=TRUE)
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

test_coda <- coda_glmnet_null_altered(x = data[,1:(ncol(data)-1)], y = data$y, niter = 1000)
test_coda$p.value
test_coda


