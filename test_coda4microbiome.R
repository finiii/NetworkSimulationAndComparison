library(coda4microbiome)
library(dplyr)

load("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations/gen_data_e-d_reps_100__samples_250_p_10_cn_50_band.RData")
gen_data1 <- result$data_1
gen_data2 <- result$data_2

n_reps <- length(gen_data2)
 i = 1

for (i in 1:n_reps) {
#for this we need a bit of a different format

#it seems like the data needs to be transformed to proportions but this also does not make too much sense right??
data1 <- as.data.frame(gen_data1[[i]])
data1 <- data1# %>% mutate_all(.funs = exp) %>% mutate(sum = rowSums(.))%>%
  #mutate(across(everything(), ~ .*sum)) %>% mutate(y = 0)
data2 <- as.data.frame(gen_data2[[i]])# %>% mutate_all(.funs = exp)  %>% mutate(sum = rowSums(.))%>%
  #mutate(across(everything(), ~ .*sum)) %>% mutate(y = 1)
data <- rbind(data1, data2) %>% mutate(y = as.factor(y))

#now we can use the coda4microbiome package

test <- coda4microbiome::coda_glmnet(x = data[,1:20], y = data$y)
#print(test)
summary(test)
}

