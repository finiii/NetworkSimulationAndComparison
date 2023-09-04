rm(list=ls())

library(coda4microbiome)
library(dplyr)

load("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations/simulation_10_edges_0.3_prob_0.3_0.9_bounds_1000_sample_size_100_repetitions.RData")
gen_data1 <- simulation$data_1
gen_data2 <- simulation$data_2

n_reps <- length(gen_data2)
 i = 1

for (i in 1:n_reps) {
#for this we need a bit of a different format

#it seems like the data needs to be transformed to proportions but this also does not make too much sense right??
data1 <- as.data.frame(gen_data1[[i]])
data1 <- data1 %>% mutate(y = 0)
data2 <- as.data.frame(gen_data2[[i]])  %>% mutate(y = 1)
data <- rbind(data1, data2) %>% mutate(y = as.factor(y))

#now we can use the coda4microbiome package

test <- coda4microbiome::coda_glmnet(x = data[,1:(ncol(data)-1)], y = data$y)
#print(test)
summary(test)
}

