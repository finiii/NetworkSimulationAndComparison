library(coda4microbiome)
library(dplyr)

load(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations/simulation1.RData")
load(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations/simulation2.RData")



data1 <- simulation1$data
data2 <- simulation2$data

#for this we need a bit of a different format

data1 <- as.data.frame(data1) %>% mutate(y = 0)
data2 <- as.data.frame(data2) %>% mutate(y = 1)
data <- rbind(data1, data2) %>% mutate(y = as.factor(y))

#now we can use the coda4microbiome package

#hier bekomme ich aber negative abunances, was ich ja nicht will
#wie kann ich generieren sodass ich keine negativen abundances bekomme?
test1 <- coda4microbiome::coda_glmnet(x = data[,1:10], y = data$y)

###

load("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations/gen_data_e-d_reps_1_zinegbin_samples_500_p_127_cn_5_cluster.RData")
gen_data1 <- gen_data
load("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations/gen_data_e-d_reps_1_zinegbin_samples_500_p_127_cn_100_cluster.RData")
gen_data2 <- gen_data


data1 <- gen_data1[[1]][[1]]
data2 <- gen_data2[[1]][[1]]

data1 <- as.data.frame(data1) %>% mutate(sum = rowSums(data1), y = 0)
data2 <- as.data.frame(data2) %>% mutate(sum = rowSums(data2), y = 1)
data <- rbind(data1, data2) %>% mutate(y = as.factor(y))


test1 <- coda4microbiome::coda_glmnet(x = data[,1:127]/data$sum, y = data$y)
