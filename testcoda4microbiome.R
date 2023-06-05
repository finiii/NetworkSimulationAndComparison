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

