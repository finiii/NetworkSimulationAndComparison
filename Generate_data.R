library(SpiecEasi)

#Petrus code
#ok das jetzt mal so ausprobieren und dann schauen ob das bei coda durchläuft...ansonsten keine ahnung

data(amgut1.filt)
#summary(amgut1.filt)
amgut1.filt <- amgut1.filt[,1:5]
depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total)) 
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

#d <- ncol(amgut1.filt.cs) #number of OTUs

d = 5

input_function <- function(e = d, graph_type, marginal, condition_number = 100, no_samples, prob = 'NULL', sd){
  
    set.seed(sd)
    graph <- SpiecEasi::make_graph(graph_type, d, e)
    Prec <- SpiecEasi::graph2prec(graph, targetCondition = condition_number)
    Cor <- cov2cor(prec2cov(Prec))
    X <- SpiecEasi::synth_comm_from_counts(amgut1.filt.cs, mar=2, distr=marginal, Sigma=Cor, n=no_samples)
    return(list(X, graph))
  
}

#Choose graph_type between 'cluster', 'band' or 'scale_free'
#Choose marginal distribution between 'zinegbin', 'negbin', 'zipois', 'pois', 'lognorm'



#Choose graph type 

graph_type = 'cluster' 

# Change this number by differences of thousands every time you run the script
s = 3728471 # seed - see below

###

### Choose condition number of precision matrix

condition_number = 100

###

### Choose marginal distributions

marginal = 'zinegbin'

### 


### Choose number of samples

no_samples = 500

###

### Choose number of repetitions

no_reps = 1

###
i = 1

gen_data <- list()
seeds <- list()

for (i in 1:no_reps){
  print('track')
  print(i)
  r <- input_function(graph_type = graph_type, marginal = marginal, condition_number = condition_number,
                      no_samples = no_samples, prob = 'NULL', sd = i+s)
  
  seeds[[i]] <- i+s
  gen_data[[i]] <- r
  
}



  
save(gen_data, seeds, file = paste("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations/gen_data_e-d_reps_", no_reps, "_" , marginal,
                              "_samples_",no_samples,"_p_",d,"_cn_",condition_number, "_", graph_type, ".RData", sep=''))



