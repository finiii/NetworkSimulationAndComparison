
library(SpiecEasi)

### it try to construct a graph not using the SpiecEasi::synth_comm_from_counts functino because I don't have the data

# this whole thing is unweighted, how can I get a weighted version?!?

# Change this number by differences of thousands every time you run the script
#sd = 3828471 # seed - see below


### Choose condition number of precision matrix

condition_number = 10


### Choose number of samples

no_samples = 250


#specify the number of edges
d = 20
e = 80


#Choose graph_type between 'cluster', 'band' or 'scale_free'
#graph_type = 'cluster' 

#einen graph bauen, das mache ich hier

    #set.seed(sd)
    graph <- SpiecEasi::make_graph(D = d, e = e, method = "erdos_renyi")

    graph_plot <- igraph::graph_from_adjacency_matrix(
  graph,
  mode = "undirected"
)

#pdf(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/graphs_SpiecEasi.pdf")
#plot(graph_plot)
#dev.off()

#for an unweighted graph, choose the theta limits to be c(1,1)
    Prec <- SpiecEasi::graph2prec(graph, targetCondition = condition_number)
    Cor <- cov2cor(prec2cov(Prec))

    round(Cor, 3)
    pracma::cond(Cor)


#how is it with the seeds bc when I set one it will always generate the same multivariate normal distr
  gen_data1 <- list()
  gen_data2 <- list()
  #seeds1 <- list()
 #seeds2 <- list()
  no_reps <- 100
  for (i in 1: no_reps){
    print(i)
    #set.seed(sd + i)
    x = SpiecEasi::rmvnorm(n=no_samples, mu=rep(0, d), Sigma=Cor)
    gen_data1[[i]] <- x
    #seeds1[[i]] <- sd + i
  }

#create a second datasets under different seeds!! 
 for (i in 1: no_reps){
    print(i)
    #set.seed(sd + i + 100000)
    y = SpiecEasi::rmvnorm(n=no_samples, mu=rep(0, d), Sigma=Cor)
    gen_data2[[i]] <- y
    #seeds2[[i]] <- sd + i + 100000
  }
  
  
  result <- list(data_1=gen_data1, data_2 = gen_data2, sigma=Cor, graph = graph)

result$sigma
#geht also auch mit kleienerer condition number!!!

pdf(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/graphs.pdf")
plot(graph_plot)
dev.off()

  #hier ändern je nachdem ob weighted oder nicht weighted
save(result,  file = paste("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_SpiecEasi/gen_data_SpiecEasi_weighted_reps_", no_reps, "_" ,
                              "_samples_",no_samples,"_d_",d, "_e_", e, "_cn_",condition_number, "_", ".RData", sep=''))



