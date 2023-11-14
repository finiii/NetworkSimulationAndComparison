load("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/pPGM_simulation_dim10_lower_bound_-0.1_upper_bound_0_n_samples_6_iterations_5000.RData")
result10 <- result
load("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/simulations_Poisson/pPGM_simulation_dim100_lower_bound_-0.1_upper_bound_0_n_samples_6_iterations_20000.RData")
result100 <- result

data10 <- result10$X_1_simulations
data100 <- result100$X_1_simulations
n_samples = 6

pdf("/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/acf_comparison_10_100.pdf")

      for (j in 1:n_samples)  {
        acf(data10[j, -c(1:2500)], main = paste0("X_1"," chain ",j, " dimension 10"))
      }

         for (j in 1:n_samples)  {
        acf(data100[j, -c(1:5000)], main = paste0("X_1", " chain ",j, "dimension 100"))
      }

dev.off()


#