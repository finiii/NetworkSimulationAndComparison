rm(list = ls())
#set the parameter values
eta_0 = rep(3, 10)
#eta_1 = eta_0
# set the dimension/length of X
dim = 10
# set the number of samples
n_samples = 1000
# Create a sample correlation matrix with negative values
# Generate a random correlation matrix with negative values
theta <- matrix(0, nrow = dim, ncol = dim)
# Fill the upper triangular part with random negative values
for (i in 1:(dim - 1)) {
  for (j in (i + 1):dim) {
    theta[i, j] <- -.05#runif(1, -.5, 0) # Random value
    theta[ j,i]=theta[i, j]
  }
}
# Copy the upper triangular values to the lower triangular part
theta[lower.tri(theta)] <- t(theta)[lower.tri(theta)]
# Set the diagonal elements to 0
diag(theta) <- 0
# Step 1: Compute eta_y
# since we have a fixed eta, we do not need to compute it
eta_y = eta_0
# Step 2: generate X_0 from the independent model with eta_y
# write an empty matrix for the samples
X_0 = matrix(NA, nrow = dim, ncol = n_samples)
for (i in 1:dim){
  eta_i = eta_y[i]
  theta_ii = 0 #theta[i,i], we have theta_ii = 0 from the def of the model
  X_0[i,] = rpois(n_samples, exp(eta_i))
}
# iterations
iterations = 100

report = numeric(iterations)

log_likelihood_values = numeric(iterations)
eta_minus_i_old = matrix(eta_0, nrow = dim, ncol = n_samples)
eta_minus_i= matrix(NA, nrow = dim, ncol = n_samples)
X_old = X_0
X_new = matrix(NA, nrow = dim, ncol = n_samples)
mean_X_new = matrix(NA, nrow = dim, ncol = iterations)
for(loop in 1:iterations){
  #step 3
  #a) compute the conditional parameters
  for(j in 1:n_samples){
    for (i in 1:dim){
      #i-th row of theta after removing the i-th column from theta
      theta_minus_i = theta[i,-i]
      eta_minus_i[i,j] = eta_minus_i_old[i,j] + theta_minus_i %*% X_old[-i,j]
      X_new[i,j] = rpois(1, exp(eta_minus_i[i,j]))
      #X_old = X_new
    }
  }
  #update X_old and eta_minus_i_minus_1
  #eta_minus_i_old = eta_minus_i
  
  #c) update
  # Step 4: caluclate the likelihood
  log_likelihood_samples = numeric(n_samples)
  for(j in 1:n_samples){
    log_likelihood_samples[j] = t(eta_0) %*% X_new[,j]+ 0.5 * t(X_new[,j]) %*% theta %*% X_new[,j]+
      t(rep(1, dim)) %*% -log(factorial(X_new[,j]))
  }
  Log_likelihood = sum(log_likelihood_samples)
  #print(likelihood)
  log_likelihood_values[loop] = Log_likelihood
  
  # Step 4: caluclate the mean 
  mean_X_new[,loop] = rowMeans(X_new)
    if (loop>1){
    report[loop] = mean(X_new)
  }

  X_old = X_new
}

pdf(file = "/dss/dsshome1/03/ga27hec2/NetworkSimulationAndComparison/mean_pPGM_later_update.pdf")
#das evtl als ggplot schreiben für alle variablen
plot(mean_X_new[1,], type = "l")
plot(mean_X_new[2,], type = "l")
plot(mean_X_new[3,], type = "l")
plot(mean_X_new[4,], type = "l")
plot(mean_X_new[5,], type = "l")
plot(mean_X_new[6,], type = "l")
plot(report, type = "l")
dev.off()

