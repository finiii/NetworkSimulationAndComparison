
rm(list = ls())

#set the parameter values
eta_0 = rep(0.5, 10)
#eta_1 = eta_0

# set the dimension/length of X
dim = 10

#  set the number of samples
n_samples = 100

# Create a sample correlation matrix with negative values
# Generate a random correlation matrix with negative values
theta <- matrix(0, nrow = dim, ncol = dim)

# Fill the upper triangular part with random negative values
for (i in 1:(dim - 1)) {
  for (j in (i + 1):dim) {
    theta[i, j] <- runif(1, -0.1, 0)  # Random value between -1 and 0
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
# die sind independent also kann ich die einzeln generieren, muss nur noch genau schauen wie
# ich bin mir gar nicht sicher ob das stimmt

# write an empty matrix for the samples
X_0 = matrix(NA, nrow = dim, ncol = n_samples)


for (i in 1:dim){
    eta_i = eta_y[i]
    theta_ii =  0 #theta[i,i], we have theta_ii = 0 from the def of the model
    X_0[i,] = rpois(n_samples, exp(eta_i))

}
X <- X_0
#write a loop that does the following code 10 times (ich brauche anderes stopping criterion)
#hier dran arbeite ich gerade!!!

likelihood_values = numeric(1000)
for(loop in 1:1000){
#step 3
 #a) compute the conditional parameters
eta_minus_i = matrix(NA, nrow = dim, ncol = n_samples)
for(j in 1:n_samples){
  for (i in 1:dim){
    #i-th row of theta after removing the i-th column from theta
    theta_minus_i = theta[i,-i]
    eta_minus_i[i,j] = eta_0[i] + theta_minus_i %*% X[-i,j]

}
}

 #b) compute X_1
 X = matrix(NA, nrow = dim, ncol = n_samples)
 #how do I generate from this distribution?
  for(j in 1:n_samples){
    for (i in 1:dim){
      X[i,j] = rpois(1, exp(eta_minus_i[i,j]))
    }
  }


#c) update


# Step 4: caluclate the likelihood
likelihood_samples = numeric(n_samples)
#ich muss das hier mit dem ziel-eta berechnen, also mit eta_0!!!
for(j in 1:n_samples){
  likelihood_samples[j] = exp(t(eta_0) %*% X[,j]+ 0.5 * t(X[,j]) %*% theta %*% X[,j]+ t(rep(1, dim)) %*% -log(factorial(X[,j])) )
}
likelihood = prod(likelihood_samples)
print(likelihood)
likelihood_values[loop] = likelihood
}

#likelihood ändert sich nicht, habe ich einen Fehler im Code????