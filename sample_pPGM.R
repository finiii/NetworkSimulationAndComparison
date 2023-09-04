
rm(list = ls())

#set the parameter values
eta_0 = 1
#eta_1 = eta_0

# set the dimension/length of X
dim = 10

#  set the number of samples
n_samples = 100

# Create a sample correlation matrix with negative values
#das kann ich dann auch noch besser auftomatisch machen, aber gerade ok
theta = matrix(c(
  1.00, -0.70, -0.20, -0.40, -0.50, -0.60, -0.30, -0.10, -0.45, -0.55,
 -0.70,  1.00, -0.60, -0.35, -0.25, -0.15, -0.20, -0.40, -0.30, -0.50,
 -0.20, -0.60,  1.00, -0.50, -0.40, -0.25, -0.15, -0.10, -0.20, -0.30,
 -0.40, -0.35, -0.50,  1.00, -0.45, -0.35, -0.10, -0.20, -0.15, -0.25,
 -0.50, -0.25, -0.40, -0.45,  1.00, -0.55, -0.60, -0.30, -0.20, -0.10,
 -0.60, -0.15, -0.25, -0.35, -0.55,  1.00, -0.45, -0.40, -0.30, -0.20,
 -0.30, -0.20, -0.15, -0.10, -0.60, -0.45,  1.00, -0.35, -0.25, -0.40,
 -0.10, -0.40, -0.10, -0.20, -0.30, -0.40, -0.35,  1.00, -0.55, -0.25,
 -0.45, -0.30, -0.20, -0.15, -0.20, -0.30, -0.25, -0.55,  1.00, -0.45
 -0.45, -0.40, -0.10, -0.35, -0.10, -0.30, -0.20, -0.25, -0.35,  1.00
), nrow = 10, byrow = TRUE)

#set the diagoal to 0
diag(theta) = 0



# Step 1: Compute eta_y
# since we have a fixed eta, we do not need to compute it 
eta_y = eta_0

# Step 2: generate X_0 from the independent model with eta_y
# die sind independent also kann ich die einzeln generieren, muss nur noch genau schauen wie
# ich bin mir gar nicht sicher ob das stimmt

# write an empty matrix for the samples
X_0 = matrix(NA, nrow = dim, ncol = n_samples)


for (i in 1:dim){
    eta_i = eta_y
    theta_ii =  0 #theta[i,i], we have theta_ii = 0 from the def of the model
    #A_ii = wie berechne ich A? A = exp(...)
    A_ii = exp(eta_i)
    X_0[i,] = rpois(n_samples, exp(eta_i))

}

#step 3
 #a) compute the conditional parameters
eta_minus_i = matrix(NA, nrow = dim, ncol = n_samples)
for(j in 1:n_samples){
  for (i in 1:dim){
    #i-th row of theta after removing the i-th column from theta
    theta_minus_i = theta[i,-i]
    eta_minus_i[i,j] = eta_0 + theta_minus_i %*% X_0[-i,j]

}
}

 #b) compute X_1
 X = matrix(NA, nrow = dim, ncol = n_samples)
 #how do I generate from this distribution?

#c) update

X = X_0
# Step 4: caluclate the likelihood
j = 1
likelihood = numeric(n_samples)
for(j in 1:n_samples){
  likelihood[j] = exp(t(eta_minus_i[,j]) %*% X[,j]+ 0.5 * t(X[,j]) %*% theta %*% X[,j]+ t(rep(1, dim)) %*% -log(factorial(X[,j])) - exp(eta_minus_i[,j]))
}

#das A stimmt hier sicher noch nicht!!!!!
