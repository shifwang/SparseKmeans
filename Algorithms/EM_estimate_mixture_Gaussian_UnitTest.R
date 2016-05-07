setwd('~/Documents/')
source('EM_estimate_mixture_Gaussian.R')
# Experiment 1
data <- matrix(c(1.1, 1, 2, 2, 3, 3, 2.1, 2.0, 1, 1.1, 3.1, 3), 6, 2)
out <- EstimateMixtureGaussian(data, 3, 0, max.iter = 100)

for (step in out$history) {
  print(step$log.likelihood + step$penalty.value)
}
# Output:
# Calculating initials using kmeans... Passed 0 min..
# iter: 1, log.likelihood: -2.35 Passed 0 min..
# iter: 2, log.likelihood: 10.71 Passed 0 min..
# iter: 3, log.likelihood: 38.29 Passed 0 min..
# iter: 4, log.likelihood: 59.03 Passed 0 min..
# iter: 5, log.likelihood: 59.03 Passed 0 min..
# iter: 6, log.likelihood: 59.03 Passed 0 min..
# iter: 7, log.likelihood: 59.03 Passed 0 min..
# iter: 8, log.likelihood: 59.03 Passed 0 min..
# iter: 9, log.likelihood: 59.03 Passed 0 min..
# iter: 10, log.likelihood: 59.03 Passed 0 min..
# iter: 11, log.likelihood: 59.03 Passed 0 min..
# quit because of absolute value Passed 0 min..

# Experiment 2
give.data <- function (SamplesPerCluster = 20, 
                       Nfeatures = 50, 
                       Nsignals = 10,
                       ClusterNumber = 3, mu = 1) {
  x <- matrix(rnorm(Nfeatures*SamplesPerCluster*ClusterNumber),ncol=Nfeatures)
  for (i in 1:ClusterNumber){
    x[(SamplesPerCluster*(i-1)+1):(SamplesPerCluster*i),1:Nsignals] <- x[(SamplesPerCluster*(i-1)+1):(SamplesPerCluster*i),1:Nsignals] +(i-2)*mu
  }
  x <- scale(x, TRUE, TRUE)
  return(x)
}
data <- give.data(SamplesPerCluster = 20, Nfeatures = 10, 
                  Nsignals = 10, ClusterNumber = 4, mu = 1)
out <- EstimateMixtureGaussian(data, 4, 0, max.iter = 300, 
                               seed = 50, verbose = F)
out$partition
out$BIC
out$noise.feature
# [1] 3 3 3 3 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 3
# [31] 2 2 2 3 2 2 2 2 2 2 4 2 4 4 1 4 4 4 4 4 4 4 4 4 4 2 4 4 4 4
# [61] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

# Experiment 3
data <- give.data(SamplesPerCluster = 20, Nfeatures = 30, 
                  Nsignals = 10, ClusterNumber = 4, mu = 1)
out <- SelectLambda(data, 4, seed = 50, verbose = T, nvals = 20)
out$partition
out$BIC
out$noise.feature
# Result: 
#   out$best.lambda
#   [1] 6.178119
