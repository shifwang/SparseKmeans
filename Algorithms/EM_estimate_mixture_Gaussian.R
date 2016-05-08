EstimateMixtureGaussian <- function(data, k, lambda,
                                    penalty.type = 'l1',
                                    rel.tol  = 1e-6,
                                    abs.tol  = 1e-7,
                                    max.iter = 300,
                                    verbose  = TRUE,
                                    record   = FALSE,
                                    seed     = 101){
  # use EM algorithm to solve penalized mixture Gaussian likelihood
  # Args:
  #   data         - n times p matrix, observation matrix
  #   k            - number of components of mixture Guassian
  #   lambda       - penalized parameter
  #   penalty.type - type of penalized function
  #   rel.tol      - relative tolerance DEFAULT: 1e-6 
  #   abs.tol      - absolute tolerance DEFAULT: 0.1 * RELTOL
  #   max.iter     - maximum number of iterations
  #   record       - boolean, whether record all the variables
  #   verbose      - boolean, whether print the process
  #
  # Returns:
  #   portion       - length K vector, mixing portion
  #   centers       - K times p matrix, the center of each Gaussian
  #   variance      - length p vector, variance of each feature
  #   noise.feature - length p boolean vector, 0 means relevant, 1 means noise
  #   partition     - length n vector, indicating the cluster
  #   history       - list, store all the intermediate variables, empty if not record
  
  set.seed(seed)  #  set seed
  # Sanity check and prepare
  n <- dim(data)[1]
  p <- dim(data)[2]
  stopifnot(k < n) #  number of components should be smaller than number of samples
  start.time <- proc.time()
  len <- 0
  history <- list()
  
  # Main part
  if (verbose) {
    len <- Print(len, 'Calculating initials using kmeans...', start.time)
  }
  
  out <- kmeans(data, k, nstart = 10)
  partition <- out$cluster
  tau <- matrix(1/k, n, k) #  TODO(Yu): give a better initial 
  centers <- out$centers
  portion <- out$size/n
  noise.feature <- rep(FALSE, p)
  variance <- diag(cov(data)) #  TODO(Yu): make sure this is correct
  
  # Intermediate variables
  log.likelihood <- list() 
  penalty.value <- list() # 
  total.obj <- list() # log.likelihood + penalty.value
  for (iter in 1:max.iter) {
    # E step: update tau
    for (i in 1:n) {
      pr <- exp(-Dist(matrix(data[i,],nrow = 1), centers, 0.5/variance))
      tmp <- portion * pr
      tau[i,] <- tmp/sum(tmp)
    }
    # M step: update other variables
    # update portion
    portion <- apply(tau, 2, sum) / n 
    # update variance
    variance <- matrix(0, p)
    for (i in 1:n) {
      for (j in 1:k) {
        variance <- variance + tau[i,j] * (data[i,] - centers[j,])^2
      }
    }
    variance <- variance / n
    stopifnot(!any(abs(variance) < 1e-10))
    # update centers
    for (j in 1:k) {
      center.hat <- apply(tau[,j] * data, 2, sum) / portion[j] / n
      lambda.hat <- lambda / (n * portion[j]) * variance
      if (penalty.type == 'l1') {
        centers[j,] <- sign(center.hat) * sapply(abs(center.hat) - lambda.hat, function(x) max(x, 0))
      } else {
        stop('Unknown penalty type in main.')
      }
    }
    # record and print
    log.likelihood[[iter]] <- ComputeLogLikelihood(data, portion, centers, variance, tau)
    penalty.value[[iter]] <- ComputePenaltyValue(centers, penalty.type, lambda)
    total.obj[[iter]] <-  log.likelihood[[iter]] + penalty.value[[iter]]
    if (record) {
      history[[iter]] <- list(log.likelihood = log.likelihood[[iter]],
                              penalty.value = penalty.value[[iter]],
                              portion = portion, 
                              centers = centers, 
                              variance = variance, tau = tau)
    }
    if (verbose) {
      len <- Print(len, sprintf('iter: %d, log.likelihood: %3.2f',iter, total.obj[[iter]]), start.time)
    }
    # quit condition check
    if (iter > 10 && abs(total.obj[[iter - 3]] - total.obj[[iter]]) < abs.tol) {
      flag <- 'quit because of absolute value'
      break
    } else if (iter == max.iter) {
      flag <- 'quit because of max iteration'
      break
    }
  }
  if (verbose) {
    len <- Print(len, flag, start.time)
  }
  # calculate noise.feature and partition
  noise.feature <- apply(centers, 2, function(x) sum(abs(x)) == 0)
  partition <- apply(tau, 1, which.max)
  BIC <- -2 * log.likelihood[[iter]] + log(n) * (k * p + k + p - 1 - sum(noise.feature)) #  TODO(Yu): BIC does not make much sense here
  return(list(portion = portion, 
              centers = centers, 
              variance = variance, 
              noise.feature = noise.feature, 
              partition = partition, 
              history = history,
              log.likelihood = log.likelihood[[iter]],
              seed = seed,
              iter = iter,
              BIC  = BIC))
}


Print <- function(len = 0, msg, start.time){
  # customized print info on the screen
  # Args:
  #   len  - the length of previous message, default is 0
  #   msg  - msg to convey
  #   start.time - the start time 
  #
  # Returns:
  #   len - the length of current message
  kReWrite <- FALSE
  if (kReWrite) {
    cat(paste(rep('\b', len), collapse = ''))
  } else {
    #  Do nothing
  }
  str <- paste(msg, sprintf('Passed %1.0f min..\n', (proc.time()[[3]] - start.time[[3]])/60))
  cat(str)
  return(nchar(str))
}

ComputeLogLikelihood <- function(data, portion, centers,
                                 variance, tau) {
  # Compute the log likelihood function for data under mixture Gaussian model
  p <- dim(data)[2]
  likelihood.first <- sum(apply(tau, 2, sum) * log(portion))
  square.dist <- Dist(data, centers, 0.5/variance)
  likelihood.second <- sum(tau * -square.dist) - sum(tau) * (p / 2 * log(pi) + sum(log(variance)))
  return(likelihood.first + likelihood.second)
}
Dist <- function(x, y, weight) {
  # Calculate the square distance between rows of x and y
  stopifnot(dim(x)[2] == dim(y)[2])
  nx <- dim(x)[1]
  ny <- dim(y)[1]
  distance <- matrix(0, dim(x)[1], dim(y)[1])
  for (i in 1:nx){
    for (j in 1:ny){
      distance[i, j] <- sum((x[i,] - y[j,])^2 * weight)
    }
  }
  return(distance)
}
ComputePenaltyValue <- function(centers, penalty.type, lambda) {
  if (penalty.type == 'l1') {
    return(lambda * sum(abs(centers)))
  } else {
    stop('Do not understand the penalty type in ComputePenaltyValue')
  }
}

SelectLambda <- function (x, k = NULL, nperms = 5, 
                          lambda.list = NULL, 
                          nvals = 10, 
                          seed = 101,
                          verbose = TRUE) 
{
  # Use Gap Statistics to select lambda for EM
  #
  # Args:
  #
  #   x  - data matrix
  #   k  - number of clusters
  #   nperms - number of permutations
  #   lambda.list - vector, lambdas to be selected
  #   nvals - number of lambdas if selected by DEFAULT
  #   seed - random generator seed
  #   verbose - whether print out details
  #
  # Returns:
  #   
  #  log.likelihood - nvals vector, log likelihood for each lambda
  #  perm.prob -  nperms * nvals matrix, log likelihood for each lambda and permuted data
  #  lambda.list - nvals vector, candidate lambdas
  #  gaps - nvals vector, gap statistics 
  #  sdgaps - nvals vector, standard variation of gap statistics
  #  best.lambda - selected lambda whose gap statistics is biggest
  if (is.null(lambda.list)) 
    lambda.list <- exp(seq(log(1.2), log(ncol(x) * 0.9), 
                       len = nvals))
  # Sanity check
  stopifnot(min(lambda.list) > 1) 
  stopifnot(length(lambda.list) >= 2) 
  stopifnot(!is.null(k)) 
  permx <- list()
  nnonzerows <- NULL
  set.seed(seed)
  len <- 0
  start.time <- proc.time()
  if (verbose) {
    len = Print(len, 'Generate permutated data.', start.time)
  }
  for (i in 1:nperms) {
    permx[[i]] <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    for (j in 1:ncol(x)) 
      permx[[i]][, j] <- sample(x[, j])
  }
  
  log.likelihood <- c()
  for (i in 1:length(lambda.list)) {
    if (verbose) {
      len <- Print(len, sprintf('Calculate log likelihood for %d-th lambda %3.2f ...', i, lambda.list[i]), start.time)
    }
    out <- EstimateMixtureGaussian(x, k, lambda = lambda.list[i], 
                                   record = FALSE, verbose = F)
    log.likelihood <- c(log.likelihood, out$log.likelihood)
  }
  perm.prob <- matrix(NA, nrow = length(lambda.list), ncol = nperms)
  for (j in 1:nperms) {
    if (verbose) 
      len <- Print(len, sprintf("Calculate Permutation %d out of %d", j, nperms), start.time)
    for (i in 1:length(lambda.list)) {
      out <- EstimateMixtureGaussian(permx[[j]], k, lambda = lambda.list[[i]], 
                                           record = FALSE, verbose = FALSE)
      perm.prob[i, j] <- out$log.likelihood
    }
  }
  gaps <- (log.likelihood - apply(perm.prob, 1, mean))
  if (verbose){
    len <- Print(len, 'Plotting Gap statistics..', start.time)
    plot(lambda.list, gaps, xlab='lambda', ylab='Gap Statistics', title('Gap statistics plot'))
  }
  out <- list(log.likelihood = log.likelihood,
              perm.prob = perm.prob, 
              lambda.list = lambda.list, 
              gaps = gaps, 
              sdgaps = apply(perm.prob, 1, sd), 
              best.lambda = lambda.list[which.max(gaps)])
  return(out)
}