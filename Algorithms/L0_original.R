# give L0 results
# parameter selected using original Gap Statistics
# Use parallel computing to speed up this algorithm
give.cluster <- function (x, K = NULL, wbounds = NULL, 
                          initial.iter = 20, maxiter = 6,
                          seed = 123, verbose = TRUE) 
{
  set.seed(seed)
  wbounds <- c(wbounds)
  out<- NULL
  Cs0 <- kmeans(x, centers = K, nstart = initial.iter)$cluster
  for (i in 1:length(wbounds)) {
    if (verbose) {
      cat('hard sparsity ', wbounds[i], ':')
    }
    Cs <- Cs0
    w <- rep(1/sqrt(ncol(x)), ncol(x))
    w.old <- rep(1, ncol(x)) 
    iter <- 0
    while (sum(w!=w.old) >10  && 
             iter < maxiter) {
      if (verbose) {
        cat(iter,' ')
      }
      iter <- iter + 1
      w.old <- w 
      w <- update.w(x, Cs, wbounds[i])
      Cs <- update.Cs(x, K, w, Cs)
    }
    if (verbose) {
      cat('\n')
    }
    out[[i]] <- list(weight = w, Cs = Cs, wbound = wbounds[i])
  }
  return(out);
}
L1<-function(x){
  return(sum(abs(x)))
}
update.Cs <- function (x, K, ws, Cs, initial.iter = 1) 
{
  if (sum(ws != 0) == 1) {
    only.one.feature <- T
  } else {
    only.one.feature <- F
  }
  if (only.one.feature) {
    z <- matrix(x[, ws != 0],ncol = 1)
  } else {
    z <- x[, ws != 0]
  }
  nrowz <- nrow(z) # number of samples
  mus <- NULL
  if (!is.null(Cs)) {
    for (k in unique(Cs)) {
      if (sum(Cs == k) > 1) {
        if (only.one.feature){
          mus <- rbind(mus, apply(matrix(z[Cs == k, ],ncol = 1), 2, mean))
        } else {
          mus <- rbind(mus, apply(z[Cs == k, ], 2, mean))
        }
      }
      if (sum(Cs == k) == 1) 
        mus <- rbind(mus, z[Cs == k, ])
    }
  }
  if (is.null(mus)) {
    km <- kmeans(z, centers = K, nstart = initial.iter)
  }
  else {
    distmat <- as.matrix(dist(rbind(z, mus)))[1:nrowz, (nrowz +1):(nrowz + K)]
    nearest <- apply(distmat, 1, which.min)
    if (length(unique(nearest)) == K) {
      km <- kmeans(z, centers = mus)
    }
    else {
      km <- kmeans(z, centers = K, nstart = initial.iter)
    }
  }
  return(km$cluster)
}
update.w <- function (x, Cs, wbound) 
{
  bcss.perfeature <- give.bcss(x, Cs)
  weight<-find.w(bcss.perfeature,wbound)
  return(weight)
}
give.bcss <- function (x, Cs,w=NULL) 
{
  wcss.perfeature <- numeric(ncol(x))
  for (k in unique(Cs)) {
    whichers <- (Cs == k)
    if (sum(whichers) > 1) 
      wcss.perfeature <- wcss.perfeature + apply(scale(x[whichers,], center = TRUE, scale = FALSE)^2, 2, sum)
  }
  bcss.perfeature <- apply(scale(x, center = TRUE, scale = FALSE)^2, 
                           2, sum) - wcss.perfeature
  #############################
  #######          ############
  #############################
  if (is.null(w))
    return(bcss.perfeature)
  else
    return(sum(bcss.perfeature*w))
}

select.bound <- function (x, K = NULL, nperms = 5, wbounds = NULL, 
                          nvals = 10, seed = 101, verbose = TRUE) 
{
  require('foreach')
  require('doParallel')
  ncores <- 3
  registerDoParallel(ncores)
  if (is.null(wbounds)) 
    wbounds <-exp(seq(log(10),log(ncol(x)/2),length.out=nvals))
  x.null <- list()#  the samples from x without cluster
  signal.feature <- list() #  a list of relevant features
  tmp <- foreach(iter = 0:nperms) %dopar% { # calculate the null samples "nperm" times
    if (iter == 0) { #  When iter == 0, no permutation performed.
      if (verbose) {
        cat('select.bound(ell_0): perform ell_0_kmeans on original data.\n')
      }
      x.cluster <- give.cluster(x, K, wbounds = wbounds, initial.iter = 1)#  carry out ell_0_kmeans
      bcss <- c()
      for (i in 1:length(x.cluster)) {#  scores of parameters for original data
        signal.feature[[i]] <- x.cluster[[i]]$weight
        bcss <- c(bcss, give.bcss(x, x.cluster[[i]]$Cs, x.cluster[[i]]$weight))
      } 
      return(list(iter, bcss, signal.feature))
    }
    if (verbose) {
        cat('select.bound(ell_0): perform ell_0_kmeans on iter ', iter, '.\n')
    }
    set.seed(seed * iter)  #  set random seed
    x.null <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    for (j in 1:ncol(x)) x.null[, j] <- sample(x[, j]) #  generate data
    perm.out <- give.cluster(x.null, K, wbounds = wbounds, initial.iter = 1) #  calculate scores for permuted data
    null.bcss <- rep(NA, length(wbounds))
    for (i in 1:length(perm.out)) {
      null.bcss[i] <- give.bcss(x.null, perm.out[[i]]$Cs,perm.out[[i]]$weight)
    }
    return(list(iter, null.bcss))
  }
  gaps <- rep(0, length(wbounds))
  for (item in tmp) {
    if (item[[1]] == 0) {
      gaps <- gaps + log(item[[2]])
      signal.feature <- item[[3]]
    } else {
      gaps <- gaps - log(item[[2]])/nperms
    }
  }
  out <- list(signal.feature = signal.feature[[which.max(gaps)]], gaps = gaps, wbounds = wbounds, 
              bestw = wbounds[which.max(gaps)])
  if (verbose) {
    plot(wbounds, gaps, xlab = 'Number of Non-zero Weights', ylab = 'Gap Statistics')
  }
  return(out)
}
find.w <- function(ci, s){
  ci[ci<0] <- 0
  s <- trunc(s)
  if (s == 0)
    stop('s cannot be zero!!!!')
  ci.sorted <- sort(ci, index.return=TRUE, decreasing=TRUE)
  x <- rep(0, length(ci))
  x[ci.sorted$ix[1:s]] <- 1
  return(x)
}
