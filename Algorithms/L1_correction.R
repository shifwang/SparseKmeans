#2013年12月16日修改，用于gapstatistics说明
#稍微修改了原文聚类的方法
#修改了每次不同参数对应的Cs的重置
KMeansSparseCluster <- function (x, K = NULL, wbounds = NULL, nstart = 20, silent = FALSE, 
                                 maxiter = 6, centers = NULL, seed = 123) 
{
  set.seed(seed)
  if (is.null(K) && is.null(centers)) 
    stop("Must provide either K or centers.")
  if (!is.null(K) && !is.null(centers)) {
    if (nrow(centers) != K) 
      stop("If K and centers both are provided, then nrow(centers) must equal K!!!")
    if (nrow(centers) == K) 
      K <- NULL
  }
  if (!is.null(centers) && ncol(centers) != ncol(x)) 
    stop("If centers is provided, then ncol(centers) must equal ncol(x).")
  if (is.null(wbounds)) 
    wbounds <- seq(1.1, sqrt(ncol(x)), len = 20)
  if (min(wbounds) <= 1) 
    stop("wbounds should be greater than 1")
  wbounds <- c(wbounds)
  out <- list()
  for (i in 1:length(wbounds)) {
    if (!is.null(K)) 
      Cs <- kmeans(x, centers = K, nstart = nstart)$cluster
    if (is.null(K)) 
      Cs <- kmeans(x, centers = centers)$cluster
    if (length(wbounds) > 1 && !silent) 
      cat(i, fill = FALSE)
    ws <- rep(1/sqrt(ncol(x)), ncol(x))
    ws.old <- rnorm(ncol(x))
    store.bcss.ws <- NULL
    niter <- 0
    while ((sum(abs(ws - ws.old))/sum(abs(ws.old))) > 1e-04 && 
             niter < maxiter) {
      if (!silent) 
        cat(niter, fill = FALSE)
      niter <- niter + 1
      ws.old <- ws
      if (!is.null(K)) {
        if (niter > 1) 
          Cs <- UpdateCs(x, K, ws, Cs)
      }
      else {
        if (niter > 1) 
          Cs <- UpdateCs(x, nrow(centers), ws, Cs)
      }
      ws <- UpdateWs(x, Cs, wbounds[i])
      store.bcss.ws <- c(store.bcss.ws, sum(GetWCSS(x, 
                                                    Cs)$bcss.perfeature * ws))
    }
    out[[i]] <- list(ws = ws, Cs = Cs, wcss = GetWCSS(x, 
                                                      Cs, ws), crit = store.bcss.ws, wbound = wbounds[i])
  }
  if (!silent) 
    cat(fill = TRUE)
  class(out) <- "KMeansSparseCluster"
  return(out)
}

UpdateWs <- function (x, Cs, l1bound) 
{
  wcss.perfeature <- GetWCSS(x, Cs)$wcss.perfeature
  tss.perfeature <- GetWCSS(x, rep(1, nrow(x)))$wcss.perfeature
  lam <- BinarySearch(-wcss.perfeature + tss.perfeature, l1bound)
  ws.unscaled <- soft(-wcss.perfeature + tss.perfeature, lam)
  return(ws.unscaled/l2n(ws.unscaled))
}


UpdateCs <- function (x, K, ws, Cs) 
{
  if (sum(ws != 0) == 1) {
    only.one.feature <- T
  } else {
    only.one.feature <- F
  }
  if (only.one.feature) {
    x <- matrix(x[, ws != 0], ncol = 1)
  } else {
    x <- x[, ws != 0]
  }
  
  z <- sweep(x, 2, sqrt(ws[ws != 0]), "*")
  nrowz <- nrow(z)
  mus <- NULL
  if (!is.null(Cs)) {
    for (k in unique(Cs)) {
      if (sum(Cs == k) > 1) {
        if (only.one.feature) {
          mus <- rbind(mus, apply(matrix(z[Cs == k, ], ncol = 1), 2, mean))
        } else {
          mus <- rbind(mus, apply(z[Cs == k, ], 2, mean))
        }
      }
        
      if (sum(Cs == k) == 1) 
        mus <- rbind(mus, z[Cs == k, ])
    }
  }
  if (is.null(mus)) {
    km <- kmeans(z, centers = K, nstart = 10)
  }
  else {
    distmat <- as.matrix(dist(rbind(z, mus)))[1:nrowz, (nrowz + 
                                                          1):(nrowz + K)]
    nearest <- apply(distmat, 1, which.min)
    if (length(unique(nearest)) == K) {
      km <- kmeans(z, centers = mus)
    }
    else {
      km <- kmeans(z, centers = K, nstart = 10)
    }
  }
  return(km$cluster)
}

l2n <- function (vec) 
{
  return(sqrt(sum(vec^2)))
}

BinarySearch <- function (argu, sumabs) 
{
  if (l2n(argu) == 0 || sum(abs(argu/l2n(argu))) <= sumabs) 
    return(0)
  lam1 <- 0
  lam2 <- max(abs(argu)) - 1e-05
  iter <- 1
  while (iter <= 15 && (lam2 - lam1) > (1e-04)) {
    su <- soft(argu, (lam1 + lam2)/2)
    if (sum(abs(su/l2n(su))) < sumabs) {
      lam2 <- (lam1 + lam2)/2
    }
    else {
      lam1 <- (lam1 + lam2)/2
    }
    iter <- iter + 1
  }
  return((lam1 + lam2)/2)
}

soft <- function (x, d) 
{
  return(sign(x) * pmax(0, abs(x) - d))
}
GetWCSS <- function (x, Cs, ws = NULL) 
{
  wcss.perfeature <- numeric(ncol(x))
  for (k in unique(Cs)) {
    whichers <- (Cs == k)
    if (sum(whichers) > 1) 
      wcss.perfeature <- wcss.perfeature + apply(scale(x[whichers, 
                                                         ], center = TRUE, scale = FALSE)^2, 2, sum)
  }
  bcss.perfeature <- apply(scale(x, center = TRUE, scale = FALSE)^2, 
                           2, sum) - wcss.perfeature
  if (!is.null(ws)) 
    return(list(wcss.perfeature = wcss.perfeature, wcss = sum(wcss.perfeature), 
                wcss.ws = sum(wcss.perfeature * ws), bcss.perfeature = bcss.perfeature))
  if (is.null(ws)) 
    return(list(wcss.perfeature = wcss.perfeature, wcss = sum(wcss.perfeature), 
                bcss.perfeature = bcss.perfeature))
}

KMeansSparseCluster.permute<-function (x, K = NULL, nperms = 1, wbounds = NULL, silent = FALSE, 
                                       nvals = 10, centers = NULL, seed = 101) 
{
  set.seed(seed)
  if (is.null(wbounds)) 
    wbounds <- exp(seq(log(1.2), log(sqrt(ncol(x)) * 0.9), 
                       len = nvals))
  if (min(wbounds) <= 1) 
    stop("Wbounds should be greater than 1, since otherwise only one weight will be nonzero.")
  if (length(wbounds) < 2) 
    stop("Wbounds should be a vector of at least two elements.")
  if (is.null(K) && is.null(centers)) 
    stop("Must provide either K or centers.")
  if (!is.null(K) && !is.null(centers)) {
    if (nrow(centers) != K) 
      stop("If K and centers both are provided, then nrow(centers) must equal K!!!")
    if (nrow(centers) == K) 
      K <- NULL
  }
  if (!is.null(centers) && ncol(centers) != ncol(x)) 
    stop("If centers is provided, then ncol(centers) must equal ncol(x).")
  permx <- list()
  nnonzerows <- NULL
  for (i in 1:nperms) {
    permx[[i]] <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    for (j in 1:ncol(x)) permx[[i]][, j] <- sample(x[, j])
  }
  tots <- NULL
  out <- KMeansSparseCluster(x, K, wbounds = wbounds, silent = silent, 
                             centers = centers)
  for (i in 1:length(out)) {
    nnonzerows <- c(nnonzerows, sum(out[[i]]$ws != 0))
    bcss <- GetWCSS(x, out[[i]]$Cs)$bcss.perfeature
    tots <- c(tots, sum(out[[i]]$ws * bcss))
  }
  permtots <- matrix(NA, nrow = length(wbounds), ncol = nperms)
  for (k in 1:nperms) {
    if (!silent) 
      cat("Permutation ", k, "of ", nperms, fill = TRUE)
    perm.out <- KMeansSparseCluster(permx[[k]], K, wbounds = wbounds, 
                                    silent = silent, centers = centers)
    for (i in 1:length(perm.out)) {
      perm.bcss <- GetWCSS(permx[[k]], perm.out[[i]]$Cs)$bcss.perfeature
      permtots[i, k] <- sum(perm.out[[i]]$ws * perm.bcss)
    }
  }
  gaps <- (log(tots) - apply(log(permtots), 1, mean))
  plot(nnonzerows,gaps,xlab='Number of Non-zero Weights-L1',ylab='Gap Statistics')
  out <- list(tots = tots, permtots = permtots, nnonzerows = nnonzerows, 
              gaps = gaps, sdgaps = apply(log(permtots), 1, sd), wbounds = wbounds, 
              bestw = wbounds[which.max(gaps)])
  if (!silent) 
    cat(fill = TRUE)
  class(out) <- "KMeansSparseCluster.permute"
  return(out)
}