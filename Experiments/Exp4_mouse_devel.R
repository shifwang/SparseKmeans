# Run experiments on mouse brain atlas data set
get.result<- function(alg, x, cluster.number = 6){
  if (alg==1){
    alg.name<-'Standard Kmeans'
    k<-kmeans(x = x, centers = cluster.number, nstart = 1)
    return(list(name=alg.name,Cs=k,weights=rep(1,ncol(x))))
  }
  else if (alg==2){
    alg.name<- 'Standard L1'
    source('Algorithms/L1_correction.R', encoding='UTF-8')
    bestl1<-KMeansSparseCluster.permute(x=x,K=cluster.number,nvals=100, silent = T)
    l1<-KMeansSparseCluster(x=x,K=cluster.number,wbounds=(bestl1$bestw), silent = T)
    #     l1<-KMeansSparseCluster(x=x,K=3,wbounds=10)
    return(list(name=alg.name,Cs=l1[[1]]$Cs,weights=l1[[1]]$ws ))
  }
  else if (alg==3){
    alg.name<- 'original L0'
    source('Algorithms/L0_original.R')
    best0 <- select.bound(x = x, K = cluster.number, nvals = 4, nperms = 2, verbose = T)
    Cs <- kmeans(x = x[, best0$signal.feature], centers = cluster.number, nstart = 2)
    #     out<-give.cluster(x=x,K=3,wbounds=50)
    return(list(name = alg.name, Cs = Cs, weights = best0$signal.feature ))
  }
  else if (alg == 4){
    alg.name = 'PCA-Kmeans'
    source('Algorithms/PCA-Kmeans.R')
    out = PCA.kmeans(x = x, K = cluster.number)
    return(list(name=alg.name,Cs=out$Cs,weights=out$weight))
  }
  else if (alg == 5) {
    alg.name <- 'EM'
    source('Algorithms/EM_estimate_mixture_Gaussian.R')
    #  out <- SelectLambda(x = x, k = cluster.number, nvals = 2, verbose = F)
    out1 <- EstimateMixtureGaussian(data = x, k = cluster.number, lambda = 1, verbose = F)
    return(list(name=alg.name, Cs=out1$partition, weights = (1 - out1$noise.feature)))
  }
}
now = proc.time()[[3]]
print <- function(num,str){
  #cat(paste(rep('\b',num),collapse = ''))
  str = paste(str,sprintf(' Passed %1.1f min..',(proc.time()[[3]] - now)/60))
  cat(str,'\n')
  return(nchar(str))
}

iter.num          = 2   # Run how many times to estimate the variation
cluster.number    = 19    # Number of Cluster
rivals.num        = 5    # Number of algorithms
verbose           = TRUE
total.exp         = c(1, iter.num)
info              = array(list(), total.exp)
true.label        = c()
library('R.matlab')
library('foreach')
library('doParallel')
ncores = 5 #  use lscpu/nproc to check availabel cpus
registerDoParallel(ncores)
tmp <- foreach(ind = 1:1) %dopar% {
  raw.data <- readMat(paste(c('../../SparseKmeans/data/dmatrix_', ind, '.mat'),collapse = ''))
  true.label <- raw.data$up.annot
  
  # Main Part
  if (verbose) 
    len <- print(len, paste(c('Run Bio Experiment ', ind, '...'), collapse = ''))
  set.seed(11)
  x <- t(raw.data$dmatrix)
  # x <- scale(x, TRUE, TRUE) #  Maybe not?
  # delete NAs, because there are features that is always the same
  x <- x[,apply(is.na(x),2,sum) == 0]
  info <- list()
  info$true.label <- raw.data$up.annot
  info$results <-  get.result(5, x, cluster.number = cluster.number)
  return(list(ind, info))
}
info = list()
for (item in tmp) {
  info[[item[[1]]]] <- item[[2]]
}
if (verbose) 
len = print(len, ' Saving to disk...')
save(info,file = 'data.RData')
system(paste(c('mv data.RData ', 'bio_EM_','$(date "+%F-%H-%M-%S").RData '),collapse = ''))
if (verbose) 
  len = print(len, 'Done.\n')
