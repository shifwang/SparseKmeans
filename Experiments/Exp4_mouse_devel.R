# Run experiments on mouse brain atlas data set
get.result<- function(alg, x, cluster.number = 6){
  if (alg==1){
    alg.name<-'Standard Kmeans'
    k<-kmeans(x = x, centers = cluster.number, nstart = 20)
    return(list(name=alg.name,Cs=k[[1]],weights=rep(1,ncol(x))))
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
    best0 <- select.bound(x = x, K = cluster.number, nvals = 10, verbose = T )
    Cs <- kmeans(x = x[, best0$signal.feature], K = cluster.number, nstart = 2 )
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
    out <- SelectLambda(x = x, k = cluster.number, nvals = 20, verbose = F)
    out1 <- EstimateMixtureGaussian(data = x, k = cluster.number, lambda = out$best.lambda, verbose = F)
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
raw.data <- readMat(paste(c('../SparseKmeans/data/dmatrix_', 1, '.mat'),collapse = ''))
true.label <- raw.data$up.annot

# Main Part
if (verbose) {
  len = 0
  len = print(len,'Setting working directory...')
}
setwd('./')
if (verbose) {
  len = print(len,' Collecting Data...')
}

library('foreach')
library('doParallel')
# ncores = 2 # Use two cores, use lscpu/nproc to check availabel cpus
# registerDoParallel(ncores)


if (verbose) 
  len <- print(len, 'Run Experiments...')
set.seed(11)
x <- raw.data$dmatrix
x <- scale(x, TRUE, TRUE) #  Maybe not?
results <-  get.result(3, x, cluster.number = cluster.number)
if (verbose) 
  len = print(len, ' Saving to disk...')
save(results,file = 'data.RData')
system(paste(c('mv data.RData ', 'bio_devel_','$(date "+%F-%H-%M-%S").RData '),collapse = ''))
if (verbose) 
  len = print(len, 'Done.\n')
