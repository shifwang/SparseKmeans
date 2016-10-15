# Run second experiment parallelly 
# data is generated from mixture Gaussian with correlated features
# mu = 0.6 or 0.7
# p  = 200, 500, 1000
give.data <- function (SamplesPerCluster=20, Nfeatures=50, Nsignals=10,ClusterNumber=3,mu=1,rho=0.1) {
  covr=matrix(NA,ncol=Nfeatures,nrow=Nfeatures)
  for (i in 1:Nfeatures)
    for (j in 1:Nfeatures){
      covr[i,j]=rho^(abs(i-j))
      if (i == j){
      	covr[i, j] = covr[i, j] + runif(1, 0.75, 1.25)
    }
  x <- matrix(NA,ncol=Nfeatures,nrow=SamplesPerCluster*ClusterNumber)
  for (i in 1:ClusterNumber){
    x[(SamplesPerCluster*(i-1)+1):(SamplesPerCluster*i),1:Nfeatures] <- mvrnorm(n=SamplesPerCluster,c(rep(i*mu,Nsignals),rep(0,Nfeatures-Nsignals)),Sigma=covr)
  }
  x <- scale(x, TRUE, TRUE)
  return(x)
}


get.result<- function(alg,x){
  if (alg==1){
    alg.name<-'Standard Kmeans'
    k<-kmeans(x=x,centers=6,nstart=20)
    return(list(name=alg.name,Cs=k[[1]],weights=rep(1,ncol(x))))
  }
  else if (alg==2){
    alg.name<- 'Standard L1'
    source('Algorithms/L1_correction.R', encoding='UTF-8')
    bestl1<-KMeansSparseCluster.permute(x=x,K=6,nvals=100, silent = T)
    l1<-KMeansSparseCluster(x=x,K=6,wbounds=(bestl1$bestw), silent = T)
    #     l1<-KMeansSparseCluster(x=x,K=3,wbounds=10)
    return(list(name=alg.name,Cs=l1[[1]]$Cs,weights=l1[[1]]$ws ))
  }
  else if (alg==3){
    alg.name<- 'original L0'
    source('Algorithms/L0_original.R')
    best0<-select.bound(x=x,K=6,nvals=100)
    out<-give.cluster(x=x,K=6,wbounds=best0$bestw)
    #     out<-give.cluster(x=x,K=3,wbounds=50)
    return(list(name=alg.name,Cs=out[[1]]$Cs,weights=out[[1]]$weight ))
  }
  else if (alg == 4){
    alg.name = 'PCA-Kmeans'
    source('Algorithms/PCA-Kmeans.R')
    out = PCA.kmeans(x = x, K = 6)
    return(list(name=alg.name,Cs=out$Cs,weights=out$weight))
  }
  else if (alg == 5) {
    alg.name <- 'EM'
    source('Algorithms/EM_estimate_mixture_Gaussian.R')
    out <- SelectLambda(x = x, k = 6, nvals = 20, verbose = F)
    out1 <- EstimateMixtureGaussian(data = x, k = 6, lambda = out$best.lambda, verbose = F)
    return(list(name=alg.name, Cs=out1$partition, weights = (1 - out1$noise.feature)))
  }
}

mydivision <- function(testx, x, weights,Cs){
  cluster<- array(NA,nrow(testx))
  for (i in 1:nrow(testx)){
    y<-array(NA,nrow(x))
    for (j in 1:nrow(x)){
      y[j]<- sum((x[j,]*weights-testx[i,]*weights)^2)
    }
    cluster[i] <- Cs[which.min(y)]
  }
  return(cluster)
}
now = proc.time()[[3]]
print <- function(num,str){
  #cat(paste(rep('\b',num),collapse = ''))
  str = paste(str,sprintf(' Passed %1.1f min..',(proc.time()[[3]] - now)/60))
  cat(str,'\n')
  return(nchar(str))
}
criteria.num      = 6    # Number of criteria
iter.num          = 50   # Run how many times to estimate the variation
SamplesPerCluster = 20   # How many samples in one cluster
Nsignals          = 50   # Signal features
ClusterNumber     = 6    # Number of Cluster
rivals.num        = 5    # Number of algorithms
Nfeatures         = c(200, 500, 1000) # Number of features
mus               = c(0.6, 0.7)
verbose           = TRUE
set.seed(11)
total.exp         = c(iter.num, length(Nfeatures), length(mus))
info              = array(list(), total.exp)
true.label        = c()
for (i in 1:ClusterNumber){
  true.label = c(true.label, rep(i,SamplesPerCluster))
}

# Main Part
if (verbose) {
  len = 0
  len = print(len,'Setting working directory...')
}
setwd('./')
if (verbose) {
  len = print(len,' Collecting Data...')
}
library('MASS')
library('foreach')
library('doParallel')
ncores = 12 # Use two cores, use lscpu/nproc to check availabel cpus
registerDoParallel(ncores)

for (mu_ind in 1:length(mus)){
  for (Nf_ind in 1:length(Nfeatures)){
    tmp <- foreach(iter = 1:iter.num) %dopar% {
      x = give.data(SamplesPerCluster,Nfeatures[Nf_ind],Nsignals,ClusterNumber,mus[mu_ind],
                    rho = 0.1)
      return(list(iter, x))
    }
    for (item in tmp) {
      info[item[[1]],Nf_ind,mu_ind][[1]]$data = item[[2]]
      info[item[[1]], Nf_ind, mu_ind][[1]]$true.label = true.label
      info[item[[1]], Nf_ind, mu_ind][[1]]$relevant.feature = c(rep(T, Nsignals),rep(F,Nfeatures[Nf_ind] - Nsignals))
    }
  }
}
if (verbose) 
  len = print(len, ' Saving to disk...')
save(info,file = 'data.RData')
system('mv data.RData $(date "+%F-%H-%M-%S-info-1").RData ')
if (verbose) 
  len = print(len, 'Run Experiments...')

for (mu_ind in 1:length(mus)){
  for (Nf_ind in 1:length(Nfeatures)){
    tmp = foreach(iter=1:iter.num) %dopar%{
      set.seed(100)
      if (verbose) 
        len = print(len, paste('Exp...','mu...',mus[mu_ind],'Nfeature...',Nfeatures[Nf_ind],'iter...',iter))
      x = info[iter, Nf_ind, mu_ind][[1]]$data
      #x = scale(x,TRUE,TRUE) #FIXME: already scaled in give.data
      results = list()
      for (alg in 1:rivals.num){
        results[[alg]] = get.result(alg,x)
      }
      return(list(iter,results))
    }
    for (item in tmp){
      info[item[[1]], Nf_ind, mu_ind][[1]]$results = item[[2]]
    }
  }
}
if (verbose) 
  len = print(len, ' Saving to disk...')
save(info,file = 'data.RData')
system('mv data.RData $(date "+%F-%H-%M-%S-info-2").RData ')
if (verbose) 
  len = print(len, 'Done.\n')

#
#weight_average=array(NA,Nfeatures)
#for (i in 1:Nfeatures){
#  weight_average[i] = mean(weights.info[i,,2])
#}
#cers=array(NA,c(iter.num,3))
#for (i in 1:iter.num){
#  cers[i,1] = CER(partition.info[,i,1],true.label)
#  cers[i,2] = CER(partition.info[,i,2],true.label)
#  cers[i,3] = CER(partition.info[,i,3],true.label)
#  #cers[i,4] = CER(partition.info[,i,4],true.label)
#}
#save(cers,partition.info,weights.info,original.data,true.label,file = 'Day??')
#boxplot(cers,names=list('k-means','sparse k-means'),ylab='CER')
#plot(1:Nfeatures,weight_average,xlab='Features',ylab='weight')

