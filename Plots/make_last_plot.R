# make plots for the last experiment
CER<- function (partition1=NULL,partition2=NULL){
  if (length(partition1)!= length(partition2))
    stop("partition1 and 2 don't have the same observations")
  mycer<-0;
  for (i in 1:(length(partition1)-1)){
    for (j in (i+1):length(partition1)){
      mycer<- mycer + abs((partition1[i]==partition1[j])-(partition2[i]==partition2[j])    )  
    }
  }
  mycer<- mycer/(length(partition1)*(length(partition1)-1)/2)
  return(mycer)
}
CER.fast<- function (partition1=NULL,partition2=NULL){
  if (length(partition1)!= length(partition2))
    stop("partition1 and 2 don't have the same observations")
  mycer <- 0;
  tmp <- table(partition1, partition2)
  tmp.x <- dim(tmp)[1]
  tmp.y <- dim(tmp)[2]
  n <- length(partition1)
  wrong.sum <- matrix(0, tmp.x, tmp.y)
  for (row in 1:tmp.x) {
    for (col in 1:tmp.y)
    wrong.sum[-row, col] <- wrong.sum[-row, col] + tmp[row, col]
  }
  for (row in 1:tmp.x){
    for (col in 1:tmp.y) {
      wrong.sum[row, -col] <- wrong.sum[row, -col] + tmp[row, col]
    }  
  }
  mycer <- sum(tmp * wrong.sum)/(n*(n-1)) 
  return(mycer)
}

load('~/Documents/kmeans/SparseKmeans/bio_pca_2016-05-20-09-36-08.RData') #  load the data

for (item in info) {
  results <- item[[2]]
  true.label <- item[[1]]
  print(CER.fast(results$Cs$cluster, true.label))
}
