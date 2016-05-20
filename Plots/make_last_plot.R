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
load('~/Documents/kmeans/SparseKmeans/bio_devel_2016-05-20-07-27-59.RData') #  load the data
results <- info[[2]]
true.label <- info[[1]]
sum(results$weights!=0)
CER(results$Cs, true.label)
