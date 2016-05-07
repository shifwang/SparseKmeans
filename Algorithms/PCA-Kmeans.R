PCA.kmeans <- function (x=NULL,K=6){
  a <- prcomp(t(x))
  a<- a[[5]][,1]
  xx<-x%*%as.matrix(a)
  y<- kmeans(xx,centers=K,nstart=20)
  return(list(Cs=y$cluster,weight=a))
}
  