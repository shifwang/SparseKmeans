# Generate plots using output file from experiments
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

#load('2016-05-08-06-33-48-info-2.RData') #  load Exp2 data
#load('2016-05-08-14-17-15-info-2.RData')  #  load Exp3 data
#load('2016-05-09-00-19-31-info-2.RData')  #  load Exp4 data log normal
load('2016-05-09-12-58-07-info-2.RData') #  load Exp4 data pois

num.iter <- dim(info)[1]
num.feature <- dim(info)[2]
num.mus <- dim(info)[3]
num.alg <- length(info[1,1,1][[1]]$results)
result.cer <- array(NA, dim = c(num.iter, num.feature, num.mus, num.alg))
result.recall <- array(NA, dim = c(num.iter, num.feature, num.mus, num.alg))
result.precision <- array(NA, dim = c(num.iter, num.feature, num.mus, num.alg))
result.Fscore <- array(NA, dim = c(num.iter, num.feature, num.mus, num.alg))
true.label <- rep(1:6,rep(20,6))
for (i1 in 1:num.iter)
  for (i2 in 1:num.feature)
    for (i3 in 1:num.mus)
      for (i4 in 1:num.alg) {
        tmp <- info[i1, i2, i3][[1]]
        result.cer[i1, i2, i3, i4] <- CER(true.label,tmp$results[[i4]]$Cs) 
        result.recall[i1, i2, i3, i4] <- sum(tmp$results[[i4]]$weights[1:50] != 0)/50
        result.precision[i1, i2, i3, i4] <- sum(tmp$results[[i4]]$weights[1:50] != 0)/sum(tmp$results[[i4]]$weights != 0)
        result.Fscore[i1, i2, i3, i4] <- 2/(result.precision[i1, i2, i3, i4]^-1 + result.recall[i1, i2, i3, i4]^-1)
      }
alg.names <- c("k-means",expression(italic(l)["1"]),expression(italic(l)["0"]), "PCA", "EM")

MakeBoxplot <- function(result, 
                        ylim = c(0,1), 
                        mus  = c(0.6, 0.7),
                        alg.names){
  # boxplot for simulaiotn 1 in sparse clustering paper
  par(mfrow = c(1,3), mai = c(.8, .8, .1, .1))
  mat <- result[, 1, 1, ]
  ylab <- as.formula(paste('mu ~ "=" ~ ',mus[1]))
  boxplot(mat, col = c(0, 0, 2, 0, 0), 
          boxwex = 0.4, ylab = as.expression(ylab), 
          ylim = ylim,
          names = alg.names)
  
  mat <- result[, 2, 1, ]
  boxplot(mat, col = c(0, 0, 2, 0, 0), 
          ylim = ylim,
          boxwex = 0.4, names = alg.names)
  
  mat <- result[, 3, 1, ]
  boxplot(mat, col = c(0, 0, 2, 0, 0),
          ylim = ylim,
          boxwex = 0.4, names = alg.names)
  
  
  par(mfrow = c(1, 3), mai = c(.8, .8, .1, .1))
  mat <- result[, 1, 2, ]
  ylab <- as.formula(paste('mu ~ "=" ~ ',mus[2]))
  boxplot(mat, col = c(0, 0, 2, 0, 0),
          boxwex = 0.4, xlab = "p = 200",
          ylab = as.expression(ylab),
          ylim = ylim,
          names = alg.names)
  
  mat <- result[, 2, 2, ]
  boxplot(mat,col=c(0, 0, 2, 0, 0),
          boxwex = 0.4, xlab = "p = 500",
          ylim = ylim,
          names=alg.names)
  
  mat <- result[, 3, 2, ]
  boxplot(mat, col = c(0, 0, 2, 0, 0),
          boxwex = 0.4, xlab = "p = 1000",
          ylim = ylim,
          names = alg.names)
}
MakeBoxplot(result.cer, 
            ylim = c(0, 0.2),
            mus = c(2, 3),
            alg.names)
MakeBoxplot(result.Fscore, 
            ylim = c(0, 1),
            mus = c(2, 3),
            alg.names)
