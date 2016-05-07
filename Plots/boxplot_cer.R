# Generate plots using output file from experiments

load('2016-05-07-16-45-09-info-2.RData')
num.iter <- dim(info)[1]
num.feature <- dim(info)[2]
num.mus <- dim(info)[3]
num.alg <- length(info[1,1,1][[1]]$results)
result.cer=array(NA, dim = c(num.iter, num.feature, num.mus, num.alg))
for (i1 in 1:num.iter)
  for (i2 in 1:num.feature)
    for (i3 in 1:num.mus)
      for (i4 in 1:num.alg) {
        tmp <- info[i1, i2, i3][[1]]
        result.cer[i1, i2, i3, i4] <- CER(tmp$true.label,tmp$results[[i4]]$Cs) 
      }
alg.names <- c("k-means",expression(italic(l)["1"]),expression(italic(l)["0"]), "PCA k-means", "P-likelihood")
#boxplot(wangyu)
#boxplot for simulaiotn 1 in sparse clustering paper
par(mfrow=c(1,3),mai=c(.8,.8,.1,.1))
mat<-cbind(result.cer[,, 1, 1])
boxplot(mat, col = c(0, 0, 2, 0, 0), 
        boxwex = 0.4, ylab = expression(mu~'='~0.6), 
        names = alg.names)

mat<-cbind(result.cer[1,3,,1:2],result.cer[1,3,,4])
boxplot(mat,col=c(0,0,2),boxwex=0.4,names=c("k-means",expression(italic(l)["1"]),expression(italic(l)["0"]), "PCA k-means", "P-likelihood"))

mat<-cbind(result.cer[1,4,,1:2],result.cer[1,4,,4])
boxplot(mat,col=c(0,0,2),boxwex=0.4,names=c("k-means",expression(italic(l)["1"]),expression(italic(l)["0"]), "PCA k-means", "P-likelihood"))

par(mfrow=c(1,3),mai=c(.8,.8,.1,.1))
mat<-cbind(result.cer[2,2,,1:2],result.cer[2,2,,4])
boxplot(mat,col=c(0,0,2),boxwex=0.4,xlab="p=200",ylab=expression(mu~'='~0.7),names=c("k-means",expression(italic(l)["1"]),expression(italic(l)["0"]), "PCA k-means", "P-likelihood"))


mat<-cbind(result.cer[2,3,,1:2],result.cer[2,3,,4])
boxplot(mat,col=c(0,0,2),boxwex=0.4,xlab="p=500",xlab="p=50",names=c("k-means",expression(italic(l)["1"]),expression(italic(l)["0"]), "PCA k-means", "P-likelihood"))

mat<-cbind(result.cer[2,4,,1:2],result.cer[2,4,,4])
boxplot(mat,col=c(0,0,2),boxwex=0.4,xlab="p=1000",names=c("k-means",expression(italic(l)["1"]),expression(italic(l)["0"]), "PCA k-means", "P-likelihood"))
