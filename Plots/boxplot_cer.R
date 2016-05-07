# Generate plots using output file from experiments

result.cer=array(NA,dim=c(5,4,20,4))
for (i1 in 1:5)
  for (i2 in 1:4)
    for (i in 1:20)
      for (j in 1:4)
      {result.cer[i1,i2,i,j]=CER(true.label,partition.info[,i1,i2,i,j])}
#boxplot(wangyu)
#boxplot for simulaiotn 1 in sparse clustering paper
par(mfrow=c(1,3),mai=c(.8,.8,.1,.1))
mat<-cbind(result.cer[1,2,,1:2],result.cer[1,2,,4])
boxplot(mat,col=c(0,0,2),boxwex=0.4,ylab=expression(mu~'='~0.6),names=c("k-means",expression(italic(l)["1"]),expression(italic(l)["0"]), "PCA k-means", "P-likelihood"))

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
