a=read.table("toplot.txt")
upper=50
a=a[1:upper,]
a$V2[10]=1000
attach(a)
V2=V2/sum(V2)
s=sqrt(sum((V1^2)*V2)/2)
plot(V1,V2,"l",xlab="depth",ylab="freq")
k=read.table("kingman1000x1000.txt")
k=k[1:upper,]
attach(k)
V2=V2/sum(V2)
s=sqrt(sum((V1^2)*V2)/2)
lines(V2~V1,lty=2,col="green")
#require(VGAM)
#points(function(x) drayleigh(x,s), 0, upper,add=T,col="red")
dev.print(pdf,"kingman.pdf")