rm(list = ls(all = TRUE))
system("rm toplot.txt")
# theta alpha c
theta=0
alpha=0
c=1.0
system(paste("./tree_hist.bin",theta,alpha,c,">> toplot.txt"))
a=read.table("toplot.txt")
#a=a[1:40,]
upper=50
a=a[1:upper,]
attach(a)
V2=V2/sum(V2)
s=sqrt(sum((V1^2)*V2)/2)
plot(V1,V2,"b",xlab="depth",ylab="freq",main=paste("N=1e5, theta=",theta,", alpha=",alpha,", c=",c,sep=""))
require(VGAM)
plot(function(x) drayleigh(x,s), 0, upper,add=T,col="red")
dev.print(pdf,paste("rayleigh_plots/","N=1e5_theta=",theta,"_alpha=",alpha,"_c=",c,".pdf",sep=""))