rm(list = ls(all = TRUE))
fn="depth_vs_n.txt"
system(paste("rm",fn))
# theta alpha c
theta=0.0
alpha=0.5
c=1.0
system(paste("./func_of_n.bin",theta,alpha,c,">>",fn))
a=read.table(fn)
#a=a[100:nrow(a),]
attach(a)
plot(V1,V2,"b",xlab="n",ylab="mean depth",main=paste("100 repeats, theta=",theta,", alpha=",alpha,", c=",c,sep=""))
dev.print(pdf,paste("depth_vs_n_plots/rep100_theta=",theta,"_alpha=",alpha,"_c=",c,".pdf",sep=""))
a=a[100:nrow(a),]
colnames(a)=c("n","d")
a$logn=log(a$n)
a$logd=log(a$d)
l=lm(logd~logn,data=a)
s=summary(l)
print(s)
attach(a)
plot(logn,logd,"b",xlab="log n",ylab="log mean depth",main=paste("100 repeats, theta=",theta,", alpha=",alpha,", c=",c,", slope=",format(s[[4]][2,1],digits=2),sep=""))
dev.print(pdf,paste("depth_vs_n_plots/loglog_rep100_theta=",theta,"_alpha=",alpha,"_c=",c,".pdf",sep=""))