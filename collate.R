prefix=c("slice_0","slice_-1","mh_-1")
labels=c("slice_ddt","slice pydt","mh pydt")
res=matrix(ncol=3,nrow=10)
times=matrix(ncol=3,nrow=10)
for (j in 1:3)
{
    for (i in 1:10)
    {
        a=read.table(paste("cluster_results/",prefix[j],"_rand",i,"_ml.txt",sep=""))
        res[i,j]=a$V4[nrow(a)]
        times[i,j]=a$V2[nrow(a)]
    }
}

