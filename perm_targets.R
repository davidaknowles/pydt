a=read.csv("target_drug_mapping.txt")
b=read.table("drug_names_perm.txt")
n=as.character(b$V1)
p=a[,n]
xLabels=as.character(a$drugnm)

pdf("targets.pdf",width=4,height=8)
image(1:length(xLabels),1:ncol(p),1-as.matrix(p),col=grey.colors(100),axes=F,xlab="",ylab="")
axis(1, at=1:length(xLabels), labels=xLabels,  las=2)
dev.off()