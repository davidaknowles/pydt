system("./ts 1")
system("./ts 0")

a=read.table("slice_results.txt")
g1=a$V1
plot(a$V1,a$V2,col="red")
b=read.table("mh_results.txt")
lines(b$V1,b$V2)