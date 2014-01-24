system("./gir 1 1000 6 > postsamples.txt")
system("./gir 0 1000 6 > priorsamples.txt")

a=read.table("postsamples.txt")
#a=read.table("ps_results.txt")
g1=a$V1
plot(ecdf(g1),col="red")
b=read.table("priorsamples.txt")
g2=b$V1
lines(ecdf(g2))
g_bar1=mean(g1)
g_bar2=mean(g2)
sigma_g1=sqrt(var(g1))
sigma_g2=sqrt(var(g2))
M1=length(g1)
M2=length(g2)
test_stat=(g_bar1-g_bar2) /    (sigma_g1^2/M1+sigma_g2^2/M2)^.5
print(pnorm(-abs(test_stat)))