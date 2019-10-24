#!usr/bin/env R

dat<-read.delim('all_pep.fsa.xyz.mcl.flt_xy.txt', sep='\t', header=F)

y=dat$V2
x=dat$V1

# core genes
a = 3388.226344382987
b = 0.9869163849684618
c = 84902.67215485258
fc <- function(n)(a * exp(-n/b) + c)

pdf("pan_curve.pdf")
par(mfrow=c(2,2))

#png(filename="curve_core_gene.png")
plot(x,y, xlab='# of genomes', ylab='# of core genes', pch=19)
lines(fc(1:193), col='red', lwd=3)
#dev.off()

# new genes per sequenced
K_s = 1370025070.0290346
Tau_s = 0.20672210045683362
TgTheta = 974.1165151383516
fs <- function(n){K_s * exp(-n / Tau_s) + TgTheta}

y=dat$V3
x=dat$V1

#png(filename="curve_new_gene.png")
#plot(x,y)
plot(x,y, xlab='# of genomes', ylab='# of new genes', pch=19)
lines(fs(1:193), col='blue', lwd=3)
#dev.off()



# pangenome size
K=35352.63563929516
r=0.361502042629523

fp <- function(n){K * n ** r}

y=dat$V4
x=dat$V1

#png(filename="curve_pan_size.png")
#plot(x,y)
plot(x,y, xlab='# of genomes', ylab='size of pan-genome', pch=19)
lines(fp(1:193), col='green', lwd=3)

dev.off()
#dev.off('pan_curve')




#library(ggplot2)
#qplot(x,y, geom='smooth', span =0.5)

