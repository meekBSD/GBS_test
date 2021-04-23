#!/usr/bin/env Rscript

setwd("D:\\docs_00\\probes\\new_panel\\bedfile_check\\coverage_plots\\")
d <- read.table("Result1.csv", header=F, sep="\t")

colnames(d) <- c("RN", "interval", "P01", "P02", "P03","P04", "P05", "N01", "N02" , "N03")


d1 <- d[order(d$P01, decreasing=F),]

pdf("res_coverageCompare.pdf" ,  width=12.5, height=7)

plot(d1$P01, pch=19, lty=1, type="l", col="grey", xlab="regions", ylab="depth", ylim=c(0,4500))


lines(d1$P02,type="l",col="grey",lwd=2, lty=2)
lines(d1$P03,type="l",col="grey",lwd=2, lty=3)

lines(d1$P04,type="l",col="grey",lwd=2, lty=4) 

lines(d1$P05,type="l",col="grey",lwd=2, lty=5) 
lines(d1$N01,type="l",col="orange",lwd=3, lty=6) 
lines(d1$N02,type="l",col="orange",lwd=3, lty=6) 
lines(d1$N03,type="l",col="orange",lwd=3, lty=6) 

legend("topright",                                    #图例位置为右上角
 legend=c("P01","P02","P03","P04", "P05", "N01", "N02", "N03"),        #图例内容
 col=c("black","black","black","black", "black","orange", "orange", "orange" ),                 #图例颜色
 lty=c(1,2,3,4,5,6,6,6),lwd=2) 

axis(1,label=seq(0, 200, 50), at=seq(0, 200, 50))

dev.off()
