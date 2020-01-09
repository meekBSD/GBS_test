#!/usr/bin/env Rscript


dat <- read.table('result_dep.txt', header=F, sep='\t')

colnames(dat) <- c('chr', 'x', 'y', 'region')

attach(dat)

#pdf("res_1.pdf", width = 25, height=6)
png("res_1.png", width = 2500, height=600)
plot(x, y, col=c("red","blue")[region])
dev.off()
detach(dat)


