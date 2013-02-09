library(KernSmooth)
library(mosaic)
library(nlme)
library(ggplot2)
library(lattice)
library(reshape2)

data(Soybean)
soy = as.data.frame(Soybean)

qplot(Time, weight, facets=.~Plot, data=soy)


xyplot(weight~Time | Plot, data=soy)


data(Spruce)
plot(Spruce)

xyplot(logSize~days | Tree, data=Spruce)

library(lme4)
data(sleepstudy)

xyplot(Reaction~Days | Subject, data=sleepstudy)


library(timecourse)
data(fruitfly)
m = rowMeans(fruitfly)
sf = apply(fruitfly,1,sd)
oi = order(m, decreasing=TRUE)

dros = fruitfly[oi[1:5],]

plot(dros[1,], ylim=c(8,16), type='l')
lines(dros[2,])
lines(dros[3,])
lines(dros[4,])
lines(dros[5,])

oi2 = order(sf, decreasing=TRUE)
dros = rbind(dros, fruitfly[oi2[c(1:3,6,9)],])

plot(dros[6,], ylim=c(5,16))
lines(dros[7,])
lines(dros[8,])
lines(dros[9,])
lines(dros[10,])


dros = rbind(dros, fruitfly[oi2[c(50:53)],])
plot(dros[11,], ylim=c(8,16), type='l')
lines(dros[12,])
lines(dros[13,])
lines(dros[14,])


dros = cbind(c(rep("group1",5), rep("group2", 5), rep("group3", 4)), dros)
write.csv(dros, "dros.csv", quote=FALSE, row.names=TRUE)

dros=read.csv('dros.csv', header=TRUE)
names(dros)[1:2] = c("gene", "group")
drosmelt = melt(dros)
drosmelt$time = sort(rep(1:12,14))
drosmelt$replicate = c(rep("A", 168), rep("B", 168), rep("C", 168))

xyplot(log2exp~time | gene, data=drosmelt)


write.csv(drosmelt, "droslong.csv", quote=FALSE, row.names=FALSE)
