## Initial analysis

library(ggplot2)

dat <- read.csv("../data/Woodyseedlings 20161129.csv")
summary(dat)

dat$surv1 <-as.numeric( dat$height16a !=0 & !is.na(dat$height16a))
dat$surv2 <- as.numeric(dat$height16b !=0 & !is.na(dat$height16b))
quad <-  as.factor(dat$code):dat$sub.code

dat.ag <-  aggregate(cbind(surv1, surv2) ~ treatment + exclosure + species + quad, data=dat, FUN=mean)
summary(dat.ag)
ggplot(data=dat.ag, aes(x=treatment, y=surv1, colour=as.factor(exclosure))) + geom_point(position='jitter')
ggplot(data=dat.ag, aes(x=treatment, y=surv1, colour=as.factor(exclosure))) + geom_violin()

## NEW
