## Initial analysis

library(ggplot2)
library(lme4)
library(tidyr)

rm(list=ls())

## Read in data
pestdat <- read.csv("../data/Pestwoodyseedlings_20170116.csv")


pestsurv <- aggregate(cbind(height_15f, height_16s) ~ 
                        site + code + quadrat + species, data=pestdat,
                      function(x) sum(!is.na(x) & x> 0))
pestdeaths <- aggregate(cbind(height_16s, height_16f) ~ 
                        site + code + quadrat + species, data=pestdat,
                      function(x) sum(!is.na(x) & x==0))
names(pestsurv)[5:6] <- c('total_16s', 'total_16f') 
names(pestdeaths)[5:6] <- c('deaths_16s', 'deaths_16f') 
pestdat.ag <- merge(pestsurv, pestdeaths, by=c('site', 'code', 'quadrat', 'species'))

pestdat.ag <-  gather(pestdat.ag, "meas", "N", total_16s:deaths_16f) %>% 
  separate(meas, into=c("meas", "census"), sep="_") %>%
  spread(meas, N)
head(pestdat.ag)

pestdat.ag$survs <- pestdat.ag$total - pestdat.ag$deaths


pestdat.ag2 <- separate(pesdat.ag2, meas
head(pestdat.ag2)



dat$surv1 <-as.numeric( dat$height16a !=0 & !is.na(dat$height16a))
dat$surv2 <- as.numeric(dat$height16b !=0 & !is.na(dat$height16b))
quad <-  as.factor(dat$code):dat$sub.code

dat.ag <-  aggregate(cbind(surv1, surv2) ~ treatment + exclosure + species + quad, data=dat, FUN=mean)
summary(dat.ag)
ggplot(data=dat.ag, aes(x=treatment, y=surv1, colour=as.factor(exclosure))) + geom_point(position='jitter')
ggplot(data=dat.ag, aes(x=treatment, y=surv1, colour=as.factor(exclosure))) + geom_violin()

## NEW
