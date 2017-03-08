## Initial analysis

library(ggplot2)
library(Matrix)
library(lme4)
library(permute)
library(lattice)
library(vegan)
library(tidyr)
library(dplyr)
library(pbkrtest)

## dat
## exclude shrubs which height shorter than 30 cm 
dat <- read.csv("Changbaishan/data/Woodyseedlings 20161215.csv")
## dat.all
## the original data
dat.all <- read.csv("Changbaishan/data/Woodyseedlings 20161214.csv")
####
## dat.gathered
## put all the height of each census into one column
## use tidyr to get longer data with all height in one column
gathered.dat.all <- gather(dat.all,heighcensus,height,height15b:height16b)

####
## analyze trees and shrubs seperately
dat.all.tre <- dat.all[dat.all$growth.form == "tree",]
dat.all.shr <- dat.all[dat.all$growth.form == "shrub",]
dat.shr <- dat[dat$growth.form == "shrub",]
## maybe dat.all is better for more shrub sample size included

summary(dat)
summary(dat.all)
summary(gathered.dat.all)
summary(dat.all.tre)
summary(dat.all.shr)

dat$site <- factor(dat$site)
dat$quadrat <- factor(dat$quadrat)
dat$tag <- factor(dat$tag)
dat$exclosure <- factor(dat$exclosure)
#quad <-  as.factor(dat$code):dat$sub.code
#dat <- cbind (dat,quad)
dat.all$site <- factor(dat.all$site)
dat.all$quadrat <- factor(dat.all$quadrat)
dat.all$tag <- factor(dat.all$tag)
dat.all$exclosure <- factor(dat.all$exclosure)
## dat.all.tre
dat.all.tre$site <- factor(dat.all.tre$site)
dat.all.tre$quadrat <- factor(dat.all.tre$quadrat)
dat.all.tre$tag <- factor(dat.all.tre$tag)
dat.all.tre$exclosure <- factor(dat.all.tre$exclosure)
## dat.all.shr
dat.all.shr$site <- factor(dat.all.shr$site)
dat.all.shr$quadrat <- factor(dat.all.shr$quadrat)
dat.all.shr$tag <- factor(dat.all.shr$tag)
dat.all.shr$exclosure <- factor(dat.all.shr$exclosure)


dat$quad1 <-  as.factor(dat$site):dat$quadrat
dat.all$quad1 <-  as.factor(dat.all$site):dat.all$quadrat
dat.all.tre$quad1 <-  as.factor(dat.all.tre$site):dat.all.tre$quadrat
dat.all.shr$quad1 <-  as.factor(dat.all.shr$site):dat.all.shr$quadrat

# remove the NAs
# this data including only existing seedlings in each census
# height15b dataset
dat.15b <- subset(dat,dat$height15b!="NA")
summary(dat.15b)
# height16a dataset
dat.16a <- subset(dat,dat$height16a !=0 & !is.na(dat$height16a))
summary(dat.16a)
# height16b dataset
dat.16b <- subset(dat,dat$height16b !=0 & !is.na(dat$height16b))

## dat.all
# height15b dataset
dat.all.15b <- subset(dat.all,dat.all$height15b!="NA")
summary(dat.all.15b)
# height16a dataset
dat.all.16a <- subset(dat.all,dat.all$height16a !=0 & !is.na(dat.all$height16a))
summary(dat.all.16a)
# height16b dataset
dat.all.16b <- subset(dat.all,dat.all$height16b !=0 & !is.na(dat.all$height16b))

## dat.all.tre
# height15b dataset
dat.all.tre.15b <- subset(dat.all.tre,dat.all.tre$height15b!="NA")
summary(dat.all.tre.15b)
# height16a dataset
dat.all.tre.16a <- subset(dat.all.tre,dat.all.tre$height16a !=0 & !is.na(dat.all.tre$height16a))
summary(dat.all.tre.16a)
# height16b dataset
dat.all.tre.16b <- subset(dat.all.tre,dat.all.tre$height16b !=0 & !is.na(dat.all.tre$height16b))

## dat.all.shr
# height15b dataset
dat.all.shr.15b <- subset(dat.all.shr,dat.all.shr$height15b!="NA")
summary(dat.all.shr.15b)
# height16a dataset
dat.all.shr.16a <- subset(dat.all.shr,dat.all.shr$height16a !=0 & !is.na(dat.all.shr$height16a))
summary(dat.all.shr.16a)
# height16b dataset
dat.all.shr.16b <- subset(dat.all.shr,dat.all.shr$height16b !=0 & !is.na(dat.all.shr$height16b))

## gathered.dat.all
# height15b dataset
gathered.dat.all.15b <- subset(gathered.dat.all,gathered.dat.all$height15b!="NA")
summary(gathered.dat.all.15b)
# height16a dataset
gathered.dat.all.16a <- subset(gathered.dat.all,gathered.dat.all$height16a !=0 & !is.na(gathered.dat.all$height16a))
summary(dat.all.shr.16a)
# height16b dataset
gathered.dat.all.16b <- subset(gathered.dat.all,gathered.dat.all$height16b !=0 & !is.na(gathered.dat.all$height16b))


#############################################
### calculate the mean and standard error ###
## code source:http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
##---------------------------------------##

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.2) # move them .1 to the left and righ

################
### survival ###
################
## observed mean and se

## transform to present or absent value
# for dat.15b dataset, survival for 15b existing seedlings
# initial status for the 15b
dat.15b$surv15b <- as.numeric(!is.na(dat.15b$height15b))
# survival status in the 16a census
dat.15b$surv16a <- as.numeric(dat.15b$height16a !=0 & !is.na(dat.15b$height16a))
# survival status in the 16b census
dat.15b$surv16b <- as.numeric(dat.15b$height16b !=0 & !is.na(dat.15b$height16b))

# for dat.16a dataset, survival for 16a existing seedlings
# initial status for the 16a
dat.16a$surv16a <- as.numeric(dat.16a$height16a !=0 & !is.na(dat.16a$height16a))
# survival status in the 16b census
dat.16a$surv16b <- as.numeric(dat.16a$height16b !=0 & !is.na(dat.16a$height16b))

# average the survival in every quadrat
# dat.15b dataset
dat.15b.surv <-  aggregate(cbind(surv16a, surv16b) ~ treatment + exclosure + site + quadrat, 
                     data=dat.15b, FUN=mean)
# dat.16a dataset
dat.16a.surv <-  aggregate(surv16b ~ treatment + exclosure + site + quadrat, 
                           data=dat.16a, FUN=mean)

ggplot(dat.15b.surv, aes(x=treatment, y=surv16a, colour=as.factor(exclosure))) + 
  geom_point() + theme_bw() +
  geom_jitter(position = position_jitter(width = .4)) +
  facet_wrap("site")

ggplot(dat.16a.surv, aes(x=treatment, y=surv16b, colour=as.factor(exclosure))) + 
  geom_point() + theme_bw() +
  geom_jitter(position = position_jitter(width = .4)) +
  facet_wrap("site")

### the null model
Nmod1 <- glmer(surv16a ~ 1 + (1|site), data=dat.15b, family = binomial)
summary(Nmod1)

Nmod2 <- glmer(surv16a ~ 1 + (1|site/quadrat), data=dat.15b, family = binomial)
summary(Nmod2)

Nmod3 <- glmer(surv16a ~ 1 + (1|species) + (1|site/quadrat), 
               data=dat.15b, family = binomial)
summary(Nmod3)

Nmod4 <- glmer(surv16b ~ 1 + (1|species) + (1|site/quadrat), 
               data=dat.16a, family = binomial)
summary(Nmod4)

## model1
m.surv16a1 <- glmer(surv16a ~ treatment*exclosure + (1|species) + (1|site/quadrat), 
                   data=dat.15b, family = binomial)
summary(m.surv16a1)

# add growth form as fitted-effect model
m.surv16a2 <- glmer(surv16a ~ treatment*exclosure + growth.form + 
                      (1|species) + (1|site/quadrat), 
                    data=dat.15b, family = binomial)
summary(m.surv16a2)
anova(m.surv16a1,m.surv16a2, test = "F")
# dat.16a dataset
# with growthform
m.surv16b <- glmer(surv16b ~ treatment*exclosure + growth.form + 
                      (1|species) + (1|site/quadrat), 
                    data=dat.16a, family = binomial(link=logit))
summary(m.surv16b)

### dat.all
# for dat.all.15b dataset, survival for 15b existing seedlings
# initial status for the 15b
dat.all.15b$surv15b <- as.numeric(!is.na(dat.all.15b$height15b))
# survival status in the 16a census
dat.all.15b$surv16a <- as.numeric(dat.all.15b$height16a !=0 & !is.na(dat.all.15b$height16a))
# survival status in the 16b census
dat.all.15b$surv16b <- as.numeric(dat.all.15b$height16b !=0 & !is.na(dat.all.15b$height16b))

# for dat..all.16a dataset, survival for 16a existing seedlings
# initial status for the 16a
dat.all.16a$surv16a <- as.numeric(dat.all.16a$height16a !=0 & !is.na(dat.all.16a$height16a))
# survival status in the 16b census
dat.all.16a$surv16b <- as.numeric(dat.all.16a$height16b !=0 & !is.na(dat.all.16a$height16b))

### dat.all.tre
# for dat.all.tre.15b dataset, survival for 15b existing seedlings
# initial status for the 15b
dat.all.tre.15b$surv15b <- as.numeric(!is.na(dat.all.tre.15b$height15b))
# survival status in the 16a census
dat.all.tre.15b$surv16a <- as.numeric(dat.all.tre.15b$height16a !=0 & !is.na(dat.all.tre.15b$height16a))
# survival status in the 16b census
dat.all.tre.15b$surv16b <- as.numeric(dat.all.tre.15b$height16b !=0 & !is.na(dat.all.tre.15b$height16b))

# for dat.16a dataset, survival for 16a existing seedlings
# initial status for the 16a
dat.all.tre.16a$surv16a <- as.numeric(dat.all.tre.16a$height16a !=0 & !is.na(dat.all.tre.16a$height16a))
# survival status in the 16b census
dat.all.tre.16a$surv16b <- as.numeric(dat.all.tre.16a$height16b !=0 & !is.na(dat.all.tre.16a$height16b))

### dat.all.shr
# for dat.all.shr.15b dataset, survival for 15b existing seedlings
# initial status for the 15b
dat.all.shr.15b$surv15b <- as.numeric(!is.na(dat.all.shr.15b$height15b))
# survival status in the 16a census
dat.all.shr.15b$surv16a <- as.numeric(dat.all.shr.15b$height16a !=0 & !is.na(dat.all.shr.15b$height16a))
# survival status in the 16b census
dat.all.shr.15b$surv16b <- as.numeric(dat.all.shr.15b$height16b !=0 & !is.na(dat.all.shr.15b$height16b))

# for dat.16a dataset, survival for 16a existing seedlings
# initial status for the 16a
dat.all.shr.16a$surv16a <- as.numeric(dat.all.shr.16a$height16a !=0 & !is.na(dat.all.shr.16a$height16a))
# survival status in the 16b census
dat.all.shr.16a$surv16b <- as.numeric(dat.all.shr.16a$height16b !=0 & !is.na(dat.all.shr.16a$height16b))

### gathered.dat.all data
# for gathered.dat.all.15b dataset, survival for 15b existing seedlings
# initial status for the 15b
gathered.dat.all.15b$surv15b <- as.numeric(!is.na(gathered.dat.all.15b$height15b))
# survival status in the 16a census
gathered.dat.all.15b$surv16a <- as.numeric(gathered.dat.all.15b$height16a !=0 & !is.na(gathered.dat.all.15b$height16a))
# survival status in the 16b census
gathered.dat.all.15b$surv16b <- as.numeric(gathered.dat.all.15b$height16b !=0 & !is.na(gathered.dat.all.15b$height16b))

# for dat.16a dataset, survival for 16a existing seedlings
# initial status for the 16a
gathered.dat.all.16a$surv16a <- as.numeric(gathered.dat.all.16a$height16a !=0 & !is.na(gathered.dat.all.16a$height16a))
# survival status in the 16b census
gathered.dat.all.16a$surv16b <- as.numeric(gathered.dat.all.16a$height16b !=0 & !is.na(gathered.dat.all.16a$height16b))



## model1
# dat.all.15b dataset
m.all.surv16a1 <- glmer(surv16a ~ treatment*exclosure + (1|species) + (1|site/quadrat), 
                    data=dat.all.15b, family = binomial)
summary(m.all.surv16a1)

# model2: add growth form as fitted-effect model
m.all.surv16a <- glmer(surv16a ~ treatment*exclosure + growth.form + 
                      (1|species) + (1|site/quadrat), 
                    data=dat.all.15b, family = binomial)
summary(m.all.surv16a)

anova(m.all.surv16a1,m.all.surv16a, test = "F")
# model2 is better

# dat.all.16a dataset
m.all.surv16b <- glmer(surv16b ~ treatment*exclosure + growth.form + 
                          (1|species) + (1|site/quadrat), 
                        data=dat.all.16a, family = binomial)

## dat.all.tre.15b dataset
m.all.tre.surv16a <- glmer(surv16a ~ treatment*exclosure + (1|species) + (1|site/quadrat), 
                        data=dat.all.tre.15b, family = binomial)
summary(m.all.tre.surv16a)
## dat.all.tre.16a dataset
m.all.tre.surv16b <- glmer(surv16b ~ treatment*exclosure + (1|species) + (1|site/quadrat), 
                           data=dat.all.tre.16a, family = binomial)
summary(m.all.tre.surv16b)

## dat.all.shr.15b dataset
m.all.shr.surv16a <- glmer(surv16a ~ treatment*exclosure + (1|species) + (1|site/quadrat), 
                           data=dat.all.shr.15b, family = binomial)
summary(m.all.shr.surv16a)
## dat.all.tre.16a dataset
m.all.shr.surv16b <- glmer(surv16b ~ treatment*exclosure + (1|species) + (1|site/quadrat), 
                           data=dat.all.shr.16a, family = binomial)
summary(m.all.shr.surv16b)


###---------------------------------###
## find outliers in details
# dat.15b dataset
# testdat1 without growth.from
testdat1 <- data.frame(res = resid(m.surv16a1), fit = fitted(m.surv16a1), 
                      trt = dat.15b$treatment, excl = dat.15b$exclosure,
                      site = dat.15b$site, quad = dat.15b$quad1,
                      grfo = dat.15b$growth.form, spec = dat.15b$species)   
# testdat2 with growth.from
testdat2 <- data.frame(res = resid(m.surv16a2), fit = fitted(m.surv16a2), 
                       trt = dat.15b$treatment, excl = dat.15b$exclosure,
                       site = dat.15b$site, quad = dat.15b$quad1,
                       grfo = dat.15b$growth.form, spec = dat.15b$species)   

# dat.16a dataset
testdat16a <- data.frame(res = resid(m.surv16b), fit = fitted(m.surv16b), 
                       trt = dat.16a$treatment, excl = dat.16a$exclosure,
                       site = dat.16a$site, quad = dat.16a$quad1,
                       grfo = dat.16a$growth.form, spec = dat.16a$species) 

###
qplot(x=fit, y=res, geom='smooth', data=testdat1) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual with whole data (dat.15b)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat1, facets=~trt) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each treatment (dat.15b)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat1, facets=~site) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each site (dat.15b)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat1, facets=~excl) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual between exclosure (dat.15b)")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat1, facets=~grfo) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual between growth.form (dat.15b)")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat1, facets=~quad) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each quadrat (dat.15b)")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat1, facets=~spec) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each species (dat.15b)")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat1, 
      colour=trt, facets=~site) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each site (dat.15b)")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat1, 
      colour=trt, facets=~excl) + theme_bw() +
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in exclosure (dat.15b)")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat1, 
      colour=trt, facets=~grfo) + theme_bw() +
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in growth.form (dat.15b)")+
  geom_hline(aes(yintercept=0)) + geom_point()

#####
# dat.15b
# with growth.form
qplot(x=fit, y=res, geom='smooth', data=testdat2) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual with whole data (dat.15b)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat16a) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual with whole data (dat.16a)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat2, facets=~trt) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each treatment (dat.15b)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat16a, facets=~trt) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each treatment (dat.16a)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat2, facets=~site) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each site (dat.15b)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat16a, facets=~site) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each site (dat.16a)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat2, facets=~excl) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual between exclosure (dat.15b)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat16a, facets=~excl) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual between exclosure (dat.16a)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat2, facets=~grfo) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual between growth.form (dat.15b)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat16a, facets=~grfo) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual between growth.form (dat.16a)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat2, facets=~quad) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each quadrat (dat.15b)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()
# not test this

qplot(x=fit, y=res, geom='smooth', data=testdat2, facets=~spec) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each species(dat.15b)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat16a, facets=~spec) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each species (dat.16a)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat2, 
      colour=trt, facets=~site) + theme_bw() +
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each site (dat.15b)")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat16a, 
      colour=trt, facets=~site) + theme_bw() +
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each site (dat.16a)")+
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat2, 
      colour=trt, facets=~excl) + theme_bw() +
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual between exclosure (dat.15b)")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat16a, 
      colour=trt, facets=~excl) + theme_bw() +
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual between exclosure (dat.16a)")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat2, 
      colour=trt, facets=~grfo) + theme_bw() +
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual between growth.form (dat.15b)")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat16a, 
      colour=trt, facets=~grfo) + theme_bw() +
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual between growth.form (dat.16a)")+
  geom_hline(aes(yintercept=0)) + geom_point()


### dat.all dataset
# model with growth.from
# dat.all.15b dataset
testdat.all.15b <- data.frame(res = resid(m.all.surv16a2), fit = fitted(m.all.surv16a2), 
                       trt = dat.all.15b$treatment, excl = dat.all.15b$exclosure,
                       site = dat.all.15b$site, quad = dat.all.15b$quad1,
                       grfo = dat.all.15b$growth.form, spec = dat.all.15b$species)   

# dat.all.16a dataset
testdat.all.16a <- data.frame(res = resid(m.all.surv16b), fit = fitted(m.all.surv16b), 
                         trt = dat.all.16a$treatment, excl = dat.all.16a$exclosure,
                         site = dat.all.16a$site, quad = dat.all.16a$quad1,
                         grfo = dat.all.16a$growth.form, spec = dat.all.16a$species) 

### plot the fitted value vs. residuals
# with growth.form
qplot(x=fit, y=res, geom='smooth', data=testdat.all.15b) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in whole data (dat.all.15b)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.16a) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in whole data (dat.all.16a)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.15b, facets=~trt) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each treatment (dat.all.15b)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.16a, facets=~trt) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each treatment (dat.all.16a)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.15b, facets=~site) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each site (dat.all.15b)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.16a, facets=~site) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each site (dat.all.16a)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.15b, facets=~excl) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual between exclosure (dat.all.15b)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.16a, facets=~excl) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual between exclosure (dat.all.16a)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.15b, facets=~grfo) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual between growth.form (dat.all.15b)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.16a, facets=~grfo) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual between growth.form (dat.all.16a)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.15b, facets=~spec) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each species (dat.all.15b)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.16a, facets=~spec) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each species (dat.all.16a)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.15b, 
      colour=trt, facets=~site) + theme_bw() +
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each site (dat.all.15b)")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.16a, 
      colour=trt, facets=~site) + theme_bw() +
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each site (dat.all.16a)")+
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.15b, 
      colour=trt, facets=~excl) + theme_bw() +
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual between exclosure (dat.all.15b)")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.16a, 
      colour=trt, facets=~excl) + theme_bw() +
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual between exclosure (dat.all.16a)")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.15b, 
      colour=trt, facets=~grfo) + theme_bw() +
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual between growth.form (dat.all.15b)")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.16a, 
      colour=trt, facets=~grfo) + theme_bw() +
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual between growth.form (dat.all.16a)")+
  geom_hline(aes(yintercept=0)) + geom_point()


#### dat.all.tre dataset
# dat.all.tre.15b dataset
testdat.all.tre.15b <- data.frame(res = resid(m.all.tre.surv16a), fit = fitted(m.all.tre.surv16a), 
                              trt = dat.all.tre.15b$treatment, excl = dat.all.tre.15b$exclosure,
                              site = dat.all.tre.15b$site, quad = dat.all.tre.15b$quad1,
                              spec = dat.all.tre.15b$species)   

# dat.all.tre.16a dataset
testdat.all.tre.16a <- data.frame(res = resid(m.all.tre.surv16b), fit = fitted(m.all.tre.surv16b), 
                              trt = dat.all.tre.16a$treatment, excl = dat.all.tre.16a$exclosure,
                              site = dat.all.tre.16a$site, quad = dat.all.tre.16a$quad1,
                              spec = dat.all.tre.16a$species) 

#### dat.all.shr dataset
# dat.all.shr.15b dataset
testdat.all.shr.15b <- data.frame(res = resid(m.all.shr.surv16a), fit = fitted(m.all.shr.surv16a), 
                                  trt = dat.all.shr.15b$treatment, excl = dat.all.shr.15b$exclosure,
                                  site = dat.all.shr.15b$site, quad = dat.all.shr.15b$quad1,
                                  spec = dat.all.shr.15b$species)   

# dat.all.shr.16a dataset
testdat.all.shr.16a <- data.frame(res = resid(m.all.shr.surv16b), fit = fitted(m.all.shr.surv16b), 
                                  trt = dat.all.shr.16a$treatment, excl = dat.all.shr.16a$exclosure,
                                  site = dat.all.shr.16a$site, quad = dat.all.shr.16a$quad1,
                                  spec = dat.all.shr.16a$species) 

### plot the fitted value vs. residuals
## dat.all.tre dataset
qplot(x=fit, y=res, geom='smooth', data=testdat.all.tre.15b) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("dat.all.tre.15b")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.tre.16a) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in whole data")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.tre.15b, facets=~trt) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each treatment")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.tre.16a, facets=~trt) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each treatment")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.tre.15b, facets=~site) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each site")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.tre.16a, facets=~site) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each site")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.tre.15b, facets=~excl) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual between exclosure")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.tre.16a, facets=~excl) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual between exclosure")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.tre.15b, facets=~spec) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each species")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.tre.16a, facets=~spec) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each species")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.tre.15b, 
      colour=trt, facets=~site) + theme_bw() +
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each site")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.tre.16a, 
      colour=trt, facets=~site) + theme_bw() +
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each site")+
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.tre.15b, 
      colour=trt, facets=~excl) + theme_bw() +
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual between exclosure")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.tre.16a, 
      colour=trt, facets=~excl) + theme_bw() +
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual between exclosure")+
  geom_hline(aes(yintercept=0)) + geom_point()

##----------------------------------##
## dat.all.shr dataset
qplot(x=fit, y=res, geom='smooth', data=testdat.all.shr.15b) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in whole data")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.shr.16a) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in whole data")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.shr.15b, facets=~trt) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each treatment")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.shr.16a, facets=~trt) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each treatment")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.shr.15b, facets=~site) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each site")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.shr.16a, facets=~site) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each site")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.shr.15b, facets=~excl) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual between exclosure")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.shr.16a, facets=~excl) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual between exclosure")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.shr.15b, facets=~spec) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each species")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.shr.16a, facets=~spec) + 
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each species")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.shr.15b, 
      colour=trt, facets=~site) + theme_bw() +
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual in each site")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.shr.16a, 
      colour=trt, facets=~site) + theme_bw() +
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual in each site")+
  geom_hline(aes(yintercept=0)) + geom_point()


qplot(x=fit, y=res, geom='smooth', data=testdat.all.shr.15b, 
      colour=trt, facets=~excl) + theme_bw() +
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("fitted vs. residual between exclosure")+
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=testdat.all.shr.16a, 
      colour=trt, facets=~excl) + theme_bw() +
  xlab("fitted")+
  ylab("resi16b")+
  ggtitle("fitted vs. residual between exclosure")+
  geom_hline(aes(yintercept=0)) + geom_point()


###########################
## get the BLUPs associated with the various grouping factors
# observed mean value of for each quadrat
surv16a.mean <- aggregate(surv16a ~ site + quadrat,
                          data=dat.15b, FUN=mean)

re.surv16a1 <- ranef(m.surv16a1)
str(re.surv16a1)
re.surv16a1.quad <- re.surv16a1$quad

surv16a1.quad <- cbind(re.surv16a1.quad, surv16a.mean)
# unite site and quadrat number together
surv16a1.quad <- unite(surv16a1.quad, quad, c(site, quadrat), sep = "")

surv16a1.quad$re.quad <- surv16a1.quad$`(Intercept)`
surv16a1.quad$surv16a1.quad <- surv16a1.quad$surv16a

ggplot(surv16a1.quad, aes(x=re.quad, y=surv16a1.quad,label = quad)) + 
  geom_point() + theme_bw() +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = quad),hjust=0, vjust=0)

### model2 with growth form

re.surv16a2 <- ranef(m.surv16a2)
str(re.surv16a2)
re.surv16a2.quad <- re.surv16a2$quad

surv16a2.quad <- cbind(re.surv16a2.quad, surv16a.mean)
# unite site and quadrat number together
surv16a2.quad <- unite(surv16a2.quad, quad, c(site, quadrat), sep = "")

surv16a2.quad$re.quad <- surv16a2.quad$`(Intercept)`
surv16a2.quad$surv16a2.quad <- surv16a2.quad$surv16a

ggplot(surv16a2.quad, aes(x=re.quad, y=surv16a2.quad,label = quad)) + 
  geom_point() + theme_bw() +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = quad), size = 3, angle = 45, hjust=0.1, vjust=0.2)

##
p <- ggplot(surv16a2.quad, aes(x=re.quad, y=surv16a2.quad,label = quad)) + 
  geom_point() + theme_bw() +
  geom_smooth(method = "lm") 

## try to solve the overlapping text problem
library("ggrepel")
p +   geom_text(aes(label = quad), size = 3, 
                angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

p + geom_point() + theme_bw() +
  ggrepel::geom_text_repel(aes(label = quad), color = "black", size = 2.5, segment.color = "grey")
# not work

###------------------###
### dat.all dataset
# observed mean value of for each quadrat
all.surv16a.mean <- aggregate(surv16a ~ site + quadrat,
                          data=dat.all.15b, FUN=mean)
all.surv16b.mean <- aggregate(surv16b ~ site + quadrat,
                              data=dat.all.16a, FUN=mean)
## dat.all.15b
re.all.surv16a <- ranef(m.all.surv16a)
str(re.all.surv16a)
re.all.surv16a.quad <- re.all.surv16a$quad

all.surv16a.quad <- cbind(re.all.surv16a.quad, all.surv16a.mean)
# unite site and quadrat number together
all.surv16a.quad <- unite(all.surv16a.quad, quad, c(site, quadrat), sep = "")

all.surv16a.quad$re.quad <- all.surv16a.quad$`(Intercept)`
all.surv16a.quad$surv16a.quad <- all.surv16a.quad$surv16a

ggplot(all.surv16a.quad, aes(x=re.quad, y=surv16a.quad,label = quad)) + 
  geom_point() + theme_bw() +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = quad), size = 3, angle = 45, hjust=0.1, vjust=0.2) + 
  ggtitle("dat.all")

## dat.all.16a
re.all.surv16b <- ranef(m.all.surv16b)
str(re.all.surv16b)
re.all.surv16b.quad <- re.all.surv16b$quad

all.surv16b.quad <- cbind(re.all.surv16b.quad, all.surv16b.mean)
# unite site and quadrat number together
all.surv16b.quad <- unite(all.surv16b.quad, quad, c(site, quadrat), sep = "")

all.surv16b.quad$re.quad <- all.surv16b.quad$`(Intercept)`
all.surv16b.quad$surv16b.quad <- all.surv16b.quad$surv16b

ggplot(all.surv16b.quad, aes(x=re.quad, y=surv16b.quad,label = quad)) + 
  geom_point() + theme_bw() +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = quad), size = 3, angle = 45, hjust=0.1, vjust=0.2) +
  ggtitle("dat.all")


### dat.all.tre dataset
# observed mean value of for each quadrat
all.tre.surv16a.mean <- aggregate(surv16a ~ site + quadrat,
                              data=dat.all.tre.15b, FUN=mean)
all.tre.surv16b.mean <- aggregate(surv16b ~ site + quadrat,
                              data=dat.all.tre.16a, FUN=mean)

## dat.all.tre.15b
re.all.tre.surv16a <- ranef(m.all.tre.surv16a)
str(re.all.tre.surv16a)
re.all.tre.surv16a.quad <- re.all.tre.surv16a$quad

all.tre.surv16a.quad <- cbind(re.all.tre.surv16a.quad, all.tre.surv16a.mean)
# unite site and quadrat number together
all.tre.surv16a.quad <- unite(all.tre.surv16a.quad, quad, c(site, quadrat), sep = "")

all.tre.surv16a.quad$re.quad <- all.tre.surv16a.quad$`(Intercept)`
all.tre.surv16a.quad$surv16a.quad <- all.tre.surv16a.quad$surv16a

ggplot(all.tre.surv16a.quad, aes(x=re.quad, y=surv16a.quad,label = quad)) + 
  geom_point() + theme_bw() +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = quad), size = 3, angle = 45, hjust=0.1, vjust=0.2) +
  ggtitle("dat.all.tre")

## dat.all.tre.16a
re.all.tre.surv16b <- ranef(m.all.tre.surv16b)
str(re.all.tre.surv16b)
re.all.tre.surv16b.quad <- re.all.tre.surv16b$quad

all.tre.surv16b.quad <- cbind(re.all.tre.surv16b.quad, all.tre.surv16b.mean)
# unite site and quadrat number together
all.tre.surv16b.quad <- unite(all.tre.surv16b.quad, quad, c(site, quadrat), sep = "")

all.tre.surv16b.quad$re.quad <- all.tre.surv16b.quad$`(Intercept)`
all.tre.surv16b.quad$surv16b.quad <- all.tre.surv16b.quad$surv16b

ggplot(all.tre.surv16b.quad, aes(x=re.quad, y=surv16b.quad,label = quad)) + 
  geom_point() + theme_bw() +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = quad), size = 3, angle = 45, hjust=0.1, vjust=0.2) +
  ggtitle("dat.all.tre")

### dat.all.shr dataset
# observed mean value of for each quadrat
all.shr.surv16a.mean <- aggregate(surv16a ~ site + quadrat,
                                  data=dat.all.shr.15b, FUN=mean)
all.shr.surv16b.mean <- aggregate(surv16b ~ site + quadrat,
                                  data=dat.all.shr.16a, FUN=mean)

## dat.all.shr.15b
re.all.shr.surv16a <- ranef(m.all.shr.surv16a)
str(re.all.shr.surv16a)
re.all.shr.surv16a.quad <- re.all.shr.surv16a$quad

all.shr.surv16a.quad <- cbind(re.all.shr.surv16a.quad, all.shr.surv16a.mean)
# unite site and quadrat number together
all.shr.surv16a.quad <- unite(all.shr.surv16a.quad, quad, c(site, quadrat), sep = "")

all.shr.surv16a.quad$re.quad <- all.shr.surv16a.quad$`(Intercept)`
all.shr.surv16a.quad$surv16a.quad <- all.shr.surv16a.quad$surv16a

ggplot(all.shr.surv16a.quad, aes(x=re.quad, y=surv16a.quad,label = quad)) + 
  geom_point() + theme_bw() +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = quad), size = 3, angle = 45, hjust=0.1, vjust=0.2) +
  ggtitle("dat.all.shr")


## dat.all.tre.16a
re.all.shr.surv16b <- ranef(m.all.shr.surv16b)
str(re.all.shr.surv16b)
re.all.shr.surv16b.quad <- re.all.shr.surv16b$quad

all.shr.surv16b.quad <- cbind(re.all.shr.surv16b.quad, all.shr.surv16b.mean)
# unite site and quadrat number together
all.shr.surv16b.quad <- unite(all.shr.surv16b.quad, quad, c(site, quadrat), sep = "")

all.shr.surv16b.quad$re.quad <- all.shr.surv16b.quad$`(Intercept)`
all.shr.surv16b.quad$surv16b.quad <- all.shr.surv16b.quad$surv16b

ggplot(all.shr.surv16b.quad, aes(x=re.quad, y=surv16b.quad,label = quad)) + 
  geom_point() + theme_bw() +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = quad), size = 3, angle = 45, hjust=0.1, vjust=0.2) +
  ggtitle("dat.all.shr")


## the AUC (Area Under the Curve of the Receiver Operator Characteristic
## If it is 0.5 it is no better than random, 0.5 < AUC < 0.7 is not good, 0.7 < AUC < 0.8 is useful, 
## 0.8 < AUC < 0.9 is good 0.9 < AUC < 1.0 is very good
# install.packages("verification")
library(verification)
roc.area(obs=getME(m.surv16a2, 'y'), pred=fitted(m.surv16a2))
roc.area(obs=getME(m.surv16b, 'y'), pred=fitted(m.surv16b))
roc.area(obs=getME(m.all.surv16a, 'y'), pred=fitted(m.all.surv16a))
roc.area(obs=getME(m.all.surv16b, 'y'), pred=fitted(m.all.surv16b))
roc.area(obs=getME(m.all.tre.surv16a, 'y'), pred=fitted(m.all.tre.surv16a))
roc.area(obs=getME(m.all.tre.surv16b, 'y'), pred=fitted(m.all.tre.surv16b))
roc.area(obs=getME(m.all.shr.surv16a, 'y'), pred=fitted(m.all.shr.surv16a))
roc.area(obs=getME(m.all.shr.surv16b, 'y'), pred=fitted(m.all.shr.surv16b))

#####
##  bootstrapping to get pval
summary(m.surv16a2)
anova(m.surv16a2)

m.surv16a0 <- update(m.surv16a2, ~.-treatment*exclosure-growth.form)
b1 <- PBmodcomp(m.surv16a2, m.surv16a0, nsim=99)

## prediction
newdata <- data.frame(x=dat.15b$treatment*(dat.15b$exclosure)+dat.15b$growth.from)

p1 <- predict(m.surv16a2,newdata = dat.15b, type="response")
p1

p1 <- predict(m.surv16a2,newdata = dat.15b, re.form=~0)
p1

plogis (predict(m.surv16a2,newdata = dat.15b, re.form=~0))


####---------------------------####
## JAN 12,2017
## consider snow cover treatment as random effect effects (slope)


