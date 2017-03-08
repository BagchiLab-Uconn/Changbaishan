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
## the original data
## maybe dat.all is better for more shrub sample size included
dat <- read.csv("Changbaishan/data/Woodyseedlings 20170112.csv")
####
## dat.gathered
## put all the height of each census into one column
## use tidyr to get longer data with all height in one column
gathered.dat <- gather(dat, heightcensus, height, height15b:height16b)


####
## analyze trees and shrubs seperately
gathered.dat.tre <- gathered.dat[gathered.dat$growth.form == "tree",]
gathered.dat.shr <- gathered.dat[gathered.dat$growth.form == "shrub",]

summary(dat)
summary(gathered.dat)
summary(gathered.dat.tre)
summary(gathered.dat.tre)


dat$site <- factor(dat$site)
dat$quadrat <- factor(dat$quadrat)
dat$tag <- factor(dat$tag)
dat$exclosure <- factor(dat$exclosure)
dat$snow <- factor(dat$snow)
#quad <-  as.factor(dat$code):dat$sub.code
#dat <- cbind (dat,quad)
gathered.dat$site <- factor(gathered.dat$site)
gathered.dat$quadrat <- factor(gathered.dat$quadrat)
gathered.dat$tag <- factor(gathered.dat$tag)
gathered.dat$exclosure <- factor(gathered.dat$exclosure)
gathered.dat$snow <- factor(gathered.dat$snow)
## dat.all.tre
gathered.dat.tre$site <- factor(gathered.dat.tre$site)
gathered.dat.tre$quadrat <- factor(gathered.dat.tre$quadrat)
gathered.dat.tre$tag <- factor(gathered.dat.tre$tag)
gathered.dat.tre$exclosure <- factor(gathered.dat.tre$exclosure)
gathered.dat.tre$snow <- factor(gathered.dat.tre$snow)
## dat.all.shr
gathered.dat.shr$site <- factor(gathered.dat.shr$site)
gathered.dat.shr$quadrat <- factor(gathered.dat.shr$quadrat)
gathered.dat.shr$tag <- factor(gathered.dat.shr$tag)
gathered.dat.shr$exclosure <- factor(gathered.dat.shr$exclosure)
gathered.dat.shr$snow <- factor(gathered.dat.shr$snow)

dat$quad1 <-  as.factor(dat$site):dat$quadrat
gathered.dat$quad1 <-  as.factor(gathered.dat$site):gathered.dat$quadrat
gathered.dat.tre$quad1 <-  as.factor(gathered.dat.tre$site):gathered.dat.tre$quadrat
gathered.dat.shr$quad1 <-  as.factor(gathered.dat.shr$site):gathered.dat.shr$quadrat


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

####---------------------------####
## JAN 12,2017
## consider snow cover treatment as random effect effects (slope)

### gathered.dat
# survival status 
gathered.dat$surv <- as.numeric(gathered.dat$height !=0 & !is.na(gathered.dat$height))
### gathered.dat.tre
gathered.dat.tre$surv <- as.numeric(gathered.dat.tre$height !=0 & !is.na(gathered.dat.tre$height))
### gathered.dat.shr
gathered.dat.shr$surv <- as.numeric(gathered.dat.shr$height !=0 & !is.na(gathered.dat.shr$height))


## mixed-effect model
# gathered.dat
m.g.all.surv0 <- glmer(surv ~ heightcensus + pesticide*exclosure + snow*exclosure + (1+snow|species) + (1|site/quadrat), 
                      data=gathered.dat, family = binomial)
summary(m.g.all.surv0)
m.g.all.surv <- glmer(surv ~ heightcensus + pesticide*exclosure + snow*exclosure + growth.form + 
                      (1+snow|species) + (1|site/quadrat), data=gathered.dat, family = binomial)
summary(m.g.all.surv)

# gathered.dat.tre
m.g.tre.surv <- glmer(surv ~ heightcensus + pesticide*exclosure + snow*exclosure + (1+snow|species) + (1|site/quadrat), 
                      data=gathered.dat.tre, family = binomial)
summary(m.g.tre.surv)
# gathered.dat.shr
m.g.shr.surv <- glmer(surv ~ heightcensus + pesticide*exclosure + snow*exclosure + (1+snow|species) + (1|site/quadrat), 
                      data=gathered.dat.shr, family = binomial)
summary(m.g.shr.surv)



###---------------------------------###
## find outliers in details
## gathered.dat
test.g.all <- data.frame(res = resid(m.g.all.surv), fit = fitted(m.g.all.surv), 
                      pest = gathered.dat$pesticide, snow = gathered.dat$snow,
                      excl = gathered.dat$exclosure, cens = gathered.dat$heightcensus,
                      site = gathered.dat$site, quad = gathered.dat$quad1,
                      grfo = gathered.dat$growth.form, spec = gathered.dat$species)   
## gathered.dat.tre
test.g.tre <- data.frame(res = resid(m.g.tre.surv), fit = fitted(m.g.tre.surv), 
                         pest = gathered.dat.tre$pesticide, snow = gathered.dat.tre$snow,
                         excl = gathered.dat.tre$exclosure, cens = gathered.dat.tre$heightcensus,
                         site = gathered.dat.tre$site, quad = gathered.dat.tre$quad1,
                         grfo = gathered.dat.tre$growth.form, spec = gathered.dat.tre$species) 

## gathered.dat.shr
test.g.shr <- data.frame(res = resid(m.g.shr.surv), fit = fitted(m.g.shr.surv), 
                         pest = gathered.dat.shr$pesticide, snow = gathered.dat.shr$snow,
                         excl = gathered.dat.shr$exclosure, cens = gathered.dat.shr$heightcensus,
                         site = gathered.dat.shr$site, quad = gathered.dat.shr$quad1,
                         grfo = gathered.dat.shr$growth.form, spec = gathered.dat.shr$species) 

###
qplot(x=fit, y=res, geom='smooth', data=test.g.all) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("fitted vs. residual with whole data (gathered.dat)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.all,facets=~snow) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("between snow trt (gathered.dat)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.all,facets=~pest) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("among pest trt (gathered.dat)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.all,facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("in each census (gathered.dat)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.all, colour=snow, facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("snow trt in each census (gathered.dat)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.all, colour=snow, facets=~site) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("snow trt in each site (gathered.dat)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.all, colour=snow, facets=~grfo) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("snow trt between growthform (gathered.dat)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.all, colour=snow, facets=~excl) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("snow trt between exclosure (gathered.dat)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


######
## gathered.dat.tre
qplot(x=fit, y=res, geom='smooth', data=test.g.tre) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("fitted vs. residual with whole data (gathered.dat.tre)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.tre,facets=~snow) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("between snow trt (gathered.dat.tre)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.tre,facets=~pest) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("among pest trt (gathered.dat.tre)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.tre,facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("in each census (gathered.dat.tre)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.tre, colour=snow, facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("snow trt in each census (gathered.dat.tre)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.tre, colour=snow, facets=~site) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("snow trt in each site (gathered.dat.tre)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.tre, colour=snow, facets=~excl) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("snow trt between exclosure (gathered.dat.tre)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

######
## gathered.dat.shr
qplot(x=fit, y=res, geom='smooth', data=test.g.shr) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("fitted vs. residual with whole data (gathered.dat.shr)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.shr,facets=~snow) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("between snow trt (gathered.dat.shr)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.shr,facets=~pest) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("among pest trt (gathered.dat.shr)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.shr,facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("in each census (gathered.dat.shr)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.shr, colour=snow, facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("snow trt in each census (gathered.dat.shr)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.shr, colour=snow, facets=~site) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("snow trt in each site (gathered.dat.shr)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.shr, colour=snow, facets=~excl) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("snow trt between exclosure (gathered.dat.shr)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()





### 
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




