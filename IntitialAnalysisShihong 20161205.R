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



####################
## Shihong play
####################

library(ggplot2)

dat <- read.csv("Changbaishan/data/Woodyseedlings 20161205.csv")
summary(dat)

dat$tag <- factor(dat$tag)
dat$exclosure <- factor(dat$exclosure)

quad <-  as.factor(dat$code):dat$sub.code
dat <- cbind (dat,quad)

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


######
# calculate the mean and standard error
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

## plot the figure
# Standard error of the mean
<<<<<<< Updated upstream

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.2) # move them .1 to the left and righ


################
### survival ###
################

# remove the NAs
# height15b dataset
summary(dat.15b)
# height16a dataset
summary(dat.16a)

=======

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.2) # move them .1 to the left and righ


################
### survival ###
################

# remove the NAs
# height15b dataset
summary(dat.15b)
# height16a dataset
summary(dat.16a)

>>>>>>> Stashed changes
## For the 15b dataset
# initial status for the 15b
dat.15b$surv0 <- as.numeric(!is.na(dat.15b$height15b))
# survival status in the 16a census
dat.15b$surv1 <- as.numeric(dat.15b$height16a !=0 & !is.na(dat.15b$height16a))
# survival status in the 16b census
dat.15b$surv2 <- as.numeric(dat.15b$height16b !=0 & !is.na(dat.15b$height16b))

# get mean survival rate for the quadarts for each treatment in different exclosure treatments
# within each of the three sites
dat.15b.surv  <- aggregate(cbind(surv0,surv1,surv2)~ treatment + exclosure + site + sub.code, data=dat.15b, FUN=sum)
#View(dat.15b.surv)
dat.15b.surv$s1 <- dat.15b.surv$surv1/dat.15b.surv$surv0
dat.15b.surv$s2 <- dat.15b.surv$surv2/dat.15b.surv$surv0

# survival rate in the 16a census
s1.cal <- summarySE(dat.15b.surv,measurevar="s1", groupvars=c("site","exclosure","treatment"))
#View(s1.cal)

ggplot(s1.cal, aes(x=treatment, y=s1, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=s1-se, ymax=s1+se), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")

# survival rate in the 16b census
summary(dat.15b.surv$s2)
<<<<<<< Updated upstream
which.max(dat.15b.surv$s2)
=======
>>>>>>> Stashed changes
s2.cal <- summarySE(dat.15b.surv,measurevar="s2", groupvars=c("site","exclosure","treatment"))
#View(s2.cal)
summary(s2.cal)

ggplot(s2.cal, aes(x=treatment, y=s2, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=s2-se, ymax=s2+se), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")

#################
### abundance ###
#################

dat$ab1 <- as.numeric(dat$height16a !=0 & !is.na(dat$height16a))
dat$ab2 <- as.numeric(dat$height16b !=0 & !is.na(dat$height16b))

dat.ab  <- aggregate(cbind(ab1, ab2) ~ treatment + site +exclosure + sub.code, data=dat, FUN=sum)
#View(dat.ab)

#dat.ab  <- aggregate(cbind(ab1, ab2) ~ treatment + site + exclosure + sub.code, data=dat, FUN=sum)

# the violin plot
# abundance in the 16a census
ab1 <- ggplot(data=dat.ab, aes(x=treatment, y=ab1,colour=as.factor(exclosure)))
ab1 <- ab1 + geom_violin()
ab1 + facet_wrap( "site" )

# abundance in the 16b census
ab2 <- ggplot(data=dat.ab, aes(x=treatment, y=ab2,colour=as.factor(exclosure)))
ab2 <- ab2 + geom_violin()
ab2 + facet_wrap( "site" )


# abundance in the 16a census
ab1.cal <- summarySE(dat.ab,measurevar="ab1", groupvars=c("site","exclosure","treatment"))
#View(ab1.cal)
# abundance in the 16b census
ab2.cal <- summarySE(dat.ab,measurevar="ab2", groupvars=c("site","exclosure","treatment"))
#View(ab2.cal)

# 16a
ggplot(ab1.cal, aes(x=treatment, y=ab1, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=ab1-se, ymax=ab1+se), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")

# 16b
ggplot(ab2.cal, aes(x=treatment, y=ab2, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=ab2-se, ymax=ab2+se), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")
<<<<<<< Updated upstream

# Use 95% confidence interval instead of SEM
ggplot(ab1.cal, aes(x=treatment, y=ab1, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=ab1-ci, ymax=ab1+ci), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")

### end of the abundance

#########################
### number of species ###
#########################

# 16a
dat.16a$no.species <-as.numeric(dat.16a$species)
dat.16a.spec  <- aggregate(no.species~ treatment + exclosure + site + sub.code, data=dat.16a, FUN=function(x)length(unique(x)))
#View(dat.16a.spec)
# 16b
dat.16b$no.species <-as.numeric(dat.16b$species)
dat.16b.spec  <- aggregate(no.species~ treatment + exclosure + site + sub.code, data=dat.16b, FUN=function(x)length(unique(x)))
#View(dat.16b.spec)

spec1.cal <- summarySE(dat.16a.spec,measurevar="no.species", groupvars=c("site","exclosure","treatment"))
spec2.cal <- summarySE(dat.16b.spec,measurevar="no.species", groupvars=c("site","exclosure","treatment"))

# 16a 
ggplot(spec1.cal, aes(x=treatment, y=no.species, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=no.species-se, ymax=no.species+se), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")

# 16b
ggplot(spec2.cal, aes(x=treatment, y=no.species, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=no.species-se, ymax=no.species+se), width=.1,position=pd) +
=======

# Use 95% confidence interval instead of SEM
ggplot(ab1.cal, aes(x=treatment, y=ab1, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=ab1-ci, ymax=ab1+ci), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")

### end of the abundance

#########################
### number of species ###
#########################

# 16a
dat.16a$no.species <-as.numeric(dat.16a$species)
dat.16a.spec  <- aggregate(no.species~ treatment + exclosure + site + sub.code, data=dat.16a, FUN=function(x)length(unique(x)))
#View(dat.16a.spec)
# 16b
dat.16b$no.species <-as.numeric(dat.16b$species)
dat.16b.spec  <- aggregate(no.species~ treatment + exclosure + site + sub.code, data=dat.16b, FUN=function(x)length(unique(x)))
#View(dat.16b.spec)

spec1.cal <- summarySE(dat.16a.spec,measurevar="no.species", groupvars=c("site","exclosure","treatment"))
spec2.cal <- summarySE(dat.16b.spec,measurevar="no.species", groupvars=c("site","exclosure","treatment"))

# 16a 
ggplot(spec1.cal, aes(x=treatment, y=no.species, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=no.species-se, ymax=no.species+se), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")

# 16b
ggplot(spec2.cal, aes(x=treatment, y=no.species, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=no.species-se, ymax=no.species+se), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")

### end of the number of species


###############
### grouwth ###
###############

# 15 existing seedlings data
head (dat.15b)
# growth change from 15b to 16a
dat.15b$heig1 <- dat.15b$height16a-dat.15b$height15b
# growth change from 15b to 16b
dat.15b$heig2 <- dat.15b$height16b-dat.15b$height15b

dat.15b.heig <- aggregate(cbind(heig1,heig2)~ treatment + exclosure + site + sub.code, data=dat.15b, FUN=mean)
heig1.cal <- summarySE(dat.15b.heig,measurevar="heig1", groupvars=c("site","exclosure","treatment"))
heig2.cal <- summarySE(dat.15b.heig,measurevar="heig2", groupvars=c("site","exclosure","treatment"))

# 16a
ggplot(heig1.cal, aes(x=treatment, y=heig1, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=heig1-se, ymax=heig1+se), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")

# 16b
ggplot(heig2.cal, aes(x=treatment, y=heig2, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=heig2-se, ymax=heig2+se), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")


####
# consider the effects of species
dat.15b.heig1 <- aggregate(cbind(heig1,heig2)~ treatment + exclosure + site + sub.code + species, data=dat.15b, FUN=mean)
heig1.1.cal <- summarySE(dat.15b.heig1,measurevar="heig1", groupvars=c("site","exclosure","treatment"))
heig1.2.cal <- summarySE(dat.15b.heig1,measurevar="heig2", groupvars=c("site","exclosure","treatment"))

# 16a
ggplot(heig1.1.cal, aes(x=treatment, y=heig1, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=heig1-se, ymax=heig1+se), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")

# 16b
ggplot(heig1.2.cal, aes(x=treatment, y=heig2, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=heig2-se, ymax=heig2+se), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")

########################
### using glms function
## incorporating the species effect

# get mean survival rate for the quadarts for each treatment in different exclosure treatments
# within each of the three sites
dat.15b.surv  <- aggregate(cbind(surv0,surv1,surv2)~ treatment + exclosure + site + sub.code + species, data=dat.15b, FUN=sum)

dat.15b.surv$s1 <- dat.15b.surv$surv1/dat.15b.surv$surv0
dat.15b.surv$s2 <- dat.15b.surv$surv2/dat.15b.surv$surv0
dat.15b.surv$d1 <- 1-dat.15b.surv$surv1/dat.15b.surv$surv0
dat.15b.surv$d2 <- 1-dat.15b.surv$surv2/dat.15b.surv$surv0

mod1 <-  glm(cbind(s1,d1) ~ species, data=dat.15b.surv, family = binomial)
mod1

dat.15b.surv$resi1 <- residuals(mod1)

ggplot(dat.15b.surv, aes(x=treatment, y=resi1, colour=as.factor(exclosure))) + 
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")

# survival rate in the 16a census
surv1.cal <- summarySE(dat.15b.surv,measurevar="resi1", groupvars=c("site","exclosure","treatment","species"))

ggplot(surv1.cal, aes(x=treatment, y=resi1, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=resi1-se, ymax=resi1+se), width=.1,position=pd) +
>>>>>>> Stashed changes
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")

<<<<<<< Updated upstream
### end of the number of species


###############
### grouwth ###
###############

# 15 existing seedlings data
head (dat.15b)
# growth change from 15b to 16a
dat.15b$heig1 <- dat.15b$height16a-dat.15b$height15b
# growth change from 15b to 16b
dat.15b$heig2 <- dat.15b$height16b-dat.15b$height15b

dat.15b.heig <- aggregate(cbind(heig1,heig2)~ treatment + exclosure + site + sub.code, data=dat.15b, FUN=mean)
heig1.cal <- summarySE(dat.15b.heig,measurevar="heig1", groupvars=c("site","exclosure","treatment"))
heig2.cal <- summarySE(dat.15b.heig,measurevar="heig2", groupvars=c("site","exclosure","treatment"))

# 16a
ggplot(heig1.cal, aes(x=treatment, y=heig1, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=heig1-se, ymax=heig1+se), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")

# 16b
ggplot(heig2.cal, aes(x=treatment, y=heig2, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=heig2-se, ymax=heig2+se), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")


####
# consider the effects of species
dat.15b.heig1 <- aggregate(cbind(heig1,heig2)~ treatment + exclosure + site + sub.code + species, data=dat.15b, FUN=mean)
heig1.1.cal <- summarySE(dat.15b.heig1,measurevar="heig1", groupvars=c("site","exclosure","treatment"))
heig1.2.cal <- summarySE(dat.15b.heig1,measurevar="heig2", groupvars=c("site","exclosure","treatment"))
=======
ggplot(dat.15b.surv, aes(x=treatment, y=resi1, colour=as.factor(exclosure))) + 
  geom_violin() +
  facet_wrap("site")


# how about incorporating the growth form effect
dat.15b.surv  <- aggregate(cbind(surv0,surv1,surv2)~ treatment + exclosure + site + sub.code + growth.form, data=dat.15b, FUN=sum)

dat.15b.surv$s1 <- dat.15b.surv$surv1/dat.15b.surv$surv0
dat.15b.surv$s2 <- dat.15b.surv$surv2/dat.15b.surv$surv0
dat.15b.surv$d1 <- 1-dat.15b.surv$surv1/dat.15b.surv$surv0
dat.15b.surv$d2 <- 1-dat.15b.surv$surv2/dat.15b.surv$surv0

mod1 <-  glm(cbind(s1,d1) ~ growth.form, data=dat.15b.surv, family = binomial)
mod1

dat.15b.surv$resi1 <- residuals(mod1)
ggplot(dat.15b.surv, aes(x=treatment, y=resi1, colour=as.factor(exclosure))) + 
  geom_violin() +
  facet_wrap("site")

# survival rate in the 16a census
surv1.cal <- summarySE(dat.15b.surv,measurevar="resi1", groupvars=c("site","exclosure","treatment","growth.form"))

ggplot(surv1.cal, aes(x=treatment, y=resi1, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=resi1-se, ymax=resi1+se), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")

ggplot(dat.15b.surv, aes(x=treatment, y=resi1, colour=as.factor(exclosure))) + 
  geom_violin() +
  facet_wrap("site")
>>>>>>> Stashed changes

# 16a
ggplot(heig1.1.cal, aes(x=treatment, y=heig1, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=heig1-se, ymax=heig1+se), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")

<<<<<<< Updated upstream
# 16b
ggplot(heig1.2.cal, aes(x=treatment, y=heig2, colour=as.factor(exclosure))) + 
  geom_errorbar(aes(ymin=heig2-se, ymax=heig2+se), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  facet_wrap("site")
=======
s1 <- glm (surv1 ~ species, data= dat.15b, family = binomial)
dat.15b$s1.resi <- residuals(s1)
#View(dat.15b)
>>>>>>> Stashed changes
