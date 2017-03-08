## Library packages
library(ggplot2)
library(Matrix)
library(lme4)
library(permute)
library(lattice)
library(vegan)
library(tidyr)
library(plyr)
library(dplyr)
library(pbkrtest)
library(verification)


######>>>>>>>>>>>>>>>>>>>>>>>>>>#########
######  Woody seedlings 
######>>>>>>>>>>>>>>>>>>>>>>>>>>#########


#### load the data
## Data set1
## Pesticide treatments analysis
## period: from 2015 fall to 2016 fall, three censuses (two periods)

w.pest <- read.csv("Changbaishan/data/Pestwoodyseedlings 20170116.csv")

## Data set2
## Snow cover removal treatment
## Period: from 2014 fall to 2016 fall, five censuses (four periods)
## Initial analysis
w.snow <- read.csv("Changbaishan/data/Snowwoodyseedlings 20170115.csv")

#####################
## data preparation #
#####################

####  gather data
## put all the height of each census into one column
## add one column, and put all survival situation into this column

## First comes to the pest data
## use tidyr::gather() to get longer data with all height and surv in one column
gathered.w.pest1 <- gather(w.pest, census_surv, heig_surv, height_15fa:surv_16fa)
## use separate() to separate the height and surv in two columns
gathered.w.pest2 <- gathered.w.pest1 %>%
  separate(census_surv,c("item","census"))
## use spread() to get the final version with height and surv in separate columns
gathered.w.pest <- gathered.w.pest2 %>%
  spread(item,heig_surv)
## the data set (gathered.pest) is the data that we will analyze

## Similar for the snow data set
gathered.w.snow1 <- gather(w.snow, census_surv, heig_surv, height_14fa:surv_16fa)
gathered.w.snow2 <- gathered.w.snow1 %>%
  separate(census_surv,c("item","census"))
gathered.w.snow <- gathered.w.snow2 %>%
  spread(item,heig_surv)
## gathered.snow for snow treatment analysis


## pest.new
# we want to test pesticide effects among fungicide, insecticide and water control
# gathered.pest.new1 <- gather(pest.new, census_surv, heig_surv, height_15f:surv_16f)
## use separate() to separate the height and surv in two columns
# gathered.pest.new2 <- gathered.pest.new1 %>%
#  separate(census_surv,c("item","census"))
## use spread() to get the final version with height and surv in separate columns
# gathered.pest.new <- gathered.pest.new2 %>%
#  spread(item,heig_surv)

### Manipulate data
## gathered.pest
gathered.w.pest$site <- factor(gathered.w.pest$site)
gathered.w.pest$quadrat <- factor(gathered.w.pest$quadrat)
gathered.w.pest$tag <- factor(gathered.w.pest$tag)
gathered.w.pest$exclosure <- factor(gathered.w.pest$exclosure)
gathered.w.pest$census <-  as.factor(gathered.w.pest$census)
gathered.w.pest$quad1 <-  as.factor(gathered.w.pest$code):gathered.w.pest$quadrat

## gathered.snow
gathered.w.snow$site <- factor(gathered.w.snow$site)
gathered.w.snow$quadrat <- factor(gathered.w.snow$quadrat)
gathered.w.snow$tag <- factor(gathered.w.snow$tag)
gathered.w.snow$census <-  as.factor(gathered.w.snow$census)
gathered.w.snow$quad1 <-  as.factor(gathered.w.snow$site):gathered.w.snow$quadrat



###############################################################################
################################## survival ################################### 
###############################################################################

# survival status 
# gathered.dat$surv <- as.numeric(gathered.dat$height !=0 & !is.na(gathered.dat$height))
# gathered.dat$surv <- as.numeric(gathered.dat$height !=0)

# see the data structure
table(gathered.w.pest$surv)
table(gathered.w.pest$species)

table(gathered.w.snow$surv)

#### set new data set without NAs
gathered.w.pest0 <- subset(gathered.w.pest,!is.na(gathered.w.pest$surv))
gathered.w.snow0 <- subset(gathered.w.snow,!is.na(gathered.w.snow$surv))


#######################
## mixed-effect model #
#######################
## gathered.pest
# null model
m.w.pest.surv0 <- glmer(surv ~ (1|site/quadrat) + (1|species),
                         data=gathered.w.pest, family = binomial)
# mixed-effect models
m.w.pest.surv1 <- glmer(surv ~ pesticide + exclosure + census + growth.form +  (1|site/quadrat) + (1|species),
                       data=gathered.w.pest, family = binomial)

m.w.pest.surv2 <- glmer(surv ~ pesticide*exclosure + census + growth.form +  (1|site/quadrat) + (1|species),
                      data=gathered.w.pest, family = binomial)

## consider the interaction between pesticide and growth.form
m.w.pest.surv3 <- glmer(surv ~ census + pesticide*exclosure + pesticide*growth.form +  (1|site/quadrat) + (1|species),
                       data=gathered.w.pest, family = binomial)

## consider the interaction among pesticide, exclosure and growth.form
m.w.pest.surv4 <- glmer(surv ~ census + pesticide*exclosure*growth.form +  (1|site/quadrat) + (1|species),
                       data=gathered.w.pest, family = binomial)

## consider random slope for the pesticide in each species
m.w.pest.surv5 <- glmer(surv ~ pesticide*exclosure + census + growth.form +  (1|site/quadrat) + (1+pesticide|species),
                          data=gathered.w.pest, family = binomial)

## interaction across pesticide, exclosure and growth.form
m.w.pest.surv6 <- glmer(surv ~ (pesticide + exclosure + growth.form)^2  + census +(1|site/quadrat) + (1+pesticide|species),
                          data=gathered.w.pest, family = binomial)

m.w.pest.surv7 <- glmer(surv ~ pesticide*exclosure + growth.form +(1|site/quadrat) + (1+pesticide|species),
                          data=gathered.w.pest, family = binomial)

m.w.pest.surv8 <- glmer(surv ~ pesticide + growth.form +(1|site/quadrat) + (1+pesticide|species),
                          data=gathered.w.pest, family = binomial)

m.w.pest.surv9 <- glmer(surv ~ pesticide +(1|site/quadrat) + (1+pesticide|species),
                          data=gathered.w.pest, family = binomial)

## compare the models
anova (m.w.pest.surv0, m.w.pest.surv1, m.w.pest.surv2, m.w.pest.surv3, m.w.pest.surv4,
       m.w.pest.surv5, m.w.pest.surv6, m.w.pest.surv7, m.w.pest.surv8, m.w.pest.surv9)
# AIC(m.w.pest.surv0, m.w.pest.surv1, m.w.pest.surv2, m.w.pest.surv3, m.w.pest.surv4, 
#    m.w.pest.surv5, m.w.pest.surv6, m.w.pest.surv7, m.w.pest.surv8, m.w.pest.surv9)

## the m.w.pest.surv1(without interaction) is better
m.w.pest.surv <- m.w.pest.surv1

## drop one level 
m.w.pest.surv.drop <- drop1(m.w.pest.surv,test="Chi")
m.w.pest.surv.drop 


################# Model validation ######################
plot(m.w.pest.surv)

w.pest.surv.resi <- resid(m.w.pest.surv)

dev.new()
op <- par(mfrow=c(2,3),mar = c(5, 4, 1, 2))
hist(w.pest.surv.resi, xlab = "Residuals", main = "")
plot(gathered.w.pest0$pesticide, w.pest.surv.resi, xlab = "Pesticide",
     ylab = "Residuals")
abline(h=0,col="blue")
plot(gathered.w.pest0$exclosure, w.pest.surv.resi, xlab = "Exclosure",
     ylab = "Residuals")
abline(h=0,col="blue")
plot(gathered.w.pest0$census, w.pest.surv.resi, xlab = "Census",
     ylab = "Residuals")
abline(h=0,col="blue")
plot(gathered.w.pest0$growth.form, w.pest.surv.resi, xlab = "Growth form",
     ylab = "Residuals")
abline(h=0,col="blue")
plot(gathered.w.pest0$site, w.pest.surv.resi, xlab = "Site",
     ylab = "Residuals")
abline(h=0,col="blue")

par(op)

## quadrats
plot(gathered.w.pest0$quadrat, w.pest.surv.resi, xlab = "Quadrat",
     ylab = "Residuals")
abline(h=0,col="blue")




#######
# relevel factor, put W (water control) as the intercept
# m.g.all.surv <- profile(glmer(surv ~ pesticide*exclosure + census + growth.form +  (1|site/quadrat) + (1|species),
#                              within(data=gathered.dat, pesticide <- relevel(pesticide, "W")), family = binomial)) 
# summary(m.g.all.surv)

gathered.w.pest$pesticide <- relevel(gathered.w.pest$pesticide, ref = "W")
m.w.pest.surv <- glmer(surv ~ pesticide*exclosure + census + growth.form +  (1|site/quadrat) + (1|species),
                      data=gathered.w.pest, family = binomial)
summary(m.w.pest.surv)


###### weight with height

m.w.pest.surv10 <- glmer(surv ~ pesticide + exclosure + census + growth.form +  (1|site/quadrat) + (1|species),
                        data=gathered.w.pest, weights=height, family = binomial)

summary(m.w.pest.surv10)





##############################################################################################
###################################### gathered.snow #########################################
##############################################################################################
m.w.snow.surv1 <- glmer(surv ~ treatment + census + growth.form +  (1|site/quadrat) + (1|species),
                          data=gathered.w.snow, family = binomial)
# model includes interaction between treatment and growth.form
m.w.snow.surv2 <- glmer(surv ~ treatment*growth.form + census + (1|site/quadrat) + (1|species),
                         data=gathered.w.snow, family = binomial)

anova(m.w.snow.surv1, m.w.snow.surv2)
# m.w.snow.surv2 is better

# consider the random slope
m.w.snow.surv3 <- glmer(surv ~ treatment*growth.form + census + (1|site/quadrat) + (1 + treatment|species),
                          data=gathered.w.snow, family = binomial)

anova (m.w.snow.surv2, m.w.snow.surv3)

# No difference between them

m.w.snow.surv <- m.w.snow.surv2
summary(m.w.snow.surv)
###---------------------------------####

#############################################################################
######>>>>>>>>>>>>>>>> Survival Porportion <<<<<<<<<<<<<<<<<<################
#############################################################################

######>>>>>>>>>>>>>>>> gathered.w.pest
### Set density as a variable, to see if density dependence exists

# Mar 3,2017
### How many individuals for each species in each quadrat
gathered.w.pest3 <- gathered.w.pest

for (i in 1:dim(gathered.w.pest3)[1]){
  if (any(is.na(gathered.w.pest3$surv[i]))) {
    gathered.w.pest3$pp[i] <- NA
  } else {gathered.w.pest3$pp[i] <- 1}
}

######## Without considering growth.form
w.pest.surv.all  <- aggregate(pp ~ pesticide + exclosure + site + quadrat + census + species, 
                              data=gathered.w.pest3, FUN=sum)

### How many survivors for each species in each quadrat
w.pest.surv.dat  <- aggregate(surv ~ pesticide + exclosure + site + quadrat + census + species, 
                              data=gathered.w.pest, FUN=sum)

w.pest.surv.dat$all <- w.pest.surv.all$pp

w.pest.surv.dat$surv.rati <- w.pest.surv.dat$surv/ w.pest.surv.dat$all
hist(w.pest.surv.dat$surv.rati)

dev.new()
par(mfrow=c(1,2))
hist(w.pest.surv.dat$surv)
hist(sqrt(w.pest.surv.dat$surv))

### Fit the model

m.w.pest.surv.rati1 <- glmer(surv.rati ~ pesticide*exclosure + census +  (1|site/quadrat) + (1|species),
                             data=w.pest.surv.dat, family = binomial)
m.w.pest.surv.rati2 <- glmer(surv.rati ~ pesticide + exclosure + census +  (1|site/quadrat) + (1|species),
                             data=w.pest.surv.dat, family = binomial)
m.w.pest.surv.rati3 <- glmer(surv.rati ~ pesticide + exclosure + census + surv + (1|site/quadrat) + (1|species),
                             data=w.pest.surv.dat, weights=surv, family = binomial)
m.w.pest.surv.rati4 <- glmer(surv.rati ~ pesticide*surv + exclosure + census + (1|site/quadrat) + (1|species),
                             data=w.pest.surv.dat, weights=surv, family = binomial)
m.w.pest.surv.rati5 <- glmer(surv.rati ~ pesticide*surv + exclosure + census + (1|site/quadrat) + (1|species),
                             data=w.pest.surv.dat, weights=surv, family = binomial, 
                             control=glmerControl(optimizer="bobyqa"))


anova(m.w.pest.surv.rati1, m.w.pest.surv.rati2)

plot(m.w.pest.surv.rati2)
plot(m.w.pest.surv.rati3)
plot(m.w.pest.surv.rati4)

qqPlot(resid(m.w.pest.surv.rati4))
qqnorm(resid(m.w.pest.surv.rati4))
summary(m.w.pest.surv.rati4)

ranef(m.w.pest.surv.rati4)

drop1(m.w.pest.surv.rati4)


######>>>>>>>>>>>>>>>>>>>>> Considering growth.form <<<<<<<<<<<<<<<<<<######

w.pest.surv.all.g  <- aggregate(pp ~ pesticide + exclosure + site + quadrat + census + growth.form + species, 
                                data=gathered.w.pest3, FUN=sum)

### How many survivors for each species in each quadrat
w.pest.surv.dat.g  <- aggregate(surv ~ pesticide + exclosure + site + quadrat + census + growth.form + species, 
                                data=gathered.w.pest, FUN=sum)

w.pest.surv.dat.g$all <- w.pest.surv.all.g$pp

w.pest.surv.dat.g$surv.rati <- w.pest.surv.dat.g$surv/ w.pest.surv.dat.g$all
hist(w.pest.surv.dat.g$surv.rati)

### Fit the model

m.w.pest.g.surv.rati1 <- glmer(surv.rati ~ pesticide*exclosure + census + growth.form + (1|site/quadrat) + (1|species),
                               data=w.pest.surv.dat.g, family = binomial)
m.w.pest.g.surv.rati2 <- glmer(surv.rati ~ pesticide + exclosure + census + growth.form +   (1|site/quadrat) + (1|species),
                               data=w.pest.surv.dat.g, family = binomial)
m.w.pest.g.surv.rati3 <- glmer(surv.rati ~ pesticide + exclosure + census + growth.form + surv + (1|site/quadrat) + (1|species),
                               data=w.pest.surv.dat.g, weights=surv, family = binomial)
m.w.pest.g.surv.rati4 <- glmer(surv.rati ~ pesticide*surv + exclosure + census + growth.form + (1|site/quadrat) + (1|species),
                               data=w.pest.surv.dat.g, weights=surv, family = binomial)
m.w.pest.g.surv.rati5 <- glmer(surv.rati ~ pesticide*surv + exclosure + census + growth.form + (1|site/quadrat) + (1|species),
                               data=w.pest.surv.dat.g, weights=surv, family = binomial, 
                               control=glmerControl(optimizer="bobyqa"))

anova(m.w.pest.g.surv.rati1, m.w.pest.g.surv.rati2)

plot(m.w.pest.g.surv.rati1)
plot(m.w.pest.g.surv.rati3)
plot(m.w.pest.g.surv.rati4)

qqPlot(resid(m.w.pest.g.surv.rati4))
qqnorm(resid(m.w.pest.g.surv.rati4))

summary(m.w.pest.g.surv.rati4)



######>>>>>>>>>>>>>>>> gathered.w.snow
gathered.w.snow3 <- gathered.w.snow

for (i in 1:dim(gathered.w.snow3)[1]){
  if (any(is.na(gathered.w.snow3$surv[i]))) {
    gathered.w.snow3$pp[i] <- NA
  } else {gathered.w.snow3$pp[i] <- 1}
}

######## Considering growth.form
w.snow.surv.all.g  <- aggregate(pp ~ treatment +  site + quadrat + census + growth.form + species, 
                                data=gathered.w.snow3, FUN=sum)

### How many survivors for each species in each quadrat
w.snow.surv.dat.g  <- aggregate(surv ~ treatment + site + quadrat + census + growth.form + species, 
                                data=gathered.w.snow, FUN=sum)

w.snow.surv.dat.g$all <- w.snow.surv.all.g$pp

w.snow.surv.dat.g$surv.rati <- w.snow.surv.dat.g$surv/ w.snow.surv.dat.g$all
hist(w.snow.surv.dat.g$surv.rati)

### Initial density distribution
dev.new()
par(mfrow=c(1,2))
hist(w.snow.surv.dat.g$surv)
hist(sqrt(w.snow.surv.dat.g$surv))


### Fit the model

m.w.snow.g.surv.rati1 <- glmer(surv.rati ~ treatment + census + growth.form + (1|site/quadrat) + (1|species),
                               data=w.snow.surv.dat.g, family = binomial)
m.w.snow.g.surv.rati2 <- glmer(surv.rati ~ treatment*growth.form  + census + (1|site/quadrat) + (1|species),
                               data=w.snow.surv.dat.g, family = binomial)

anova(m.w.snow.g.surv.rati1, m.w.snow.g.surv.rati2)

### fit initial density to the model
m.w.snow.g.surv.rati3 <- glmer(surv.rati ~ treatment + census + growth.form + surv + (1|site/quadrat) + (1|species),
                               data=w.snow.surv.dat.g, weights=surv, family = binomial)
m.w.snow.g.surv.rati4 <- glmer(surv.rati ~ treatment*growth.form + census + surv + (1|site/quadrat) + (1|species),
                               data=w.snow.surv.dat.g, weights=surv, family = binomial)
m.w.snow.g.surv.rati5 <- glmer(surv.rati ~ treatment*surv + census + growth.form + (1|site/quadrat) + (1|species),
                               data=w.snow.surv.dat.g, weights=surv, family = binomial)

m.w.snow.g.surv.rati6 <- glmer(surv.rati ~ treatment*surv + census + growth.form + (1|site/quadrat) + (1|species),
                               data=w.snow.surv.dat.g, weights=surv, family = binomial,
                               control=glmerControl(optimizer="bobyqa"))

anova(m.w.snow.g.surv.rati3, m.w.snow.g.surv.rati4, m.w.snow.g.surv.rati5, m.w.snow.g.surv.rati6)
anova(m.w.snow.g.surv.rati1, m.w.snow.g.surv.rati4)

plot(m.w.snow.g.surv.rati1)
plot(m.w.snow.g.surv.rati4)

### Transform initial density with sqrt, fit this in the models
m.w.snow.g.surv.rati7 <- glmer(surv.rati ~ treatment + census + growth.form + sqrt(surv) + (1|site/quadrat) + (1|species),
                               data=w.snow.surv.dat.g, weights=surv, family = binomial)
m.w.snow.g.surv.rati8 <- glmer(surv.rati ~ treatment*growth.form + census + sqrt(surv) + (1|site/quadrat) + (1|species),
                               data=w.snow.surv.dat.g, weights=surv, family = binomial)
m.w.snow.g.surv.rati9 <- glmer(surv.rati ~ treatment*sqrt(surv) + census + growth.form + (1|site/quadrat) + (1|species),
                               data=w.snow.surv.dat.g, weights=surv, family = binomial)

m.w.snow.g.surv.rati10 <- glmer(surv.rati ~ treatment*sqrt(surv) + census + growth.form + (1|site/quadrat) + (1|species),
                               data=w.snow.surv.dat.g, weights=surv, family = binomial,
                               control=glmerControl(optimizer="bobyqa"))

anova(m.w.snow.g.surv.rati7, m.w.snow.g.surv.rati8, m.w.snow.g.surv.rati9, m.w.snow.g.surv.rati10)

plot(m.w.snow.g.surv.rati8)
qqnorm(resid(m.w.snow.g.surv.rati8))
summary(m.w.snow.g.surv.rati8)





#######################################################
####### Check the residuals with fitted value
#######################################################
## gathered.w.pest
test.w.pest <- data.frame(res = resid(m.w.pest.surv), fit = fitted(m.w.pest.surv), 
                      pest = gathered.w.pest0$pesticide, excl = gathered.w.pestt0$exclosure,
                      cens = gathered.w.pest0$census, site = gathered.w.pest0$site,
                      quad = gathered.w.pest0$quad1,
                      grfo = gathered.w.pest0$growth.form, spec = gathered.w.pest0$species)

## gathered.w.snow
test.w.snow <- data.frame(res = resid(m.w.snow.surv), fit = fitted(m.w.snow.surv), 
                         cens = gathered.w.snow0$census, site = gathered.w.snow0$site, 
                         grfo = gathered.w.snow0$growth.form, spec = gathered.w.snow0$species,
                         quad = gathered.w.snow0$quad1, tret = gathered.w.snow0$treatment)

## gathered.pest
#qplot(x=fit, y=res, geom='smooth', data=test.pest.all) + 
#  xlab("fitted")+
#  ylab("resi")+
#  ggtitle("the whole data (gathered.pest)")+ theme_bw() +
#  geom_hline(aes(yintercept=0)) + geom_point() +
#  geom_text(aes(label = quad), size = 3, 
#            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.w.pest) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("the whole data (gathered.w.pest)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.w.pest,facets=~pest) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("among pesticide trt (gathered.w.pest)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.w.pest,facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("in each census (gathered.w.pest)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.w.pest, colour=pest, facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("pesticide trt in each census (gathered.w.pest)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.w.pest, colour=pest, facets=~site) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("pesticide trt in each site (gathered.w.pest)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.w.pest, colour=pest, facets=~grfo) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("pesticide trt between growthform (gathered.w.pest)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.w.pest, colour=pest, facets=~excl) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("pesticidetrt between exclosure (gathered.w.pest)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)
##-------------------------------##


######>>>>>>>>>>>>>> gathered.snow
#qplot(x=fit, y=res, geom='smooth', data=test.w.snow) + 
#  xlab("fitted")+
#  ylab("resi")+
#  ggtitle("the whole data (gathered.w.snow)")+ theme_bw() +
#  geom_hline(aes(yintercept=0)) + geom_point() +
#  geom_text(aes(label = quad), size = 3, 
#            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.w.snow) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("the whole data (gathered.w.snow)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() 

qplot(x=fit, y=res, geom='smooth', data=test.w.snow,facets=~tret) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("among snow trt (gathered.w.snow)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.w.snow,facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("in each census (gathered.w.snow)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.w.snow, colour=tret, facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("snow trt in each census (gathered.w.snow)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.w.snow, colour=tret, facets=~site) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("snow trt in each site (gathered.w.snow)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.w.snow, colour=tret, facets=~grfo) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("snow trt between growthform (gathered.w.snow)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)
##-------------------------------##


## the AUC (Area Under the Curve of the Receiver Operator Characteristic
roc.area(obs=getME(m.g.all.surv, 'y'), pred=fitted(m.g.all.surv))
roc.area(obs=getME(m.snow.all.surv, 'y'), pred=fitted(m.snow.all.surv))
roc.area(obs=getME(m.pest.new.all.surv, 'y'), pred=fitted(m.pest.new.all.surv))





#####################################################
## gathered.w.pest
## analyze growth.form (trees and shrubs) seperately
#####################################################

## need to refresh the codes

gathered.dat.tre <- gathered.dat[gathered.dat$growth.form == "tree",]
gathered.dat.shr <- gathered.dat[gathered.dat$growth.form == "shrub",]

summary(gathered.dat.tre)
summary(gathered.dat.tre)

## dat.all.tre
gathered.dat.tre$site <- factor(gathered.dat.tre$site)
gathered.dat.tre$code <- factor(gathered.dat.tre$code)
gathered.dat.tre$quadrat <- factor(gathered.dat.tre$quadrat)
gathered.dat.tre$tag <- factor(gathered.dat.tre$tag)
gathered.dat.tre$exclosure <- factor(gathered.dat.tre$exclosure)

## dat.all.shr
gathered.dat.shr$site <- factor(gathered.dat.shr$site)
gathered.dat.shr$code <- factor(gathered.dat.shr$code)
gathered.dat.shr$quadrat <- factor(gathered.dat.shr$quadrat)
gathered.dat.shr$tag <- factor(gathered.dat.shr$tag)
gathered.dat.shr$exclosure <- factor(gathered.dat.shr$exclosure)

## combine the code and quadrat number as real quadrat
gathered.dat.tre$quad1 <-  as.factor(gathered.dat.tre$code):gathered.dat.tre$quadrat
gathered.dat.shr$quad1 <-  as.factor(gathered.dat.shr$code):gathered.dat.shr$quadrat

## non NAs data set
gathered.dat.tre0 <- subset(gathered.dat.tre,!is.na(gathered.dat.tre$surv))
gathered.dat.shr0 <- subset(gathered.dat.shr,!is.na(gathered.dat.shr$surv))

#######################
## fitting in the mixed-effect model
# gathered.dat.tre
m.pest.tre.surv <- glmer(surv ~ pesticide*exclosure + census + (1|site/quadrat) + (1|species),
                      data=gathered.dat.tre, family = binomial)
summary(m.pest.tre.surv)
# gathered.dat.shr
m.pest.shr.surv <- glmer(surv ~ pesticide*exclosure + census + (1|site/quadrat) + (1|species),
                         data=gathered.dat.shr, family = binomial)
summary(m.pest.shr.surv)

## the AUC (Area Under the Curve of the Receiver Operator Characteristic
roc.area(obs=getME(m.pest.tre.surv, 'y'), pred=fitted(m.pest.tre.surv))
roc.area(obs=getME(m.pest.shr.surv, 'y'), pred=fitted(m.pest.shr.surv))


## gathered.dat.tre
test.pest.tre <- data.frame(res = resid(m.pest.tre.surv), fit = fitted(m.pest.tre.surv), 
                         pest = gathered.dat.tre0$pesticide,
                         excl = gathered.dat.tre0$exclosure, cens = gathered.dat.tre0$census,
                         site = gathered.dat.tre0$site, quad = gathered.dat.tre0$quad1,
                         grfo = gathered.dat.tre0$growth.form, spec = gathered.dat.tre0$species) 

## gathered.dat.shr
test.pest.shr <- data.frame(res = resid(m.pest.shr.surv), fit = fitted(m.pest.shr.surv), 
                         pest = gathered.dat.shr0$pesticide,
                         excl = gathered.dat.shr0$exclosure, cens = gathered.dat.shr0$census,
                         site = gathered.dat.shr0$site, quad = gathered.dat.shr0$quad1,
                         grfo = gathered.dat.shr0$growth.form, spec = gathered.dat.shr0$species) 
######
## gathered.dat.tre
qplot(x=fit, y=res, geom='smooth', data=test.pest.tre) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("the whole data (gathered.dat.tre)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.pest.tre,facets=~pest) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("among pest trt (gathered.dat.tre)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.pest.tre,facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("in each census (gathered.dat.tre)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.pest.tre, colour=pest, facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("pesticide trt in each census (gathered.dat.tre)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.pest.tre, colour=pest, facets=~site) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("pesticide trt in each site (gathered.dat.tre)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.pest.tre, colour=pest, facets=~excl) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("pesticidetrt between exclosure (gathered.dat.tre)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)
##-------------------------------##

######
## gathered.dat.shr
qplot(x=fit, y=res, geom='smooth', data=test.pest.shr) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("the whole data (gathered.dat.shr)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.pest.shr,facets=~pest) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("between pesticide trt (gathered.dat.shr)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.pest.shr,facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("in each census (gathered.dat.shr)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.pest.shr, colour=pest, facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("pesticide trt in each census (gathered.dat.shr)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.pest.shr, colour=pest, facets=~site) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("pesticide trt in each site (gathered.dat.shr)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.pest.shr, colour=pest, facets=~excl) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("pesticide trt between exclosure (gathered.dat.shr)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

##-------------------------------------##


##################### Overdispersion? ######################
## gathered.pest
plot(m.w.pest.g.surv.rati4)
chisq.pest <- sum(resid(m.w.pest.g.surv.rati4, type='pearson')^2)
chisq.pest/df.residual(m.w.pest.g.surv.rati4)  ## 0.92, Less than 1
## significantly so?
1-pchisq(chisq.pest, df.residual(m.w.pest.g.surv.rati4)) # 1,No

## gathered.snow
#plot(m.snow.all.surv)
#chisq.snow <- sum(resid(m.snow.all.surv, type='pearson')^2)
#chisq.snow/df.residual(m.snow.all.surv)  ## 0.97, Less than 1
## significantly so?
#1-pchisq(chisq.snow, df.residual(m.snow.all.surv)) # 0.8, No




###########################
## get the BLUPs associated with the various predictors
# observed mean value of for each quadrat
#pest.surv.mean <- aggregate(surv ~ site + quadrat + census,
#                       data=gathered.pest, FUN=mean)
pest.surv.mean <- aggregate(surv ~ site + quadrat + pesticide,
                            data=gathered.pest, FUN=mean)

re.pest.surv <- ranef(m.pest.all.surv)
str(re.pest.surv)





########## Random effects distribution ###########
re.pest.surv <- ranef(m.pest.all.surv,postVar = TRUE)
str(re.pest.surv)
dotplot(re.pest.surv,postVar = TRUE)[["quadrat:site"]]
dotplot(re.pest.surv,postVar = TRUE)[["species"]]
dotplot(re.pest.surv,postVar = TRUE)[["site"]]

##----------------------------------------------##






#pest.surv.quad <- cbind(re.pest.surv.quad, pest.surv.mean)
# unite site and quadrat number together
#pest.surv.quad <- unite(pest.surv.quad, quad, c(site, quadrat), sep = "")

#pest.surv.quad$re.pest.quad <- pest.surv.quad$`(Intercept)`
#pest.surv.quad$pest.surv.quad <- pest.surv.quad$surv

#ggplot(pest.surv.quad, aes(x=re.pest.quad, y=pest.surv.quad,label = quad)) + 
#  geom_point() + theme_bw() +
#  geom_smooth(method = "lm") + 
#  geom_text(aes(label = quad),hjust=0, vjust=0)

### model2 with growth form

#re.surv16a2 <- ranef(m.surv16a2)
#str(re.surv16a2)
#re.surv16a2.quad <- re.surv16a2$quad

#surv16a2.quad <- cbind(re.surv16a2.quad, surv16a.mean)
# unite site and quadrat number together
#surv16a2.quad <- unite(surv16a2.quad, quad, c(site, quadrat), sep = "")

#surv16a2.quad$re.quad <- surv16a2.quad$`(Intercept)`
#surv16a2.quad$surv16a2.quad <- surv16a2.quad$surv16a

#ggplot(surv16a2.quad, aes(x=re.quad, y=surv16a2.quad,label = quad)) + 
#  geom_point() + theme_bw() +
#  geom_smooth(method = "lm") + 
#  geom_text(aes(label = quad), size = 3, angle = 45, hjust=0.1, vjust=0.2)


## try to solve the overlapping text problem
## use labels
#library("ggrepel")
#p +   geom_text(aes(label = quad), size = 3, 
#                angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

#p + geom_point() + theme_bw() +
#  ggrepel::geom_text_repel(aes(label = quad), color = "black", size = 2.5, segment.color = "grey")
# it does not work


##--------------------------##
### Normality of residuals
hist(residuals(m.pest.all.surv))
qqnorm(residuals(m.pest.all.surv))

hist(residuals(m.snow.all.surv))
qqnorm(residuals(m.snow.all.surv))
##-----------------------##





######
## Extract the predicted value from the glmer model
## only including the fixed-effect

## all the results in all census
new.w.pest <-  expand.grid(pesticide =c("C", "F", "I", "W"), growth.form= c("tree","shrub"),
                         exclosure=c("0", "1"), census= c("16sp","16fa"))
# fit.pest <- predict(m.pest.all.surv, newdata = new.pest, type = "link",re.form=~0)
# fit.pest1 <- predict(m.pest.all.surv, newdata = new.pest, type = "response",re.form=~0)
# fit.pest2 <- predict(m.pest.all.surv, newdata = new.pest, type = "response",re.form=~(1|site/quadrat) + (1|species))
# pred.pest <- cbind (new.pest, fit.pest)

new.w.pest$quadrat <- gathered.w.pest$quadrat[1]
new.w.pest$site <- gathered.w.pest$site[1]
new.w.pest$species <- gathered.w.pest$species[1]

##################################
## get the confidence intervals

################## merTools::predictInterval #####################
library(httpuv)
library(merTools)
library(foreach)

predInt.w.pest <- predictInterval(m.w.pest.surv, newdata = new.pest, n.sims = 999, level = 0.95,
                                type = "linear.prediction", stat = "mean", which="fixed")

predprob.w.pest <- within(predInt.w.pest, {
  Predprob <- plogis(fit)
  LL <- plogis(lwr)
  UL <- plogis(upr)
})

pred.w.pest.dat <- cbind(new.w.pest, predprob.w.pest)


ggplot(pred.w.pest.dat, aes(x = pesticide, y = Predprob, ymin = LL, ymax = UL,
                          group=interaction(census, exclosure), col=census, shape=exclosure)) + 
  geom_pointrange() + geom_line(aes(y=Predprob, lty=exclosure), position=pd, size=0.8) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("growth.form") + theme_bw()

ggplot(pred.w.pest.dat, aes(x = pesticide, y = Predprob, ymin = LL, ymax = UL,
                          group=interaction(census, exclosure), col=census, shape=exclosure)) + 
  geom_pointrange(position=pd) + 
  #geom_line(aes(y=Predprob, lty=exclosure), position=pd, size=0.8) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("growth.form") + theme_bw()


############### Alternatively, the probability prediction ###################
predInt.w.pest.prob <- predictInterval(m.w.pest.surv, newdata = new.w.pest, n.sims = 999, level = 0.95,
                                     type = "probability", stat = "mean", which="fixed")
predInt.w.pest.prob.dat <- cbind(new.w.pest, predInt.w.pest.prob)

ggplot(predInt.w.pest.prob.dat, aes(x = pesticide, y = fit, ymin = lwr, ymax = upr,
                                  group=interaction(census, exclosure), col=census, shape=exclosure)) + 
  geom_pointrange(position=pd) + 
  #geom_line(aes(y=Predprob, lty=exclosure), position=pd, size=0.8) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("growth.form") + theme_bw()

### the confidence bars are too wide
### If the approach of predictInterval cause this
### Use another bootstrap method

################ bootMer ########################

predFunc <- function(.){ predict(., newdata=new.w.pest, re.form=~0) }

bootfit.w.pest <- bootMer(m.w.pest.surv, FUN=function(x)predict(x, new.w.pest, re.form=~0), nsim=99)

new.w.pest$lci <- apply(bootfit.w.pest$t, 2, quantile, 0.025) 
new.w.pest$uci <- apply(bootfit.w.pest$t, 2, quantile, 0.975) 
new.w.pest$survpred <- predict(m.w.pest.surv, newdata= new.w.pest, re.form=~0)

pred.w.pest <- within(new.w.pest, {
  Predprob <- plogis(survpred)
  LL <- plogis(lci)
  UL <- plogis(uci)
})


## Plot the results, using ggplot2
# The errorbars overlapped, so use position_dodge to move them horizontally

ggplot(pred.w.pest, aes(x = pesticide, y = Predprob, ymin = LL, ymax = UL,
                      group=interaction(census, exclosure), col=census, shape=exclosure)) + 
  geom_pointrange(position=pd) + 
  #geom_line(aes(y=Predprob, lty=exclosure), position=pd, size=0.8) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("growth.form") + theme_bw() +
  ylab("Survival")

ggplot(pred.w.pest, aes(x = pesticide, y = Predprob, ymin = LL, ymax = UL,
                      group=interaction(growth.form, exclosure), col=growth.form, shape=exclosure)) + 
  geom_pointrange(position=pd) + 
  #geom_line(aes(y=Predprob, lty=exclosure), position=pd, size=0.8) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("census") + theme_bw() +
  ylab("Survival")





#################
### SNOW DATA ###
#################
new.w.snow <- expand.grid(treatment = c("S", "C"), growth.form= c("tree","shrub"),
                        census= c("15sp","15fa","16sp","16fa"))

new.w.snow$quadrat <- gathered.w.snow$quadrat[1]
new.w.snow$site <- gathered.w.snow$site[1]
new.w.snow$species <- gathered.w.snow$species[1]

################## merTools::predictInterval #####################

predInt.w.snow <- predictInterval(m.w.snow.surv, newdata = new.w.snow , n.sims = 999, level = 0.95,
                                type = "linear.prediction", stat = "mean", which="fixed")

predIntprob.w.snow <- within(predInt.w.snow, {
  Predprob <- plogis(fit)
  LL <- plogis(lwr)
  UL <- plogis(upr)
})

predInt.w.snow.dat <- cbind(new.w.snow, predIntprob.w.snow)

ggplot(predInt.w.snow.dat, aes(x = treatment, y = Predprob, ymin = LL, ymax = UL,
                          color=as.factor(census))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("growth.form") + theme_bw()


#################### bootMer #########################

bootfit.w.snow <- bootMer(m.w.snow.surv, FUN=function(x)predict(x, new.w.snow, re.form=~0), nsim=99)

new.w.snow$lci <- apply(bootfit.w.snow$t, 2, quantile, 0.025) 
new.w.snow$uci <- apply(bootfit.w.snow$t, 2, quantile, 0.975) 
new.w.snow$survpred <- predict(m.w.snow.surv, newdata= new.w.snow, re.form=~0)

pred.w.snow <- within(new.w.snow, {
  Predprob <- plogis(survpred)
  LL <- plogis(lci)
  UL <- plogis(uci)
})


## Plot the results, using ggplot2
# The errorbars overlapped, so use position_dodge to move them horizontally

ggplot(pred.w.snow, aes(x = treatment, y = Predprob, ymin = LL, ymax = UL,
                      color=as.factor(census))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("growth.form") + theme_bw() +
  ylab("Survival")

ggplot(pred.w.snow, aes(x = treatment, y = Predprob, ymin = LL, ymax = UL,
                      color=as.factor(growth.form))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("census") + theme_bw() +
  ylab("Survival")

ggplot(pred.w.snow, aes(x = census, y = Predprob, ymin = LL, ymax = UL,
                      color=as.factor(treatment))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("growth.form") + theme_bw() +
  ylab("Survival")

##--------------------------------------------##




#################################################
## As we consider the census in the fixed effect, we get two different results in each census
## However, we want to know the overall trend acorss the census
## We need to figure out how to get the overall value
## Let's try the contrasts function, which I don't really understand

contrasts(gathered.w.pest$census) <- "contr.sum"
contrasts(gathered.w.pest$census) <- contr.sum(3)
contrasts(gathered.w.pest$census) <- contr.sum(levels(gathered.w.pest$census))
options(contrasts = c("contr.sum", "contr.poly"))
contrasts(gathered.w.pest$census)

# gathered.pest$census <- relevel(gathered.pest$census, ref="15f")

m.w.pest.surv.cont <- glmer(surv ~ pesticide*exclosure + census + growth.form +  (1|site/quadrat) + (1|species),
                         data=gathered.w.pest, family = binomial)
summary(m.w.pest.all.surv.cont)
## Warning: contrasts dropped from factor census due to missing levels
## I guess it's the census 15f, because the survival in this census is considered to be NAs
## HERE I CNA NOT USE THE CONTRASTS TO COMBINE THE RESULTS OF TWO CENSUSES




##################################################################
###### Survival analysis                               
###### Proportion of survival
###### All the treatements
###### Weighted by inital density
##################################################################

##########>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>########################
### Investigate survival seedlings in census 16fa when them existed in 16sp
### proportion of survival

w <- read.csv("Changbaishan/data/Woodyseedlings 20170301.csv")

# 16sp dataset for existing seedlings
w.16sp <- subset(w, w$heig_16sp !=0 & !is.na(w$heig_16sp))

w.16sp$surv.16sp <- as.numeric(w.16sp$heig_16sp !=0 & !is.na(w.16sp$heig_16sp))
w.16sp$surv.16fa <- as.numeric(w.16sp$heig_16fa !=0 & !is.na(w.16sp$heig_16fa))

## get mean survival rate for each quadrat for each treatment in different exclosure treatments
w.16sp.surv  <- aggregate(cbind(surv.16sp,surv.16fa)~ treatment + exclosure + site + quadrat, data=w.16sp, FUN=sum)

#View(dat.15b.surv)
w.16sp.surv$surv.rati <- w.16sp.surv$surv.16fa/w.16sp.surv$surv.16sp
summary(w.16sp.surv$surv.rati)

hist(w.16sp.surv$surv.rati)
hist(asin(w.16sp.surv$surv.rati))
# hist(log10(1-w.16sp.surv$surv.rati))

dev.new()
par(mfrow=c(1,2))
hist(w.16sp.surv$surv.rati,main = "Histogram of survival")
hist(asin(w.16sp.surv$surv.rati),main = "Histogram of survival(arcsin)")


### Since it's difficult to transform the proportion of survival to be normal distribution
### Next step, try get values for each species
## get mean survival rate for each quadrat for each treatment in different exclosure treatments
w.16sp.surv1  <- aggregate(cbind(surv.16sp,surv.16fa)~ treatment + exclosure + site + quadrat + species, data=w.16sp, FUN=sum)

#View(dat.15b.surv)
w.16sp.surv1$surv.rati <- w.16sp.surv1$surv.16fa/ (w.16sp.surv1$surv.16sp)

hist(w.16sp.surv1$surv.rati)
hist(asin(w.16sp.surv1$surv.rati))
# hist(log10(1-w.16sp.surv1$surv.rati))

### it's worse

par(mfrow=c(1,2))
hist(w.16sp.surv1$surv.rati,main = "Histogram of survival for species")
hist(log10(1-w.16sp.surv1$surv.rati),main = "Histogram of survival(log10(1-x)) for species")


#####>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
########## Only tree species?

w.16sp.tree <- w.16sp[w.16sp$growth.form=="tree",]
## how many species
length(unique(w.16sp.tree$species))
w.16sp.tree$surv.16sp <- as.numeric(w.16sp.tree$heig_16sp !=0 & !is.na(w.16sp.tree$heig_16sp))
w.16sp.tree$surv.16fa <- as.numeric(w.16sp.tree$heig_16fa !=0 & !is.na(w.16sp.tree$heig_16fa))

### First, for each species
## get mean survival rate for each quadrat for each treatment in different exclosure treatments
w.16sp.tree.surv  <- aggregate(cbind(surv.16sp,surv.16fa)~ treatment + exclosure + site + quadrat + species, data=w.16sp.tree, FUN=sum)

w.16sp.tree.surv$surv.rati <- w.16sp.tree.surv$surv.16fa/w.16sp.tree.surv$surv.16sp
summary(w.16sp.tree.surv$surv.rati)

dev.new()
par(mfrow=c(2,2))
hist(w.16sp.tree.surv$surv.rati, main="")
hist(asin(w.16sp.tree.surv$surv.rati), main="")
hist((w.16sp.tree.surv$surv.rati)^3, main="")
hist(exp(w.16sp.tree.surv$surv.rati), main="")

## Hard to transform to normal distribution

### Second, for whole species
## get mean survival rate for each quadrat for each treatment in different exclosure treatments
w.16sp.tree.surv1  <- aggregate(cbind(surv.16sp,surv.16fa)~ treatment + exclosure + site + quadrat, data=w.16sp.tree, FUN=sum)

w.16sp.tree.surv1$surv.rati <- w.16sp.tree.surv1$surv.16fa/(w.16sp.tree.surv1$surv.16sp)
summary(w.16sp.tree.surv1$surv.rati)

par(mfrow=c(2,2))
hist(w.16sp.tree.surv1$surv.rati, main="")
hist(asin(w.16sp.tree.surv1$surv.rati), main="")
hist((w.16sp.tree.surv1$surv.rati)^3, main="")
hist(exp(w.16sp.tree.surv1$surv.rati), main="")



################ Two dorminant species
### Fraxinus mandschurica
### Tilia amurensis
frax <- w[w$species == "Fraxinus mandschurica",]
tili <- w[w$species == "Tilia amurensis",]

# 16sp dataset for existing seedlings
frax.16sp <- subset(frax, frax$heig_16sp !=0 & !is.na(frax$heig_16sp))

frax.16sp$surv.16sp <- as.numeric(frax.16sp$heig_16sp !=0 & !is.na(frax.16sp$heig_16sp))
frax.16sp$surv.16fa <- as.numeric(frax.16sp$heig_16fa !=0 & !is.na(frax.16sp$heig_16fa))

## get mean survival rate for each quadrat for each treatment in different exclosure treatments
frax.16sp.surv  <- aggregate(cbind(surv.16sp,surv.16fa)~ treatment + exclosure + site + quadrat, data=frax.16sp, FUN=sum)

#View(dat.15b.surv)
frax.16sp.surv$surv.rati <- frax.16sp.surv$surv.16fa/frax.16sp.surv$surv.16sp
summary(frax.16sp.surv$surv.rati)

hist(frax.16sp.surv$surv.rati)

# 16sp dataset for existing seedlings
tili.16sp <- subset(tili, tili$heig_16sp !=0 & !is.na(tili$heig_16sp))

tili.16sp$surv.16sp <- as.numeric(tili.16sp$heig_16sp !=0 & !is.na(tili.16sp$heig_16sp))
tili.16sp$surv.16fa <- as.numeric(tili.16sp$heig_16fa !=0 & !is.na(tili.16sp$heig_16fa))

## get mean survival rate for each quadrat for each treatment in different exclosure treatments
tili.16sp.surv  <- aggregate(cbind(surv.16sp,surv.16fa)~ treatment + exclosure + site + quadrat, data=tili.16sp, FUN=sum)

#View(dat.15b.surv)
tili.16sp.surv$surv.rati <- tili.16sp.surv$surv.16fa/tili.16sp.surv$surv.16sp
summary(tili.16sp.surv$surv.rati)

hist(tili.16sp.surv$surv.rati)
hist(asin(tili.16sp.surv$surv.rati))
### both of these two species experienced high survival from spring to fall

par(mfrow = c(1,2))
hist(frax.16sp.surv$surv.rati, main="Histogram of Fraxius")
hist(tili.16sp.surv$surv.rati, main="Histogram of Tilia")



####>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
### fit in the models
### First, the whole data

m.w.16sp.surv1 <- lmer(asin(surv.rati)~ treatment*exclosure + (1|site/quadrat),
                       data=w.16sp.surv)

plot(m.w.16sp.surv1)

m.w.16sp.surv2 <- lmer(asin(surv.rati)~ treatment + exclosure + (1|site/quadrat),
                       data=w.16sp.surv)
plot(m.w.16sp.surv2)

anova(m.w.16sp.surv1, m.w.16sp.surv2)

drop1(m.w.16sp.surv2)


library(car)
qqPlot(m.w.16sp.surv2)
qqnorm(m.w.16sp.surv2)

summary(m.w.16sp.surv2)

w.16sp.surv$treatment <- relevel(w.16sp.surv$treatment, ref="W")
m.w.16sp.surv <- lmer(asin(surv.rati)~ treatment + exclosure + (1|site/quadrat),
                      data=w.16sp.surv)
summary(m.w.16sp.surv)

m.w.16sp.surv3 <- lmer(asin(surv.rati)~ treatment + exclosure + (1|site/quadrat),
                       data=w.16sp.surv, control=glmerControl(optimizer="bobyqa"))
summary(m.w.16sp.surv3)

### weighted with initial density of the whole quadrat
m.w.16sp.surv4 <- lmer(asin(surv.rati)~ treatment + exclosure + (1|site/quadrat),
                       data=w.16sp.surv, weights=surv.16sp)
summary(m.w.16sp.surv4)

m.w.16sp.surv5 <- lmer(asin(surv.rati)~ treatment + (1|site/quadrat),
                       data=w.16sp.surv)

anova(m.w.16sp.surv, m.w.16sp.surv5)

m.w.16sp.surv6 <- lmer(asin(surv.rati)~ treatment + (1|site/quadrat),
                       data=w.16sp.surv, weights=surv.16sp)
summary(m.w.16sp.surv6)


plot(m.w.16sp.surv)
plot(m.w.16sp.surv4)
plot(m.w.16sp.surv6)

dev.new()
par(mfrow=c(1,2))
library(car)
qqPlot(resid(m.w.16sp.surv), main="Quantile-Comparison plot for survival")
qqPlot(resid(m.w.16sp.surv4),
       main="Quantile-Comparison plot for survival weighted")
qqPlot(resid(m.w.16sp.surv6),
       main="Quantile-Comparison plot for survival weighted(trt)")


############### Not transform the survival
m.w.16sp.surv7 <- lmer(surv.rati~ treatment*exclosure + (1|site/quadrat),
                       data=w.16sp.surv)
m.w.16sp.surv8 <- lmer(surv.rati~ treatment + exclosure + (1|site/quadrat),
                       data=w.16sp.surv)
m.w.16sp.surv9 <- lmer(surv.rati~ treatment*exclosure + surv.16sp + (1|site/quadrat),
                       data=w.16sp.surv)
m.w.16sp.surv10 <- lmer(surv.rati~ treatment + exclosure + surv.16sp + (1|site/quadrat),
                        data=w.16sp.surv)
m.w.16sp.surv11 <- lmer(surv.rati~ treatment*surv.16sp + exclosure + (1|site/quadrat),
                        data=w.16sp.surv)
m.w.16sp.surv12 <- lmer(surv.rati~ treatment + exclosure + surv.16sp + (1|site/quadrat),
                        data=w.16sp.surv, weights=surv.16sp)
m.w.16sp.surv13 <- lmer(surv.rati~ treatment*surv.16sp + exclosure + (1|site/quadrat),
                        data=w.16sp.surv, weights=surv.16sp)

anova(m.w.16sp.surv7, m.w.16sp.surv8)
anova(m.w.16sp.surv9, m.w.16sp.surv10)
anova(m.w.16sp.surv7, m.w.16sp.surv9)
anova(m.w.16sp.surv7, m.w.16sp.surv11)

plot(m.w.16sp.surv7)
plot(m.w.16sp.surv8)
plot(m.w.16sp.surv9)
plot(m.w.16sp.surv10)
plot(m.w.16sp.surv11)
plot(m.w.16sp.surv12)
plot(m.w.16sp.surv13)

library(car)
qqPlot(resid(m.w.16sp.surv7))
qqPlot(resid(m.w.16sp.surv8))
qqPlot(resid(m.w.16sp.surv9))
qqPlot(resid(m.w.16sp.surv10))
qqPlot(resid(m.w.16sp.surv11))
qqPlot(resid(m.w.16sp.surv12))
qqPlot(resid(m.w.16sp.surv13))



####>>>>>>>>>>>>>>>>>>>>>> For each species

############### Not transform the survival
m.w.16sp.surv.spec1 <- lmer(surv.rati~ treatment*exclosure + (1|site/quadrat) + (1|species),
                            data=w.16sp.surv1)
m.w.16sp.surv.spec2 <- lmer(surv.rati~ treatment + exclosure + (1|site/quadrat) + (1|species),
                            data=w.16sp.surv1)
m.w.16sp.surv.spec3 <- lmer(surv.rati~ treatment*exclosure + surv.16sp + (1|site/quadrat) + (1|species),
                            data=w.16sp.surv1)
m.w.16sp.surv.spec4 <- lmer(surv.rati~ treatment + exclosure + surv.16sp + (1|site/quadrat) + (1|species),
                            data=w.16sp.surv1)
m.w.16sp.surv.spec5 <- lmer(surv.rati~ treatment*surv.16sp  + exclosure + (1|site/quadrat) + (1|species),
                            data=w.16sp.surv1)
m.w.16sp.surv.spec6 <- lmer(surv.rati~ treatment*exclosure + surv.16sp + (1|site/quadrat) + (1|species),
                            data=w.16sp.surv1, weights=surv.16sp)
m.w.16sp.surv.spec7 <- lmer(surv.rati~ treatment + exclosure + surv.16sp + (1|site/quadrat) + (1|species),
                            data=w.16sp.surv1, weights=surv.16sp)
m.w.16sp.surv.spec8 <- lmer(surv.rati~ treatment*surv.16sp  + exclosure + (1|site/quadrat) + (1|species),
                            data=w.16sp.surv1, weights=surv.16sp)

anova(m.w.16sp.surv.spec1, m.w.16sp.surv.spec2)
anova(m.w.16sp.surv.spec1, m.w.16sp.surv.spec3)
anova(m.w.16sp.surv.spec1, m.w.16sp.surv.spec4)
anova(m.w.16sp.surv.spec1, m.w.16sp.surv.spec5)
anova(m.w.16sp.surv.spec6, m.w.16sp.surv.spec7)
anova(m.w.16sp.surv.spec7, m.w.16sp.surv.spec8)


plot(m.w.16sp.surv.spec1)
plot(m.w.16sp.surv.spec2)
plot(m.w.16sp.surv.spec3)
plot(m.w.16sp.surv.spec4)
plot(m.w.16sp.surv.spec5)
plot(m.w.16sp.surv.spec6)
plot(m.w.16sp.surv.spec7)

library(car)
qqPlot(resid(m.w.16sp.surv.spec1))
qqPlot(resid(m.w.16sp.surv.spec3))
qqPlot(resid(m.w.16sp.surv.spec5))
qqPlot(resid(m.w.16sp.surv.spec6))
qqPlot(resid(m.w.16sp.surv.spec7))
qqPlot(resid(m.w.16sp.surv.spec8))

######## glmer with binomial distribution
w.16sp.surv1$treatment <- relevel(w.16sp.surv1$treatment, ref="W") 
m.w.16sp.surv.spec11 <- glmer(surv.rati~ treatment*exclosure + (1|site/quadrat) + (1|species),
                              data=w.16sp.surv1, family="binomial")
m.w.16sp.surv.spec12 <- glmer(surv.rati~ treatment + exclosure + (1|site/quadrat) + (1|species),
                              data=w.16sp.surv1, family="binomial")
m.w.16sp.surv.spec13 <- glmer(surv.rati~ treatment*exclosure + surv.16sp + (1|site/quadrat) + (1|species),
                              data=w.16sp.surv1, family="binomial")
m.w.16sp.surv.spec14 <- lmer(surv.rati~ treatment + exclosure + surv.16sp + (1|site/quadrat) + (1|species),
                             data=w.16sp.surv1)
m.w.16sp.surv.spec15 <- lmer(surv.rati~ treatment*surv.16sp  + exclosure + (1|site/quadrat) + (1|species),
                             data=w.16sp.surv1)
m.w.16sp.surv.spec16 <- lmer(surv.rati~ treatment*exclosure + surv.16sp + (1|site/quadrat) + (1|species),
                             data=w.16sp.surv1, weights=surv.16sp)
m.w.16sp.surv.spec17 <- lmer(surv.rati~ treatment + exclosure + surv.16sp + (1|site/quadrat) + (1|species),
                             data=w.16sp.surv1, weights=surv.16sp)
m.w.16sp.surv.spec18 <- glmer(surv.rati~ treatment*surv.16sp  + exclosure + (1|site/quadrat) + (1|species),
                              data=w.16sp.surv1, weights=surv.16sp, family="binomial")

summary(m.w.16sp.surv.spec11)
summary(m.w.16sp.surv.spec18)

plot(m.w.16sp.surv.spec11)
plot(m.w.16sp.surv.spec18)

qqPlot(resid(m.w.16sp.surv.spec11))
qqPlot(resid(m.w.16sp.surv.spec18))



##---------------------------###
### Second, only tree species
### Whole in quadrat
m.w.16sp.tree.surv1 <- lmer(asin(surv.rati)~ treatment*exclosure + (1|site/quadrat),
                            data=w.16sp.tree.surv)
m.w.16sp.tree.surv2 <- lmer(asin(surv.rati)~ treatment + exclosure + (1|site/quadrat),
                            data=w.16sp.tree.surv)
m.w.16sp.tree.surv3 <- lmer(asin(surv.rati)~ treatment + (1|site/quadrat),
                            data=w.16sp.tree.surv)
anova(m.w.16sp.tree.surv1, m.w.16sp.tree.surv2)
anova(m.w.16sp.tree.surv1, m.w.16sp.tree.surv3)

plot(m.w.16sp.tree.surv1)

m.w.16sp.tree.surv4 <- lmer(asin(surv.rati)~ treatment*exclosure + (1|site/quadrat),
                            data=w.16sp.tree.surv, weights=surv.16sp)

plot(m.w.16sp.tree.surv4)

w.16sp.tree.surv$treatment <- relevel(w.16sp.tree.surv$treatment, ref="W")
m.w.16sp.tree.surv <- lmer(asin(surv.rati)~ treatment*exclosure + (1|site/quadrat),
                           data=w.16sp.tree.surv, weights=surv.16sp)
summary(m.w.16sp.tree.surv)


###### Without transforming surv.rati
w.16sp.tree.surv$treatment <- relevel(w.16sp.tree.surv$treatment, ref="W") 
m.w.16sp.tree.surv5 <- glmer(surv.rati~ treatment*exclosure + (1|site/quadrat),
                            data=w.16sp.tree.surv, family="binomial")
m.w.16sp.tree.surv6 <- glmer(surv.rati~ treatment + exclosure + (1|site/quadrat),
                             data=w.16sp.tree.surv, family="binomial")
m.w.16sp.tree.surv7 <- glmer(surv.rati~ treatment + (1|site/quadrat),
                             data=w.16sp.tree.surv, family="binomial")
m.w.16sp.tree.surv8 <- glmer(surv.rati~ treatment + exclosure + (1|site/quadrat),
                             data=w.16sp.tree.surv, weights=surv.16sp, family="binomial")
m.w.16sp.tree.surv9 <- glmer(surv.rati~ treatment*exclosure + (1|site/quadrat),
                             data=w.16sp.tree.surv, weights=surv.16sp, family="binomial")


anova(m.w.16sp.tree.surv5, m.w.16sp.tree.surv6)
anova(m.w.16sp.tree.surv5, m.w.16sp.tree.surv7)
anova(m.w.16sp.tree.surv5, m.w.16sp.tree.surv8)
anova(m.w.16sp.tree.surv5, m.w.16sp.tree.surv9)

plot(m.w.16sp.tree.surv5)
qqnorm(resid(m.w.16sp.tree.surv5))



########################### bootstrap ################################
### Whole species
### both treatment and exclosure
new.w.16sp.surv <- expand.grid(treatment = c("C", "F","I","W","S"), exclosure = c("0","1"))

bootfit.w.16sp.surv <- bootMer(m.w.16sp.surv4, 
                               FUN=function(x)predict(x, new.w.16sp.surv, re.form=~0), nsim=999)

new.w.16sp.surv$lci <- apply(bootfit.w.16sp.surv$t, 2, quantile, 0.025) 
new.w.16sp.surv$uci <- apply(bootfit.w.16sp.surv$t, 2, quantile, 0.975) 
new.w.16sp.surv$survpred <- predict(m.w.16sp.surv4, newdata= new.w.16sp.surv , re.form=~0)

pred.w.16sp.surv <- within(new.w.16sp.surv, {
  Predprob <- plogis(survpred)
  LL <- plogis(lci)
  UL <- plogis(uci)
})

ggplot(pred.w.16sp.surv, aes(x = exclosure, y = Predprob, ymin = LL, ymax = UL,
                             color=as.factor(treatment))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) + theme_bw() +
  ylab("Survival")


### treatment
new.w.16sp.surv1 <- expand.grid(treatment = c("C", "F","I","W","S"))

bootfit.w.16sp.surv1 <- bootMer(m.w.16sp.surv6, 
                                FUN=function(x)predict(x, new.w.16sp.surv1, re.form=~0), nsim=999)

new.w.16sp.surv1$lci <- apply(bootfit.w.16sp.surv1$t, 2, quantile, 0.025) 
new.w.16sp.surv1$uci <- apply(bootfit.w.16sp.surv1$t, 2, quantile, 0.975) 
new.w.16sp.surv1$survpred <- predict(m.w.16sp.surv6, newdata= new.w.16sp.surv1, re.form=~0)

pred.w.16sp.surv1 <- within(new.w.16sp.surv1, {
  Predprob <- plogis(survpred)
  LL <- plogis(lci)
  UL <- plogis(uci)
})

ggplot(pred.w.16sp.surv1, aes(x = treatment, y = Predprob, ymin = LL, ymax = UL)) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) + theme_bw() +
  ylab("Survival")

### tree species
### both treatment and exclosure
new.w.16sp.tree.surv <- expand.grid(treatment = c("C", "F","I","W","S"), exclosure = c("0","1"))

bootfit.w.16sp.tree.surv <- bootMer(m.w.16sp.tree.surv, 
                                    FUN=function(x)predict(x, new.w.16sp.tree.surv, re.form=~0), nsim=999)

new.w.16sp.tree.surv$lci <- apply(bootfit.w.16sp.tree.surv$t, 2, quantile, 0.025) 
new.w.16sp.tree.surv$uci <- apply(bootfit.w.16sp.tree.surv$t, 2, quantile, 0.975) 
new.w.16sp.tree.surv$survpred <- predict(m.w.16sp.tree.surv, newdata= new.w.16sp.tree.surv , re.form=~0)

pred.w.16sp.tree.surv <- within(new.w.16sp.tree.surv, {
  Predprob <- plogis(survpred)
  LL <- plogis(lci)
  UL <- plogis(uci)
})

ggplot(pred.w.16sp.tree.surv, aes(x = exclosure, y = Predprob, ymin = LL, ymax = UL,
                                  color=as.factor(treatment))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) + theme_bw() +
  ylab("Survival")









