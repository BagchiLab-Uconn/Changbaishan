## Data set1
## Pesticide treatments analysis
## period: form 2015 full to 2016 full, three census (two census periods)
## Initial analysis

## Library packages
library(ggplot2)
library(Matrix)
library(lme4)
library(permute)
library(lattice)
library(vegan)
library(tidyr)
library(dplyr)
library(pbkrtest)
library(verification)

## load the data
pest <- read.csv("Changbaishan/data/Pestwoodyseedlings 20170116.csv")
snow <- read.csv("Changbaishan/data/Snowwoodyseedlings 20170115.csv")

## Just include F, I and W
# pest.new <- read.csv("Changbaishan/data/Pestwoodyseedlings (F,I & W) 20170116.csv")

#####################
## data preparation #
#####################

## gather data
## put all the height of each census into one column
## add one column, and put all survival situation into this column

## First comes to the pest data
## use tidyr::gather() to get longer data with all height and surv in one column
gathered.pest1 <- gather(pest, census_surv, heig_surv, height_15f:surv_16f)
## use separate() to separate the height and surv in two columns
gathered.pest2 <- gathered.pest1 %>%
  separate(census_surv,c("item","census"))
## use spread() to get the final version with height and surv in separate columns
gathered.pest <- gathered.pest2 %>%
  spread(item,heig_surv)
## the data set (gathered.pest) is the data that we will analyze

## Similar for the snow data set
gathered.snow1 <- gather(snow, census_surv, heig_surv, height_15f:surv_16f)
gathered.snow2 <- gathered.snow1 %>%
  separate(census_surv,c("item","census"))
gathered.snow <- gathered.snow2 %>%
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
gathered.pest$site <- factor(gathered.pest$site)
gathered.pest$quadrat <- factor(gathered.pest$quadrat)
gathered.pest$tag <- factor(gathered.pest$tag)
gathered.pest$exclosure <- factor(gathered.pest$exclosure)
gathered.pest$census <-  as.factor(gathered.pest$census)
gathered.pest$quad1 <-  as.factor(gathered.pest$code):gathered.pest$quadrat

## gathered.snow
gathered.snow$site <- factor(gathered.snow$site)
gathered.snow$quadrat <- factor(gathered.snow$quadrat)
gathered.snow$tag <- factor(gathered.snow$tag)
gathered.snow$census <-  as.factor(gathered.snow$census)
gathered.snow$quad1 <-  as.factor(gathered.snow$site):gathered.snow$quadrat

## pest.new
#gathered.pest.new$site <- factor(gathered.pest.new$site)
#gathered.pest.new$code <- factor(gathered.pest.new$code)
#gathered.pest.new$quadrat <- factor(gathered.pest.new$quadrat)
#gathered.pest.new$tag <- factor(gathered.pest.new$tag)
#gathered.pest.new$exclosure <- factor(gathered.pest.new$exclosure)
#gathered.pest.new$quad1 <-  as.factor(gathered.pest.new$code):gathered.pest.new$quadrat

################
### survival ###
################

# survival status 
# gathered.dat$surv <- as.numeric(gathered.dat$height !=0 & !is.na(gathered.dat$height))
# gathered.dat$surv <- as.numeric(gathered.dat$height !=0)

# see the data structure
table(gathered.pest$surv)
table(gathered.pest$species)

table(gathered.snow$surv)

#### set new data set without NAs
gathered.pest0 <- subset(gathered.pest,!is.na(gathered.pest$surv))
gathered.snow0 <- subset(gathered.snow,!is.na(gathered.snow$surv))
# gathered.pest.new0 <- subset(gathered.pest.new,!is.na(gathered.pest.new$surv))

#######################
## mixed-effect model #
#######################
## gathered.pest
m.pest.all.surv <- glmer(surv ~ pesticide*exclosure + census + growth.form +  (1|site/quadrat) + (1|species),
                      data=gathered.pest, family = binomial)
summary(m.pest.all.surv)

## consider the interaction between pesticide and growth.form
m.pest.all.surv1 <- glmer(surv ~ census + pesticide*exclosure + pesticide*growth.form +  (1|site/quadrat) + (1|species),
                       data=gathered.pest, family = binomial)
summary(m.pest.all.surv1)

## consider the interaction among pesticide, exclosure and growth.form
m.pest.all.surv2 <- glmer(surv ~ census + pesticide*exclosure*growth.form +  (1|site/quadrat) + (1|species),
                       data=gathered.pest, family = binomial)
summary(m.pest.all.surv2)

## consider random slope for the pesticide in each species
m.pest.all.surv3 <- glmer(surv ~ pesticide*exclosure + census + growth.form +  (1|site/quadrat) + (1+pesticide|species),
                          data=gathered.pest, family = binomial)
summary(m.pest.all.surv3)

## interaction across pesticide, exclosure and growth.form
m.pest.all.surv4 <- glmer(surv ~ (pesticide + exclosure + growth.form)^2  + census +(1|site/quadrat) + (1+pesticide|species),
                          data=gathered.pest, family = binomial)
summary(m.pest.all.surv4)

## compare the models
anova (m.pest.all.surv, m.pest.all.surv1, m.pest.all.surv2, m.pest.all.surv3, m.pest.all.surv4)
## the m.pest.all.surv is the model we want

#######
# the default intercept is C (control), but in this model, we want to know if there is difference
# between pesticide and water 
# relevel factor, put W (water control) as the intercept
# m.g.all.surv <- profile(glmer(surv ~ pesticide*exclosure + census + growth.form +  (1|site/quadrat) + (1|species),
#                              within(data=gathered.dat, pesticide <- relevel(pesticide, "W")), family = binomial)) 
# summary(m.g.all.surv)

gathered.pest$pesticide <- relevel(gathered.pest$pesticide, ref = "W")
m.pest.all.surv <- glmer(surv ~ pesticide*exclosure + census + growth.form +  (1|site/quadrat) + (1|species),
                      data=gathered.pest, family = binomial)
summary(m.pest.all.surv)


##################
## gathered.snow #
##################
m.snow.all.surv0 <- glmer(surv ~ treatment + census + growth.form +  (1|site/quadrat) + (1|species),
                          data=gathered.snow, family = binomial)
summary(m.snow.all.surv0)

# model includes interaction between treatment and growth.form
m.snow.all.surv <- glmer(surv ~ treatment*growth.form + census + (1|site/quadrat) + (1|species),
                         data=gathered.snow, family = binomial)
summary(m.snow.all.surv)
# this one is better

# consider the random slope
m.snow.all.surv2 <- glmer(surv ~ treatment*growth.form + census + (1|site/quadrat) + (1 + treatment|species),
                          data=gathered.snow, family = binomial)
summary(m.snow.all.surv2)

anova (m.snow.all.surv, m.snow.all.surv1,m.snow.all.surv2)

# m.snow.all.surv is the model we want to use

# gathered.pest.new
# m.pest.new.all.surv <- glmer(surv ~ pesticide*exclosure + census + growth.form +  (1|site/quadrat) + (1|species),
#                      data=gathered.pest.new, family = binomial)
# summary(m.pest.new.all.surv)


###---------------------------------####
## Check the residuals with fitted value
## gathered.pest
test.pest.all <- data.frame(res = resid(m.pest.all.surv), fit = fitted(m.pest.all.surv), 
                      pest = gathered.dat0$pesticide, excl = gathered.dat0$exclosure,
                      cens = gathered.dat0$census, site = gathered.dat0$site,
                      quad = gathered.dat0$quad1,
                      grfo = gathered.dat0$growth.form, spec = gathered.dat0$species)

## gathered.snow
test.snow.all <- data.frame(res = resid(m.snow.all.surv), fit = fitted(m.snow.all.surv), 
                         cens = gathered.snow0$census, site = gathered.snow0$site, 
                         grfo = gathered.snow0$growth.form, spec = gathered.snow0$species,
                         quad = gathered.snow0$quad1, tret = gathered.snow0$treatment)

## gathered.pest
qplot(x=fit, y=res, geom='smooth', data=test.pest.all) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("the whole data (gathered.pest)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.pest.all,facets=~pest) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("among pesticide trt (gathered.pest)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.pest.all,facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("in each census (gathered.pest)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.pest.all, colour=pest, facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("pesticide trt in each census (gathered.pest)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.pest.all, colour=pest, facets=~site) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("pesticide trt in each site (gathered.pest)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.pest.all, colour=pest, facets=~grfo) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("pesticide trt between growthform (gathered.pest)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.pest.all, colour=pest, facets=~excl) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("pesticidetrt between exclosure (gathered.pest)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)
##-------------------------------##

## gathered.snow
qplot(x=fit, y=res, geom='smooth', data=test.snow.all) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("the whole data (gathered.snow)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.snow.all,facets=~tret) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("among snow trt (gathered.snow)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.snow.all,facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("in each census (gathered.snow)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.snow.all, colour=tret, facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("snow trt in each census (gathered.snow)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.snow.all, colour=tret, facets=~site) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("snow trt in each site (gathered.snow)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

qplot(x=fit, y=res, geom='smooth', data=test.snow.all, colour=tret, facets=~grfo) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("snow trt between growthform (gathered.snow)")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() +
  geom_text(aes(label = quad), size = 3, 
            angle = 45, check_overlap = TRUE, hjust=0, vjust=0)
##-------------------------------##

## the AUC (Area Under the Curve of the Receiver Operator Characteristic
roc.area(obs=getME(m.g.all.surv, 'y'), pred=fitted(m.g.all.surv))
roc.area(obs=getME(m.snow.all.surv, 'y'), pred=fitted(m.snow.all.surv))
roc.area(obs=getME(m.pest.new.all.surv, 'y'), pred=fitted(m.pest.new.all.surv))


################
## gathered.pest
## analyze growth.form (trees and shrubs) seperately
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


###########################
## get the BLUPs associated with the various predictors
# observed mean value of for each quadrat
surv.mean <- aggregate(surv ~ site + quadrat + census,
                       data=gathered.dat, FUN=mean)

re.surv <- ranef(m.g.all.surv)
str(re.surv)
re.surv.quad <- re.surv$quad


surv.quad <- cbind(re.surv.quad, surv.mean)
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
## use labels
library("ggrepel")
p +   geom_text(aes(label = quad), size = 3, 
                angle = 45, check_overlap = TRUE, hjust=0, vjust=0)

p + geom_point() + theme_bw() +
  ggrepel::geom_text_repel(aes(label = quad), color = "black", size = 2.5, segment.color = "grey")
# it does not work


##--------------------------##
### Normality of residuals
hist(residuals(m.pesr.all.surv))
qqnorm(residuals(m.pest.all.surv))

hist(residuals(m.snow.all.surv))
qqnorm(residuals(m.snow.all.surv))
##-----------------------##


#############################################
### calculate the mean and standard error ###
#############################################
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



######################################
### the observational survival rate ##
## Averaged by pesticide, exclosurem, census and growth.form
gathered.pest.surv  <- aggregate(surv ~ pesticide + exclosure + census + growth.form, 
                                 data=gathered.pest, FUN=mean)

#######
## The fitted vaule
## Also averaged by pesticide, exclosurem, census and growth.form
fit  <- fitted (m.pest.all.surv)
a1  <- cbind (gathered.pest0,fit)
a2 <- aggregate(fit ~ pesticide + exclosure + census + growth.form, 
                data=a1, FUN=mean)
a2

## The observed and fitted values are similar in each group


######
## Extract the predicted value from the glmer model
## only including the fixed-effect

## First, just see the results for the last census (16f)
new.pest.16f <-  expand.grid(pesticide =c("C", "F", "I", "W"), growth.form= c("tree","shrub"),
                             exclosure=c("0", "1"), census= "16f")
# pred.pest.16f <- predict(m.pest.all.surv, newdata = new.pest.16f, type = "response", re.form=~0)
fit.pest.16f <- predict(m.pest.all.surv, newdata = new.pest.16f, type = "link",re.form=~0)
pred.pest.16f <- cbind (new.pest.16f, fit.pest.16f)

new.pest.16f$quadrat <- gathered.pest$quadrat[1]
new.pest.16f$site <- gathered.pest$site[1]
new.pest.16f$species <- gathered.pest$species[1]

## get the confidence intervals
library(httpuv)
library(merTools)
library(foreach)

predInt.pest.16f <- predictInterval(m.pest.all.surv, newdata = new.pest.16f, n.sims = 999, level = 0.95,
                                    type = "linear.prediction", stat = "mean", which="fixed")
# predInt.pest.16f <- cbind(new.pest.16f, predInt.pest.16f)

predprob.pest.16f <- within(predInt.pest.16f, {
  Predprob <- plogis(fit)
  LL <- plogis(lwr)
  UL <- plogis(upr)
})

pred.pest.16f.dat <- cbind(new.pest.16f, predprob.pest.16f)

## Plot the results, using ggplot2
# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.4) # move them .2 to the left and righ

ggplot(pred.pest.16f.dat, aes(x = pesticide, y = Predprob, ymin = LL, ymax = UL,
                              color=as.factor(exclosure))) + 
  geom_pointrange() +  
  facet_wrap("growth.form") + theme_bw()

#######
## all the results in all census
new.pest <-  expand.grid(pesticide =c("C", "F", "I", "W"), growth.form= c("tree","shrub"),
                         exclosure=c("0", "1"), census= c("16s","16f"))
# pred.pest.16f <- predict(m.pest.all.surv, newdata = new.pest.16f, type = "response", re.form=~0)
fit.pest <- predict(m.pest.all.surv, newdata = new.pest, type = "link",re.form=~0)
pred.pest <- cbind (new.pest, fit.pest)

new.pest$quadrat <- gathered.pest$quadrat[1]
new.pest$site <- gathered.pest$site[1]
new.pest$species <- gathered.pest$species[1]

predInt.pest <- predictInterval(m.pest.all.surv, newdata = new.pest, n.sims = 999, level = 0.95,
                                type = "linear.prediction", stat = "mean", which="fixed")

predprob.pest <- within(predInt.pest, {
  Predprob <- plogis(fit)
  LL <- plogis(lwr)
  UL <- plogis(upr)
})

pred.pest.dat <- cbind(new.pest, predprob.pest)


ggplot(pred.pest.dat, aes(x = pesticide, y = Predprob, ymin = LL, ymax = UL,
                          group=interaction(census, exclosure), col=census, shape=exclosure)) + 
  geom_pointrange() + geom_line(aes(y=Predprob, lty=exclosure), position=pd, size=0.8) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("growth.form") + theme_bw()

ggplot(pred.pest.dat, aes(x = pesticide, y = Predprob, ymin = LL, ymax = UL,
                          group=interaction(census, exclosure), col=census, shape=exclosure)) + 
  geom_pointrange(position=pd) + 
  #geom_line(aes(y=Predprob, lty=exclosure), position=pd, size=0.8) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("growth.form") + theme_bw()


### SNOW DATA ###
new.snow <- expand.grid(treatment = c("S", "C"), growth.form= c("tree","shrub"),
                        census= c("15s","16s","16s","16f"))

new.snow$quadrat <- gathered.snow$quadrat[1]
new.snow$site <- gathered.snow$site[1]
new.snow$species <- gathered.snow$species[1]
pred.snow <- predict(m.snow.all.surv, newdata = new.snow, type = "response",re.form=~0)
pred.snow <- cbind (new.snow, pred.snow)

predInt.snow <- predictInterval(m.snow.all.surv, newdata = new.snow , n.sims = 999, level = 0.95,
                                type = "linear.prediction", stat = "mean", which="fixed")

predprob.snow <- within(predInt.snow, {
  Predprob <- plogis(fit)
  LL <- plogis(lwr)
  UL <- plogis(upr)
})

pred.snow.dat <- cbind(new.snow, predprob.snow)

ggplot(pred.snow.dat, aes(x = treatment, y = Predprob, ymin = LL, ymax = UL,
                          color=as.factor(census))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("growth.form") + theme_bw()


###################
## As we consider the census in the fixed effect, we get two different results in each census
## However, we want to know the overall trend acorss the census
## We need to figure out how to get the overall value
## Let's try the contrasts function, which I don't really understand

contrasts(gathered.pest$census) <- "contr.sum"
contrasts(gathered.pest$census) <- contr.sum(3)
contrasts(gathered.pest$census) <- contr.sum(levels(gathered.pest$census))
options(contrasts = c("contr.sum", "contr.poly"))
contrasts(gathered.pest$census)

# gathered.pest$census <- relevel(gathered.pest$census, ref="15f")

m.pest.all.surv <- glmer(surv ~ pesticide*exclosure + census + growth.form +  (1|site/quadrat) + (1|species),
                         data=gathered.pest, family = binomial)
summary(m.pest.all.surv)
## Warning: contrasts dropped from factor census due to missing levels
## I guess it's the census 15f, because the survival in this census is considered to be NAs
## HERE I CNA NOT USE THE CONTRASTS TO COMBINE THE RESULTS OF TWO CENSUSES


##-----------------------------##
## Other vegetation metric analysis

####################
## no. of species ##
####################
## think about how to select with or without NA
## below is not correct
gathered.dat.spec1 <- subset(gathered.dat,gathered.dat$height !=0 & !is.na(gathered.dat$height))
gathered.dat.spec1$no.species <-as.numeric(gathered.dat.spec1$species)
gathered.dat.spec  <- aggregate(no.species ~ treatment + exclosure + site + quadrat + census, 
                                data=gathered.dat, FUN=function(x)length(unique(x)))
# model
g.spec <- lmer(no.species ~ census + treatment*exclosure +(1|site/quadrat),
               data=gathered.dat.spec)
summary(g.spec)

## residuals vs. fittd values
test.g.spec <- data.frame(res = resid(g.spec), fit = fitted(g.spec), 
                          trt = gathered.dat.spec$treatment, excl = gathered.dat.spec$exclosure,
                          site = gathered.dat.spec$site, cens = gathered.dat.spec$census)

###
qplot(x=fit, y=res, geom='smooth', data=test.g.spec) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("no. of species")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.spec, facets=~trt) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("no. of species")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.spec, colour=trt, facets=~site) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("no. of species in each site")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.g.spec, colour=trt, facets=~cens) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("no. of species in each census")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()