### Survival analysis
### Apr 23, 2018
### By Shihong Jia

## Analysis 1: Survival analysis
## new recruits in the spring census, and then check the status (survive or die) in the fall census



getwd()
# setwd("E:/R/My code/Git/Changbaishan/New Pest Analysis")
# setwd("E:/CBSME/Surv")


### load libraries
library(ggplot2)
library(cowplot)
library(lme4)
library(pbkrtest)
library(tidyr)
library(dplyr)
library(optimx)
library(broom)
library(glmmTMB)



### Data

### read data
dat <- read.csv("data/pest_new.dat20.csv")
str(dat)

dat$pesticide <- factor(dat$pesticide, levels=c('W', 'F', 'I'))
dat$census <- factor(dat$census, levels=c('15fa', '16sp', '16fa', '17sp', '17fa'))
dat$exclosure <- factor(dat$exclosure)
dat$site <- factor(dat$site)
dat$quad.unique <- factor(dat$quad.unique)



# this analysis only consider tree species and the survival of seedlings from spring to fall
pest.surv.dat <- subset(dat, growth.form=='tree' & (census=='15fa'|census=='16fa'|census=='17fa'))
pest.surv.dat <- pest.surv.dat[,-c(1,6:9,17,22:24,26:32)]
summary(pest.surv.dat)




### all tree species pooled
## previous analyses had shown that there was no significant interaction effect 
## between exclosure and other variables, so did the seedling density
mod.surv1 <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + exclosure + lg_A.het + lg_S.con + lg_S.het + 
	 census + site +(1+census|quad.unique) + (1|sp.), data=pest.surv.dat, family=binomial)

## do not consider exclosure in the following analyses
mod.surv2 <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + lg_S.con + lg_S.het + 
	 census + site + (1+census|quad.unique) + (1|sp.), data=pest.surv.dat, family=binomial)

# do not consider census as the random-slope effect
mod.surv3 <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + lg_S.con + lg_S.het + 
	 census + site +(1|quad.unique) + (1|sp.), data=pest.surv.dat, family=binomial)

# do not include census and site as fixed-effect, exclude census in the random-slope effect
mod.surv4 <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + lg_S.con + lg_S.het + 
	 (1|quad.unique) + (1|sp.), data=pest.surv.dat, family=binomial)


anova(mod.surv1, mod.surv2, mod.surv3, mod.surv4)

summary(mod.surv2, correlation=FALSE)

summary(mod.surv3, correlation=FALSE)

# mod.surv2 is the best model according to the AIC values.
# mod.tree.surv41.2 and mod.tree.surv41.3 show very similar results.
# we do not consider exclosure in the further analysis.



## glmmTMB model?
modTMB.surv2 <- glmmTMB(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + lg_S.con + lg_S.het + 
                     census + site + (1+census|quad.unique) + (1|sp.), 
                     data=pest.surv.dat, family=betabinomial)

summary(modTMB.surv2, correlation=FALSE)


# the glmmTMB shows similar results compared with glmer.





##---------------- Prediction ----------------------##
## In the previous analysis, we found that there are some big outliers of A.con
hist(pest.surv.dat$A.con)
table(subset(pest.surv.dat, A.con > 3000)$pesticide)

# there are 42 individuals which have the conspecific adult basal area larger than 3000 m2
# there is no seedling in fungicide treatment.
# in our prediction, we would like to see the relationship with the common x sclae.

pest.surv.dat1 <- subset(pest.surv.dat, A.con < 3000)
dim(pest.surv.dat1); 



# re-fit the model
mod.surv2_rf <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + lg_S.con + lg_S.het + 
	 census + site + (1+census|quad.unique) + (1|sp.), data=pest.surv.dat1, family=binomial)
summary(mod.surv2_rf, correlation=FALSE)

# excluding the bigest adults improve the model fitting




##--------------------- bootstrap -------------------##
surv.newdat1 <-  with(pest.surv.dat1, 
		expand.grid(lg_A.con=seq(min(lg_A.con), max(lg_A.con), length=100),
		pesticide=c("W", "F", "I"), 
		lg_height0=mean(lg_height0), lg_A.het=mean(lg_A.het),
		lg_S.con=mean(lg_S.con), lg_S.het=mean(lg_S.het)))
dim(surv.newdat1)


## prediction
surv.newdat1$mean.pred <- predict(mod.surv2_rf, surv.newdat1, re.form=~0, type='response')

# not work
# because we also include census and site as fixed-effect in the model, 
# but we are no interested in these two variables in our prediction process


# one possible solution is exclude the census and site as the fixed-effect and refit a new model
# alternatively, we can adopt a new method to generate the predicted values.




##---------------- lappy approach ------------------------##

# First, let's start with the mod.surv2_tf


## https://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/
## test the pesticide treatments together
# temporary data 
surv.tmpdat1 <- pest.surv.dat1[, c("lg_height0", "pesticide", "lg_A.con", "lg_A.het",
		 "lg_S.con", "lg_S.het", "census", "site", "quad.unique", "sp.")]
dim(surv.tmpdat1)


surv.jvalues1 <- with(pest.surv.dat1, seq(from = min(lg_A.con), to = max(lg_A.con), length.out = 100))


# calculate predicted probabilities and store in a list 
surv.biprobs1 <- lapply(levels(pest.surv.dat1$pesticide), function(treat) 
	{surv.tmpdat1$pesticide[] <- treat 
	lapply(surv.jvalues1, function(j) {
	surv.tmpdat1$lg_A.con <- j 
	# predict(mod.surv2_rf, newdata = surv.tmpdat1, type = "response", re.form=~0) 
	# include the random-effect
	predict(mod.surv2_rf, newdata = surv.tmpdat1, type = "response") 
  }) 
})


# meet an error
# I suppose this may due to the random-slope effect
# Although I can use re.form=~0 in the prediction, I should consider the random-effect in the predition
# because I noticed that the one without considering random-effect performs poorly in the recruitment
# (please check this in the recruitment analysis)


# are there any solutions?
# I get an idea that remove the random-slope from the model given the results are very similar
mod.surv3_rf <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + lg_S.con + lg_S.het + 
                        census + site + (1|quad.unique) + (1|sp.), data=pest.surv.dat1, family=binomial)
summary(mod.surv3_rf, correlation=FALSE)


surv.biprobs1 <- lapply(levels(pest.surv.dat1$pesticide), function(treat) {
  surv.tmpdat1$pesticide[] <- treat 
  lapply(surv.jvalues1, function(j) {
    surv.tmpdat1$lg_A.con <- j 
    predict(mod.surv3_rf, newdata = surv.tmpdat1, type = "response") 
  }) 
})



# get means and quartiles for all jvalues for each level of Pesticide 
surv.plotdat1 <- lapply(surv.biprobs1, function(X) {
	 surv.temp <- t(sapply(X, function(x) {
	 c(M=mean(x), quantile(x, c(.25, .75))) 
	}))
	surv.temp <- as.data.frame(cbind(surv.temp, surv.jvalues1))
	colnames(surv.temp) <- c("Predicted", "Lower", "Upper", "Lg_A.con") 
	return(surv.temp)
})


# collapse to one data frame 
surv.plotdat1 <- do.call(rbind, surv.plotdat1) 

# add pesticide treatment accordingly
surv.plotdat1$pesticide <- factor(rep(levels(pest.surv.dat1$pesticide), each = length(surv.jvalues1)))

# show first few rows 
head(surv.plotdat1)



# back to the original data of x
mean_logx <- mean(log(pest.surv.dat1$A.con + 1))
sd_logx <- sd(log(pest.surv.dat1$A.con + 1))
surv.plotdat1$x <- exp((surv.plotdat1$Lg_A.con*sd_logx) + mean_logx)-1

# relevel the pesticide treatment
surv.plotdat1$pesticide <- factor(surv.plotdat1$pesticide, levels=c("W","F","I"))
summary(surv.plotdat1)


p.surv.pred3 <- ggplot() +
  geom_line(data=surv.plotdat1, aes(x=x, y=Predicted, colour=pesticide),size=0.8) +
  geom_ribbon(data=surv.plotdat1, aes(x=x, ymin = Lower, ymax = Upper, fill = pesticide), alpha = 0.2) +
  scale_color_manual(values = c('gray40','orangered','limegreen'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','limegreen'),
	labels=c("Control","Fungicide","Insecticide")) +
  labs(x="Conspecific adult density (m2)", y="Survival", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide") +
  theme(legend.title=element_blank(), legend.background = element_rect(fill="lightgray"),
		legend.position ="right") + theme_bw()+
  ylim(c(0.5,1.0))





##--------------------- the observed values ----------------------##

## keep the lg_A.con vary with the dataset but other variables constant

## generate k sets of simulated y at these new points.
ob.surv.newdat <- with(pest.surv.dat1, 
                  expand.grid(lg_A.con=lg_A.con, pesticide=levels(pesticide), 
                              lg_height0=mean(lg_height0), lg_A.het=mean(lg_A.het),
                              lg_S.con=mean(lg_S.con), lg_S.het=mean(lg_S.het)))

dim(ob.surv.newdat)


# use re.form=~0 and on link scale
ob.surv.newdat$pred <- predict(mod.surv3_rf, newdata=ob.surv.newdat, re.form=~0, type='link')

# not work if I include the census and site in the model

# It's a little tedious now, I don't want to confuse you.
# I just want to show you more details
# Together above, I think the model should be mod.surv4
# the one without considering census and site as fixed-effect and remove census from the random-effect too.
# I hope this time we will get there





##------------------ FINAL SHOOT --------------------##
mod.surv4_rf <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + lg_S.con + lg_S.het + 
                     (1|quad.unique) + (1|sp.), data=pest.surv.dat1, family=binomial)


# prediction
surv.biprobs4 <- lapply(levels(pest.surv.dat1$pesticide), function(treat) {
  surv.tmpdat1$pesticide[] <- treat 
  lapply(surv.jvalues1, function(j) {
    surv.tmpdat1$lg_A.con <- j 
    predict(mod.surv4_rf, newdata = surv.tmpdat1, type = "response") 
  }) 
})



# get means and quartiles for all jvalues for each level of Pesticide 
surv.plotdat4 <- lapply(surv.biprobs4, function(X) {
  surv.temp <- t(sapply(X, function(x) {
    c(M=mean(x), quantile(x, c(.25, .75))) 
  }))
  surv.temp <- as.data.frame(cbind(surv.temp, surv.jvalues1))
  colnames(surv.temp) <- c("Predicted", "Lower", "Upper", "Lg_A.con") 
  return(surv.temp)
})


# collapse to one data frame 
surv.plotdat4 <- do.call(rbind, surv.plotdat4) 

# add pesticide treatment accordingly
surv.plotdat4$pesticide <- factor(rep(levels(pest.surv.dat1$pesticide), each = length(surv.jvalues1)))



# back to the original data of x
surv.plotdat4$x <- exp((surv.plotdat4$Lg_A.con*sd_logx) + mean_logx)-1

# relevel the pesticide treatment
surv.plotdat4$pesticide <- factor(surv.plotdat4$pesticide, levels=c("W","F","I"))
summary(surv.plotdat4)


p.surv.pred4 <- ggplot() +
  geom_line(data=surv.plotdat4, aes(x=x, y=Predicted, colour=pesticide),size=0.8) +
  geom_ribbon(data=surv.plotdat4, aes(x=x, ymin = Lower, ymax = Upper, fill = pesticide), alpha = 0.2) +
  scale_color_manual(values = c('gray40','orangered','limegreen'),
                     labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','limegreen'),
                    labels=c("Control","Fungicide","Insecticide")) +
  labs(x="Conspecific adult density (m2)", y="Survival", 
       colour="Pesticide",linetype="Pesticide", fill="Pesticide") +
  theme(legend.title=element_blank(), legend.background = element_rect(fill="lightgray"),
        legend.position ="right") + theme_bw()+
  ylim(c(0.5,1.0))


# how about the one with census and site as fixed-effect
p.surv.pred3
# great, they are so similar





#### Observed values

# use re.form=~0 and on link scale
ob.surv.newdat$pred <- predict(mod.surv4_rf, newdata=ob.surv.newdat, re.form=~0, type='link')


# set a new dataset for ploting
p.pest.surv.dat1 <- pest.surv.dat1


# calculate the residuals from the model
p.pest.surv.dat1$resd <- resid(mod.surv4_rf)


# give the residuals to the new data frame according to the lg_A.con values
ob.surv.newdat$resd <- p.pest.surv.dat1[p.pest.surv.dat1$lg_A.con %in% ob.surv.newdat$lg_A.con,]$resd


# the observed values are generated by adding the residuals to the predicted values
ob.surv.newdat$obse_lk <- ob.surv.newdat$pred + ob.surv.newdat$resd

# back transform the observed values
ob.surv.newdat$obse <- plogis(ob.surv.newdat$obse_lk)


# replace the pesticide
ob.surv.newdat$pesticide <- p.pest.surv.dat1[p.pest.surv.dat1$lg_A.con %in% ob.surv.newdat$lg_A.con,]$pesticide


# back to the original data of x
ob.surv.newdat$x <- exp((ob.surv.newdat$lg_A.con*sd_logx) + mean_logx)-1
summary(ob.surv.newdat)


ggplot() + 
  geom_line(data=surv.plotdat4, aes(x = x, y = Predicted, color=pesticide),size=1) +
  geom_ribbon(data=surv.plotdat4, aes(x=x, ymin = Lower, ymax = Upper, fill = pesticide), alpha = 0.2)+
  geom_point(data=ob.surv.newdat,aes(x=x, y = obse, color=pesticide)) +
  scale_color_manual(values = c('gray40','orangered','darkblue'),
                     labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','darkblue'),
                    labels=c("Control","Fungicide","Insecticide")) +
  labs(x="Conspecific adult density (m2)", y="Survival", 
       colour="Pesticide",linetype="Pesticide", fill="Pesticide")



## bin the x to 300 density
range(ob.surv.newdat$x)
ob.surv.newdat$A.con_cat <- cut(ob.surv.newdat$x, seq(0, 3000, 250))

obse.surv_cat <- summarise(group_by(ob.surv.newdat, A.con_cat, pesticide), 
                           mean = mean(obse), se=sd(obse)/sqrt(length(obse)))
summary(obse.surv_cat)

# change the A.con_cat to continuous data
obse.surv_cat$x <- rep(seq(125,2875,250), each=3)


p.ob.pr.surv <- ggplot() +
  geom_line(data=surv.plotdat4, aes(x=x, y=Predicted, colour=pesticide),size=0.8) +
  geom_ribbon(data=surv.plotdat4, aes(x=x, ymin = Lower, ymax = Upper, fill = pesticide), alpha = 0.2) +
  geom_point(data=obse.surv_cat,aes(x=x, y = mean, color=pesticide),
             size=2,position=position_dodge(width = 115)) +
  geom_errorbar(data=obse.surv_cat,aes(x=x, ymin=mean-se, ymax=mean+se,color=pesticide),
                width=0.1, position=position_dodge(width = 115)) +
  scale_color_manual(values = c('gray40','orangered','limegreen'),
                     labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','limegreen'),
                    labels=c("Control","Fungicide","Insecticide")) +
  labs(x="Conspecific adult density (m2)", y="Survival", 
       colour="Pesticide",linetype="Pesticide", fill="Pesticide") +
  theme(legend.title=element_blank(), legend.background = element_rect(fill="lightgray"),
        legend.position ="right") + theme_bw()+
  ylim(c(0.5,1.0))


