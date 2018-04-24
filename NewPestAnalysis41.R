### Survival analysis
### Apr 19, 2018
### By Shihong Jia

## Analysis 1: Survival analysis
## use survival from spring to fall
## 15 (no 107 plot) 16 and 17 census
## Extract coefficients
## Calculate conditional R square
## individual species > 80 individuals [old]
## most common and less common[new]
## plot regeneration and survival together

## Revisit the conspecific neighboring adults issue
## All seedlings
## Add site as fixed-effect


## NewPestAnalysis30: prediction with both fixed and random-effect factors
## NewPestAnalysis31: prediction with only fixed effect factors
## NewPestAnalysis35: A.con > 0, glmmTMB and prediction
## Redo the simulated regression figures (back transform, 
## replace the observed line with points)
## include the seedling densities
## new method for prediction



getwd()
setwd("E:/R/My code/Git/Changbaishan/New Pest Analysis")
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


### Analysis One

## pest_new.dat20: new recruited seedlings, check the survival status in the following census

### read data
pest_new.dat20 <- read.csv("data/pest_new.dat20.csv")

pest_new.dat20$pesticide <- factor(pest_new.dat20$pesticide, levels=c('W', 'F', 'I'))
pest_new.dat20$census <- factor(pest_new.dat20$census, levels=c('15fa', '16sp', '16fa', '17sp', '17fa'))
pest_new.dat20$exclosure <- factor(pest_new.dat20$exclosure)
pest_new.dat20$site <- factor(pest_new.dat20$site)


## season
unique(pest_new.dat20$census)
pest_new.dat20$season <- ifelse(pest_new.dat20$census=='16sp'| pest_new.dat20$census=='17sp', 'sp', 'fa')
pest_new.dat20$season <- factor(pest_new.dat20$season)
table(pest_new.dat20$season)
table(pest_new.dat20$growth.form)
table(pest_new.dat20$growth.form, pest_new.dat20$season)
table(pest_new.dat20$census)
table(pest_new.dat20[pest_new.dat20$growth.form=='tree',]$sp.)


# initial height
hist(pest_new.dat20$height0)
hist(log(pest_new.dat20$height0))
pest_new.dat20$height0_log <- scale(log(pest_new.dat20$height0))
summary(pest_new.dat20)

pest_new.dat20$A.con_log <- as.numeric(pest_new.dat20$A.con_log)
pest_new.dat20$A.het_log <- as.numeric(pest_new.dat20$A.het_log)
pest_new.dat20$A.sum_log <- as.numeric(pest_new.dat20$A.sum_log)
pest_new.dat20$con.seedling_log <- as.numeric(pest_new.dat20$con.seedling_log)
pest_new.dat20$het.seedling_log <- as.numeric(pest_new.dat20$het.seedling_log)
pest_new.dat20$height0_log <- as.numeric(pest_new.dat20$height0_log)
str(pest_new.dat20)




## survival 
## Since bootMer doesn't work well with glmerControl, 
# we need to drop the glmerControl, and dropping this doesn't change the model result


#### data
pest.new.surv.dat27 <- subset(pest_new.dat20, growth.form=='tree' & season=='fa')
pest.new.surv.dat27 <- pest.new.surv.dat27[,-c(5:8,17:18,20:26,32)]
summary(pest.new.surv.dat27)
write.csv(pest.new.surv.dat27, "results/pest.new.surv.dat.csv")


## A.con distribution
hist(pest.new.surv.dat27$lg_A.con)
hist(subset(pest.new.surv.dat27, A.con>0)$lg_A.con)
dim(pest.new.surv.dat27); dim(subset(pest.new.surv.dat27, A.con>0))
table(pest.new.surv.dat27$sp.)



## all tree species pooled
## include all variables (i.e. exclosure and seedling density)
## previous analyses had shown that there was no significant interaction effect 
## between exclosure and other variables, so did the seedling density
mod.tree.surv41.1 <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + exclosure + lg_A.het + lg_S.con + lg_S.het + 
	 census + site +(1+census|quad.unique) + (1|sp.), data=pest.new.surv.dat27, family=binomial)
summary(mod.tree.surv41.1, correlation=FALSE)


## do not consider exclosure in the following analyses
mod.tree.surv41.2 <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + lg_S.con + lg_S.het + 
	 census + site + (1+census|quad.unique) + (1|sp.), data=pest.new.surv.dat27, family=binomial)
summary(mod.tree.surv41.2, correlation=FALSE)


# do not consider census as the random-slope effect
mod.tree.surv41.3 <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + lg_S.con + lg_S.het + 
	 census + site +(1|quad.unique) + (1|sp.), data=pest.new.surv.dat27, family=binomial)
summary(mod.tree.surv41.3, correlation=FALSE)

# do not include census and site as fixed-effect, exclude census in the random-slope effect
mod.tree.surv41.4 <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + lg_S.con + lg_S.het + 
	 (1|quad.unique) + (1|sp.), data=pest.new.surv.dat27, family=binomial)
summary(mod.tree.surv41.4, correlation=FALSE)


anova(mod.tree.surv41.1, mod.tree.surv41.2, mod.tree.surv41.3, mod.tree.surv41.4)


# include exclosure do not change the result, we do not consider it further
# mod.tree.surv41.2 and mod.tree.surv41.2 show very similar results, 
# althought mod.tree.surv41.2 is better according to the AIC values.





## remove some big outliers of A.con
a <-  subset(pest.new.surv.dat27, A.con >= 3000)
range(a$A.con)
pest.new.surv.dat40 <- subset(pest.new.surv.dat27, A.con < 3000)
dim(pest.new.surv.dat40); dim(pest.new.surv.dat27)
# remove 42 individuals totally


# re-fit the model
mod.tree.surv40 <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + lg_S.con + lg_S.het + 
	 census + site + (1|quad.unique) + (1|sp.), data=pest.new.surv.dat40, family=binomial)
summary(mod.tree.surv40, correlation=FALSE)

# exclude the big trees improve the model fitting


# glmmTMB model
modTMB.tree.surv40 <- glmmTMB(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + lg_S.con + lg_S.het + 
	 census + site + (1|quad.unique) + (1|sp.), data=pest.new.surv.dat40, family=betabinomial)
summary(modTMB.tree.surv40, correlation=FALSE)



# exclude census and site in the model
mod.tree.surv41 <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + lg_S.con + lg_S.het + 
	 (1|quad.unique) + (1|sp.), data=pest.new.surv.dat40, family=binomial)
summary(mod.tree.surv41, correlation=FALSE)





#################

### calculate the confidence intervals of the predicted values

#################

p.tree.surv.dat41 <- pest.new.surv.dat40


##--------------------- bootstrap -------------------##
surv.newdat41 <-  with(p.tree.surv.dat41, 
		expand.grid(lg_A.con=seq(min(lg_A.con), max(lg_A.con), length=100),
		pesticide=c("W", "F", "I"), 
		lg_height0=mean(lg_height0), lg_A.het=mean(lg_A.het),
		lg_S.con=mean(lg_S.con), lg_S.het=mean(lg_S.het)))


## prediction
## 95% CI of the predicted values
## bootstrap
boot.surv41 <- bootMer(mod.tree.surv41, FUN=function(x)
		predict(x, surv.newdat41, re.form=~0, type='response'), nsim=99)

surv.newdat41$lci.pred <- apply(boot.surv41$t, 2, quantile, probs = 0.025) 
surv.newdat41$uci.pred <- apply(boot.surv41$t, 2, quantile, probs = 0.975)
# surv.newdat41$mean.pred <- apply(boot.surv41$t, 2, quantile, probs = 0.5)
surv.newdat41$mean.pred <- predict(mod.tree.surv41, surv.newdat41, re.form=~0, type='response')



# back to the original data of x
p.tree.surv.dat41$logx_s <- scale(log(p.tree.surv.dat41$A.con + 1))
mean_logx <- mean(log(p.tree.surv.dat41$A.con + 1))
sd_logx <- sd(log(p.tree.surv.dat41$A.con + 1))
surv.newdat41$x <- exp((surv.newdat41$lg_A.con*sd_logx) + mean_logx)-1
summary(surv.newdat41)



# log and then standardize the x 
ggplot(surv.newdat41) +
  geom_line(aes(x=lg_A.con, y=mean.pred, colour =pesticide),size=0.8) +
  scale_color_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  labs(x="Conspecific adult density (m2)", y="Survival", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")



# the original scale with 95% CI
ggplot(surv.newdat41) +
  geom_line(aes(x=x, y=mean.pred, colour=pesticide),size=0.8) +
  geom_ribbon(aes(x=x, ymin = lci.pred, ymax = uci.pred, fill = pesticide), alpha = 0.3) +
  scale_color_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  labs(x="Conspecific adult density (m2)", y="Survival", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")







#### add the observed values
## generate k sets of simulated y at these new points.
ob.newdat <- with(p.tree.surv.dat40, 
		expand.grid(lg_A.con=lg_A.con, pesticide=levels(pesticide), 
		lg_height0=mean(lg_height0), lg_A.het=mean(lg_A.het),
		lg_S.con=mean(lg_S.con), lg_S.het=mean(lg_S.het)))

dim(ob.newdat)

ob.newdat$pred <- predict(mod.tree.surv41, newdata=ob.newdat, re.form=~0, type='link')
p.tree.surv.dat41$resd <- resid(mod.tree.surv41)
str(p.tree.surv.dat41)

# select the lg_A.con in the original data
ob.newdat$resd <- p.tree.surv.dat41[p.tree.surv.dat41$lg_A.con %in% ob.newdat$lg_A.con,]$resd
summary(ob.newdat$resd)

ob.newdat$obse_lk <- ob.newdat$pred + ob.newdat$resd
# ob.newdat$obse_exp <- exp(ob.newdat$obse_lk)
ob.newdat$obse <- plogis(ob.newdat$obse_lk)

# another method to transform data
# require(boot)
# ob.newdat$obse <- inv.logit(ob.newdat$obse_lk)
summary(ob.newdat)



# replace the pesticide
ob.newdat$Pesticide <- p.tree.surv.dat41[p.tree.surv.dat41$lg_A.con %in% ob.newdat$lg_A.con,]$pesticide


# back to the original data of x
ob.newdat$x <- exp((ob.newdat$lg_A.con*sd_logx) + mean_logx)-1
summary(ob.newdat); dim(ob.newdat)


ggplot() + 
	geom_line(data=surv.plotdat2, aes(x = x, y = PredictedProbability, color=Pesticide),size=1) +
	geom_ribbon(data=surv.plotdat2, aes(x=x, ymin = Lower, ymax = Upper, fill = Pesticide), alpha = 0.2)+
	geom_point(data=ob.newdat,aes(x=x, y = obse, color=Pesticide)) +
	scale_color_manual(values = c('gray40','orangered','darkblue'),
		labels=c("Control","Fungicide","Insecticide")) + 
  	scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  		labs(x="Conspecific adult density (m2)", y="Survival", 
  	colour="Pesticide",linetype="Pesticide", fill="Pesticide")




## bin the x to 300 density
range(ob.newdat$x)
ob.newdat$A.con_cat <- cut(ob.newdat$x, seq(0, 3000, 250))

obse.surv_cat <- summarise(group_by(ob.newdat, A.con_cat, Pesticide), 
			mean = mean(obse), se=sd(obse)/sqrt(length(obse)))
summary(obse.surv_cat)

# change the A.con_cat to continuous data
obse.surv_cat$x <- rep(seq(125,2875,250), each=3)


p.ob.pr.surv <- ggplot(surv.newdat41) +
  geom_line(aes(x=x, y=mean.pred, colour=pesticide),size=0.8) +
  geom_ribbon(aes(x=x, ymin = lci.pred, ymax = uci.pred, fill = pesticide), alpha = 0.2) +
  geom_point(data=obse.surv_cat,aes(x=x, y = mean, color=Pesticide),
		size=2,position=position_dodge(width = 115)) +
  geom_errorbar(data=obse.surv_cat,aes(x=x, ymin=mean-1.96*se, ymax=mean+1.96*se,color=Pesticide),
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

ggsave(plot=p.ob.pr.surv, 
	filename="results/fig.cbs20180421/Fig.1-survival-prediction&observation-overall.pdf", 
	dpi=300,width=8,height=6,useDingbats=FALSE)









##---------------- lappy approach ------------------------##
## this method can include census and site as fixed-effect, also census as a random slope

# try to include all seedlings
mod.tree.surv41.1 <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + 
	 lg_S.con + lg_S.het + census + site + (1+census|quad.unique) + (1|sp.),
	 data=pest.new.surv.dat27, family=binomial)
summary(mod.tree.surv41.1, correlation=FALSE)

# not work for this model


## https://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/
## test the pesticide treatments together
# temporary data 
surv.tmpdat41 <- pest.new.surv.dat27[, c("lg_height0", "pesticide", "lg_A.con", "lg_A.het",
		 "lg_S.con", "lg_S.het", "census", "site", "quad.unique", "sp.")]
dim(surv.tmpdat41)


surv.jvalues41 <- with(pest.new.surv.dat27, seq(from = min(lg_A.con), to = max(lg_A.con), length.out = 100))


# calculate predicted probabilities and store in a list 
surv.biprobs41 <- lapply(levels(pest.new.surv.dat27$pesticide), function(treat) 
	{surv.tmpdat41$pesticide[] <- treat 
	lapply(surv.jvalues41, function(j) {
	surv.tmpdat41$lg_A.con <- j 
	# predict(mod.tree.surv41.1, newdata = surv.tmpdat41, type = "response", re.form=~0) 
	# include the random-effects
	predict(mod.tree.surv41.1, newdata = surv.tmpdat41, type = "response") 
  }) 
})


# get means and quartiles for all jvalues for each level of Pesticide 
surv.plotdat2 <- lapply(surv.biprobs41, function(X) {
	 surv.temp <- t(sapply(X, function(x) {
	 c(M=mean(x), quantile(x, c(.25, .75))) 
	}))
	surv.temp <- as.data.frame(cbind(surv.temp, surv.jvalues))
	colnames(surv.temp) <- c("Predicted", "Lower", "Upper", "Lg_A.con_s") 
	return(surv.temp)
})


# collapse to one data frame 
surv.plotdat2 <- do.call(rbind, surv.plotdat2) 

# add pesticide treatment
surv.plotdat2$pesticide <- factor(rep(levels(p.tree.surv.dat40$pesticide), each = length(surv.jvalues)))

# show first few rows 
head(surv.plotdat2)


# back to the original data of x
p.tree.surv.dat40$logx_s <- scale(log(p.tree.surv.dat40$A.con + 1))
mean_logx <- mean(log(p.tree.surv.dat40$A.con + 1))
sd_logx <- sd(log(p.tree.surv.dat40$A.con + 1))
surv.plotdat2$x <- exp((surv.plotdat2$Lg_A.con_s*sd_logx) + mean_logx)-1
surv.plotdat2$pesticide <- factor(surv.plotdat2$pesticide, levels=c("W","F","I"))
summary(surv.plotdat2)


ggplot() +
  geom_line(data=surv.plotdat2, aes(x=x, y=Predicted, colour=pesticide),size=0.8) +
  geom_ribbon(data=surv.plotdat2, aes(x=x, ymin = Lower, ymax = Upper, fill = pesticide), alpha = 0.2) +
  geom_point(data=obse.surv_cat,aes(x=x, y = mean, color=Pesticide),
		size=2,position=position_dodge(width = 115)) +
  geom_errorbar(data=obse.surv_cat,aes(x=x, ymin=mean-se, ymax=mean+se,color=Pesticide),
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







###----------------- Individual species analysis ---------------------###
# For individual species, to fit the model better, we exclude those without conspecific neighboring adult


dim(pest.new.surv.dat27); dim(pest.new.surv.dat27.1)
# Exclude 349 seedlings in total

## how many individuals in each species
ind.sp.pest.surv <- table(pest.new.surv.dat27$sp.)
# write.csv(ind.sp.pest.surv, "results/survival-individual in each species.csv")

# combing two censuses with the same season
sp.surv.counts <- as.data.frame(table(pest.new.surv.dat27$sp., pest.new.surv.dat27$pesticide))
colnames(sp.surv.counts) <- c('sp.', 'pesticide', 'no.quad')
sp.surv.counts <- sp.surv.counts %>% spread(pesticide, no.quad) 

# species selected
# at least 10 quadrats in each treatment
sel.sp.surv <- unique(sp.surv.counts$sp.[sp.surv.counts$W >= 10 & 
                                                 (sp.surv.counts$F >= 10 |
                                                    sp.surv.counts$I >= 10)])
# just include seedlings more than 100?
# The confidence intervals are too wide in ACMO
# sp.surv.sel2 <- c('ABNE', 'ACBA','ACMO', 'ACPS', 'ACTE', 'FRMA', 'PIKO', 'QUMO', 'TIAM')
# sp.surv.sel4 <- c('ACBA', 'ACPS', 'FRMA', 'TIAM')

sp.surv.sel2 <- c('ABNE', 'ACBA', 'ACPS', 'ACTE', 'FRMA', 'TIAM')


## glmmTMB model

## remove seedling density
surv_mods41 <- sapply(sp.surv.sel2, function(sp){
  print(sp)
  spdat <- filter(pest.new.surv.dat40, sp.==sp)
  spdat <- droplevels(spdat)
  
  # model
  mod <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + (1|quad.unique), 
               data=spdat, family=binomial)
}, simplify=FALSE)

lapply(surv_mods41, function(x) summary(x, correlation=FALSE))



##################################
# predict with the glmmTMB model results

## prediction
surv.mods41 <- sapply(sp.surv.sel2, function(sp){
  print(sp)
  spdat <- filter(pest.new.surv.dat40,sp.==sp)
  spdat <- droplevels(spdat)
  
  # model
  mod <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + (1|quad.unique), 
               data=spdat, family=binomial)
  p.surv.dat <- spdat
  
  # calculate the confidence intervals
  boot.surv.dat <- data.frame(p.surv.dat, 
	predict(mod, interval='confidence', type='link', se.fit=TRUE))
  
  # results
  p.surv.dat
  
}, simplify=FALSE)


###----------------------------
## plot species with more than 80 individuals
sel.pred.surv.dat36 <- do.call(rbind.data.frame, surv.TMBmods36)
str(sel.pred.surv.dat36)


## plot one species by one species
## ABNE
ABNE.surv.dat36 <- subset(sel.pred.surv.dat36, sp.=='ABNE')

p.ABNE.surv36 <- ggplot(data=ABNE.surv.dat36) +
  geom_smooth(method='lm', se = FALSE, 
              aes(x=lg_A.con, y=surv1, colour =pesticide, linetype=pesticide)) +
  scale_color_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  scale_linetype_manual(values = c("dashed","dashed","dashed"),
	labels=c("Control","Fungicide","Insecticide")) +
  geom_smooth(method='lm', 
              aes(x=lg_A.con, y=fit, colour=pesticide, fill=pesticide),alpha=0.4) +
  labs(x="", y="", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")+
  coord_cartesian(ylim=c(0.1,1.0))   


## ACBA
ACBA.surv.dat36 <- subset(sel.pred.surv.dat36, sp.=='ACBA')

p.ACBA.surv36 <- ggplot(data=ACBA.surv.dat36) +
  geom_smooth(method='lm', se = FALSE, 
              aes(x=lg_A.con, y=surv1, colour =pesticide, linetype=pesticide)) +
  scale_color_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  scale_linetype_manual(values = c("dashed","dashed","dashed"),
	labels=c("Control","Fungicide","Insecticide")) +
  geom_smooth(method='lm', 
              aes(x=lg_A.con, y=fit, colour=pesticide, fill=pesticide),alpha=0.4) +
  labs(x="", y="", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")+
  coord_cartesian(ylim=c(0.1,1.0))   


## ACPS
ACPS.surv.dat36 <- subset(sel.pred.surv.dat36, sp.=='ACPS')

p.ACPS.surv36 <- ggplot(data=ACPS.surv.dat36) +
  geom_smooth(method='lm', se = FALSE, 
              aes(x=lg_A.con, y=surv1, colour =pesticide, linetype=pesticide)) +
  scale_color_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  scale_linetype_manual(values = c("dashed","dashed","dashed"),
	labels=c("Control","Fungicide","Insecticide")) +
  geom_smooth(method='lm', 
              aes(x=lg_A.con, y=fit, colour=pesticide, fill=pesticide),alpha=0.4) +
  labs(x="", y="", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")+
  coord_cartesian(ylim=c(0.1,1.0))   


## ACTE
ACTE.surv.dat36 <- subset(sel.pred.surv.dat36, sp.=='ACTE')

p.ACTE.surv36 <- ggplot(data=ACTE.surv.dat36) +
  geom_smooth(method='lm', se = FALSE, 
              aes(x=lg_A.con, y=surv1, colour =pesticide, linetype=pesticide)) +
  scale_color_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  scale_linetype_manual(values = c("dashed","dashed","dashed"),
	labels=c("Control","Fungicide","Insecticide")) +
  geom_smooth(method='lm', 
              aes(x=lg_A.con, y=fit, colour=pesticide, fill=pesticide),alpha=0.4) +
  labs(x="", y="", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")+
  coord_cartesian(ylim=c(0.1,1.0))   


## FRMA
FRMA.surv.dat36 <- subset(sel.pred.surv.dat36, sp.=='FRMA')

p.FRMA.surv36 <- ggplot(data=FRMA.surv.dat36) +
  geom_smooth(method='lm', se = FALSE, 
              aes(x=lg_A.con, y=surv1, colour =pesticide, linetype=pesticide)) +
  scale_color_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  scale_linetype_manual(values = c("dashed","dashed","dashed"),
	labels=c("Control","Fungicide","Insecticide")) +
  geom_smooth(method='lm', 
              aes(x=lg_A.con, y=fit, colour=pesticide, fill=pesticide),alpha=0.4) +
  labs(x="", y="", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")+
  coord_cartesian(ylim=c(0.1,1.0))   


## TIAM
TIAM.surv.dat36 <- subset(sel.pred.surv.dat36, sp.=='TIAM')

p.TIAM.surv36 <- ggplot(data=TIAM.surv.dat36) +
  geom_smooth(method='lm', se = FALSE, 
              aes(x=lg_A.con, y=surv1, colour =pesticide, linetype=pesticide)) +
  scale_color_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  scale_linetype_manual(values = c("dashed","dashed","dashed"),
	labels=c("Control","Fungicide","Insecticide")) +
  geom_smooth(method='lm', 
              aes(x=lg_A.con, y=fit, colour=pesticide, fill=pesticide),alpha=0.4) +
  labs(x="", y="", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")+
  coord_cartesian(ylim=c(0.1,1.0))   


## plot the results of individual species together
require(gridExtra)
require(grid)

p.ABNE.surv36.2 <- p.ABNE.surv36 +   theme(legend.position = "none")
p.ACBA.surv36.2 <- p.ACBA.surv36 +   theme(legend.position = "none")
p.ACPS.surv36.2 <- p.ACPS.surv36 +   theme(legend.position = "none")
p.ACTE.surv36.2 <- p.ACTE.surv36 +   theme(legend.position = "none")
p.FRMA.surv36.2 <- p.FRMA.surv36 +   theme(legend.position = "none")
p.TIAM.surv36.2 <- p.TIAM.surv36 +   theme(legend.background = element_rect(fill="lightgray"),
							legend.position = c(0.8,0.85))

# add labels
p.ABNE.surv36.2 <- arrangeGrob(p.ABNE.surv36.2, top=textGrob("(a) A. nephrolepis", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))
p.ACBA.surv36.2 <- arrangeGrob(p.ACBA.surv36.2, top=textGrob("(b) A. barbinerve", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))
p.ACPS.surv36.2 <- arrangeGrob(p.ACPS.surv36.2, top=textGrob("(c) A. pseudo-sieboldianum", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))
p.ACTE.surv36.2 <- arrangeGrob(p.ACTE.surv36.2, top=textGrob("(d) A. tegmentosum", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))
p.FRMA.surv36.2 <- arrangeGrob(p.FRMA.surv36.2, top=textGrob("(e) F. mandshurica", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))
p.TIAM.surv36.2 <- arrangeGrob(p.TIAM.surv36.2, top=textGrob("(f) T. amurensis", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))
												   
pdf("results/fig.cbs20180324/Fig.8-Survival-prediction(glmmTMB)-individualspecies.pdf",
	height=9, width=12, useDingbats=FALSE)													   

grid.arrange(p.ABNE.surv36.2, p.ACBA.surv36.2, p.ACPS.surv36.2, 
	       p.ACTE.surv36.2, p.FRMA.surv36.2, p.TIAM.surv36.2, 
                  ncol=3, nrow = 2, vp=viewport(width=0.9, height=0.7),
                  left = textGrob("Survival", rot = 90, vjust =0.2, hjust=0.1,gp=gpar(fontsize=14)),
                  bottom = textGrob("Conspecific adult density (m2)", vjust=0.2,gp=gpar(fontsize=14)),
                  layout_matrix = rbind(c(1,2,3), c(4,5,6)),
                  widths = c(2.5, 2.5, 2.5), heights = c(2.5, 2.5))
dev.off()






#####################################################

### Extract coefficients 

#####################################################
summary(modTMB.tree.surv36)

res.TMBmod <- function(mod){
	betas <- fixef(mod)$cond[-1]
	se <- sqrt(diag(vcov(mod)$cond)[-1])
	zval <- betas/se
	pval <- 2*pnorm(abs(zval), lower.tail=FALSE)
	lower <- betas + qnorm(.025)*se
	upper <- betas + qnorm(.975)*se
	res.dat <- as.data.frame(cbind(betas, se, zval, pval, lower, upper))
}


## overall results
res.mod.tree.surv36 <- res.TMBmod(modTMB.tree.surv36)

## do not present the census and site
res.mod.tree.surv36 <- res.mod.tree.surv36[-c(6:9),]

## format the data frame
res.surv.dat <- function(data){
  names(data) <- c("betas", "se", "zval", "pval", "ci.lb", "ci.ub")
  rownames(data) <- c("Height", "Fungicide", "Insecticide","A.con", "A.het",
	"Fungicide:A.con","Insecticide:A.con")

  # re-order the rows, set the A.con at the top
  data <- data[c(4,2,3,6,7,5,1),]
  return(data)
}

res.mod.tree.surv.dat36.1 <- res.surv.dat(res.mod.tree.surv36)



##-------------------------------
## function to plot significant and insignificant effect
## the function to plot coefficients and the associated 95% CI
## The black circles indicate significant effects (P < 0.05), 
## grey circles signify marginally significant effects (0.05 < P < 0.1) 
## and white circles mean no significance.

PlotCoefCi1 <- function(model, yval){
  
  for(i in 1:length(model$betas)) {
  # error bars
  segments(x0=model$ci.lb[i], x1=model$ci.ub[i], y0=yval[i], y1=yval[i])
  
  # fill points according to p-values 
  if (model$pval[i] < 0.05) {
    points(y=yval[i], x=model$betas[i], pch=21, cex = 2, bg = "black")
  } else if (model$pval[i] >= 0.05 &  model$pval[i] < 0.1) {
    points(y=yval[i], x=model$betas[i], pch=21, cex = 2, bg = "gray")
  } else {
    points(y=yval[i], x=model$betas[i], pch=21, cex = 2, bg = "white") 
  }
  }
} # End of function PlotCoefInt()



##----------------- individual species ---------------------##\

res.mod.sp.surv36 <- lapply(surv_TMBmods36, function(x) res.TMBmod(x))
str(res.mod.sp.surv36)

res.mod.sp.surv36 <- do.call(rbind.data.frame, res.mod.sp.surv36)

## renames the rownames
res.mod.sp.surv36$parameter <- rep(c("Height", "Fungicide", "Insecticide",
	"A.con", "A.het","Fungicide:A.con","Insecticide:A.con"), times=6)
res.mod.sp.surv36$sp. <- rep(c("A. nephrolepis","A. barbinerve ","A. pseudo-sieboldianum",
				"A. tegmentosum","F. mandshurica","T. amurensis"), each=7)
names(res.mod.sp.surv36) <- c("betas", "se", "zval", "pval", "ci.lb", "ci.ub", "parameter", "sp.")



# get the idea of pval
res.mod.sp.surv36$tf <- ifelse(res.mod.sp.surv36$pval <.05, TRUE, FALSE)

# keep the original order appears in the data
res.mod.sp.surv36.1 <- res.mod.sp.surv36
res.mod.sp.surv36.1 <- transform(res.mod.sp.surv36.1, parameter=factor(parameter, levels=c("A.con", "Fungicide", "Insecticide",
	"Fungicide:A.con","Insecticide:A.con","A.het","Height")))

# re-order the parameters, set A.con at the top
res.mod.sp.surv36 <- res.mod.sp.surv36[c(4,2,3,6,7,5,1),]


res.mod.sp.surv36$sp_f = factor(res.mod.sp.surv36$sp., 
					levels=c("A. nephrolepis","A. barbinerve ","A. pseudo-sieboldianum",
					"A. tegmentosum","F. mandshurica","T. amurensis"))
str(res.mod.sp.surv36)

p.mod.sp.surv36 <- ggplot(data=res.mod.sp.surv36, 
	aes(y=parameter, x=betas, xmin=ci.lb, xmax=ci.ub)) + 
	labs(y="", x="Estimated Coefficients") + 
	geom_errorbarh(height = 0,size=0.3) + 
	geom_point(aes(shape=sig,fill=sig), size=2) + 
	scale_shape_manual(values=c(21,21,21))+
	scale_fill_manual(values = c("black", "gray", "white"))+
	geom_vline(xintercept = 0, linetype=2, color="darkgrey")+
	facet_wrap(~sp_f) +
	theme(legend.position = "none")  

ggsave(plot=p.mod.sp.surv36, 
	filename="results/fig.cbs20180324/Fig.S6-Survival-coefficients (glmmTMB)-individual-species.pdf", 
	dpi=300,width=10,height=6,useDingbats=FALSE)






############################
# Conditional R square
############################
## Citation: Nakagawa, S. & Schielzeth, H. (2013) A general and simple method for obtaining
## R2 from generalized linear mixed-effects models. Methods in Ecology and Evolution, 4, 133¨C142.

# VarF gets the variance in fitted values
# mF@X returns fixed effect design matrix

summary(mod.tree.surv27.4, correlation=FALSE)

# VarF <- var(as.vector(fixef(mod.tree.surv27.4) %*% t(mod.tree.surv27.4@X)))

library(arm)

cond.R2 <- function(){
(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1])/(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + VarCorr(mF)$Unit[1] + log(1 + 1/exp(as.numeric(fixef(m0)))))
}

## Failed.
## It is hard to implement the function to glmer model


## change to a new appraoch
## https://ecologyforacrowdedplanet.wordpress.com/2013/08/27/r-squared-in-mixed-models-the-easy-way/
require(MuMIn)
R2.modTMB.tree.surv30.3 <- r.squaredGLMM(modTMB.tree.surv30.3)

# FAILED



#####################

## Most common and less commom species

#####################

# Separate seedlings to two groups: abundant and less abundant species
# I consider TIAM and FRMA as most common species,
# the rest of species as less common species

pest.moco.surv.dat41 <- subset(pest.new.surv.dat40, sp.=="TIAM"|sp.=="FRMA")

mod.moco.surv41 <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het + 
			lg_S.con + lg_S.het + (1|quad.unique), 
				data=pest.moco.surv.dat41, family=binomial)
summary(mod.moco.surv41 , correlation=FALSE)



## less common
pest.leco.surv.dat41 <- subset(pest.new.surv.dat40, sp.!="TIAM"|sp.!="FRMA")

mod.leco.surv41 <- glmer(surv1 ~ lg_height0 + pesticide*lg_A.con + lg_A.het +
			 lg_S.con + lg_S.het + (1|quad.unique) + (1|sp.),
			 data=pest.leco.surv.dat41, family=binomial)

summary(mod.leco.surv41, correlation=FALSE)



##--------------------##
## prediction

p.moco.surv.dat41 <- pest.moco.surv.dat41


moco.surv.newdat <-  with(pest.moco.surv.dat41, 
		expand.grid(lg_A.con=seq(min(lg_A.con), max(lg_A.con), length=100),
		pesticide=c("W", "F", "I"), 
		lg_height0=mean(lg_height0), lg_A.het=mean(lg_A.het),
		lg_S.con=mean(lg_S.con), lg_S.het=mean(lg_S.het)))


leco.surv.newdat <-  with(pest.leco.surv.dat41, 
		expand.grid(lg_A.con=seq(min(lg_A.con), max(lg_A.con), length=100),
		pesticide=c("W", "F", "I"), 
		lg_height0=mean(lg_height0), lg_A.het=mean(lg_A.het),
		lg_S.con=mean(lg_S.con), lg_S.het=mean(lg_S.het)))


## bootstrap
boot.moco.surv <- bootMer(mod.moco.surv41, FUN=function(x)
		predict(x, moco.surv.newdat, re.form=~0, type='response'), nsim=99)

boot.leco.surv <- bootMer(mod.leco.surv41, FUN=function(x)
		predict(x, leco.surv.newdat, re.form=~0, type='response'), nsim=99)


surv.newdat41$lci.pred <- apply(boot.surv41$t, 2, quantile, probs = 0.025) 
surv.newdat41$uci.pred <- apply(boot.surv41$t, 2, quantile, probs = 0.975)
# surv.newdat41$mean.pred <- apply(boot.surv41$t, 2, quantile, probs = 0.5)
surv.newdat41$mean.pred <- predict(mod.tree.surv41, surv.newdat41, re.form=~0, type='response')



# back to the original data of x
p.tree.surv.dat41$logx_s <- scale(log(p.tree.surv.dat41$A.con + 1))
mean_logx <- mean(log(p.tree.surv.dat41$A.con + 1))
sd_logx <- sd(log(p.tree.surv.dat41$A.con + 1))
surv.newdat41$x <- exp((surv.newdat41$lg_A.con*sd_logx) + mean_logx)-1
summary(surv.newdat41)




# plot

p.pest.moco.surv36 <- ggplot(data=pest.moco.surv.dat36) +
  geom_smooth(method='lm', se = FALSE, 
	aes(x=lg_A.con, y=surv1, colour =pesticide, linetype=pesticide)) +
  scale_color_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  scale_linetype_manual(values = c("dashed","dashed","dashed"),
	labels=c("Control","Fungicide","Insecticide")) +
  geom_smooth(method='lm', 
              aes(x=lg_A.con, y=fit, colour=pesticide, fill=pesticide),alpha=0.4) +
  labs(x="", y="", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")


### less abundant species
p.pest.leco.surv36 <- ggplot(data=pest.leco.surv.dat36) +
  geom_smooth(method='lm', se = FALSE, 
              aes(x=lg_A.con, y=surv1, colour =pesticide, linetype=pesticide)) +
  scale_color_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  scale_linetype_manual(values = c("dashed","dashed","dashed"),
	labels=c("Control","Fungicide","Insecticide")) +
  geom_smooth(method='lm', 
              aes(x=lg_A.con, y=fit, colour=pesticide, fill=pesticide),alpha=0.4) +
  labs(x="", y="", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")  


## plot results of abundant and less abundant species together
## set the yaix as the same for two graphs

p.pest.moco.surv36.2 <- p.pest.moco.surv36 + 
				coord_cartesian(ylim=c(0.4,1.0)) +
 				theme(legend.position = "none")
p.pest.leco.surv36.2 <- p.pest.leco.surv36 + 
				coord_cartesian(ylim=c(0.4,1.0)) +
				theme(legend.background = element_rect(fill="lightgray"),
				legend.position = c(0.8,0.8))


# add labels
p.pest.moco.surv36.2 <- arrangeGrob(p.pest.moco.surv36.2, top=textGrob("(a) Most common", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))
p.pest.leco.surv36.2 <- arrangeGrob(p.pest.leco.surv36.2, top=textGrob("(b) Less common", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))

pdf("results/fig.cbs20180324/Fig.6-Pest-survival-prediction-common.pdf",
	height=7, width=10,useDingbats=FALSE)													   

grid.arrange(p.pest.moco.surv36.2, p.pest.leco.surv36.2,  
	nrow=1, ncol =2,vp=viewport(width=0.9, height=0.7),
                  left = textGrob("Survival", rot = 90, vjust =0.2, hjust=0.1, 
                                  gp=gpar(fontsize=14)),
                  bottom = textGrob("Conspecific adult density (m2)", vjust=0.2, 
                                    gp=gpar(fontsize=14)),
                  layout_matrix = rbind(c(1,2)),
                  widths = c(2.5, 2.5), heights = 2.5)
dev.off()




##----------------##
## Coefficients

## most common
res.mod.tree.moco.surv36 <- res.TMBmod(modTMB.tree.moco.surv36)

## do not present the census and site
res.mod.tree.moco.surv36 <- res.mod.tree.moco.surv36[-c(6:9),]
res.mod.tree.moco.surv.dat36 <- res.surv.dat(res.mod.tree.moco.surv36)

## less common
res.mod.tree.leco.surv36 <- res.TMBmod(modTMB.tree.leco.surv36)

## do not present the census and site
res.mod.tree.leco.surv36 <- res.mod.tree.leco.surv36[-c(6:9),]
res.mod.tree.leco.surv.dat36 <- res.surv.dat(res.mod.tree.leco.surv36)



## Plot the coefficients
# Figure parameters
pdf('results/fig.cbs20180324/Fig.S4-coefficients-pest-survival-commonness (glmmTMB).pdf',
	width=15,height=6, useDingbats=FALSE)

par(mfrow=c(1,2),lwd = 2, las = 1, mgp = c(2.5, .5, 0), tcl = -0.3,bty='n', mar = c(4, 12, 3.5, 1))

# Most common
plot(x = 1:4, y = 1:4, type = 'n',yaxt = "n",bty="o",
     ylim = c(-.3, 12.3), xlim = c(-1, 1),
     xlab = "",ylab = "", 
     cex.lab = 1.2, cex.axis = 1.2)
# xlab
title(xlab = "Estimated coefficients", cex.lab=1.5)

# vertical line at 0
abline(v = 0, lty = 2, col = "darkgrey")

axis(side = 2, rev(c(0,2,4,6,8,10,12)),tck=-0.02, labels = FALSE)


## The betas are not centered within the 95% CI.
PlotCoefCi1(model=res.mod.tree.moco.surv.dat36,yval = rev(c(0,2,4,6,8,10,12)))

# Plot x-axis labels  
text(y=rev(c(0,2,4,6,8,10,12)), x=-1.14,
     labels=c("Height","Fungicide", "Insecticide","A.con", "A.het",
	"Fungicide:A.con","Insecticide:A.con"),
     adj=1,xpd = TRUE,las=0,cex=1.5)

# main axis name
mtext(text = "(a) Most common", side = 3, line = 1.7, adj=0, at=-2.0, cex = 1.5, las = 0)


# Less common
plot(x = 1:4, y = 1:4, type = 'n',yaxt = "n",bty="o",
     ylim = c(-.3, 12.3), xlim = c(-1, 1),
     xlab = "",ylab = "", 
     cex.lab = 1.2, cex.axis = 1.2)
# xlab
title(xlab = "Estimated coefficients", cex.lab=1.5)

# vertical line at 0
abline(v = 0, lty = 2, col = "darkgrey")

axis(side = 2, rev(c(0,2,4,6,8,10,12)),tck=-0.02, labels = FALSE)

PlotCoefCi1(model=res.mod.tree.leco.surv.dat36,yval = rev(c(0,2,4,6,8,10,12)))

# Plot y-axis labels  
text(y=rev(c(0,2,4,6,8,10,12)), x=-1.14,
     labels=c("Height","Fungicide", "Insecticide","A.con", "A.het",
	"Fungicide:A.con","Insecticide:A.con"),
     adj=1,xpd = TRUE,las=0,cex=1.5)

# main axis name
mtext(text = "(b) Less common", side = 3, line = 1.7, adj=0, at=-1.8, cex = 1.5, las = 0)

dev.off() # end of the Fig. S4




###########################

## plot overall regeneration and survival together

###########################

## 1. Coefficients
## read regeneration data 
res.mod.tree.rege.dat18 <- read.csv("results/res.mod.tree.rege.dat18.csv")

# plot
pdf('results/fig.cbs20180324/Fig.1-coefficients-pest-regeneration&survival-overall.pdf',
	width=15, height=6, useDingbats=FALSE)

par(mfrow=c(1,2),lwd = 2, las = 1, mgp = c(2.5, .5, 0), tcl = -0.3,bty='n', mar = c(4, 12, 3, 1))

# Regeneration
plot(x = 1:4, y = 1:4, type = 'n',yaxt = "n",bty="o",
     ylim = c(-.6, 14.6), xlim = c(-.5, 1),
     xlab = "",ylab = "", 
     cex.lab = 1.2, cex.axis = 1)
# xlab
title(xlab = "Estimated coefficients", cex.lab=1.5)

# vertical line at 0
abline(v = 0, lty = 2, col = "darkgrey")

# axis at bottom
axis(side = 2, rev(c(0,2,4,6,8,10,12,14)),tck=-0.02, labels = FALSE)


## The betas are not centered within the 95% CI.
PlotCoefCi1(model=res.mod.tree.rege.dat18,yval = rev(c(0,2,4,6,8,10,12,14)))

# Plot y-axis labels  
text(y=rev(c(0,2,4,6,8,10,12,14)), x=-.6,
     labels=c("Fungicide", "Insecticide","A.con", "A.het",
	"Fungicide:A.con","Insecticide:A.con","Fungicide:A.het","Insecticide:A.het"),
     adj=1,xpd = TRUE,las=0,cex=1.5)

# main axis name
mtext(text = "(a) ", side = 3, line = 1.5, adj=0, at=-1.0, cex = 1.5, las = 0)


# Survival
plot(x = 1:4, y = 1:4, type = 'n',yaxt = "n",bty="o",
     ylim = c(-.3, 12.3), xlim = c(-.5, 1.0),
     xlab = "",ylab = "", 
     cex.lab = 1.2, cex.axis = 1.2)
# xlab
title(xlab = "Estimated coefficients", cex.lab=1.5)

# vertical line at 0
abline(v = 0, lty = 2, col = "darkgrey")

# axis at bottom
axis(side = 2, rev(c(0,2,4,6,8,10,12)),tck=-0.02, labels = FALSE)


## The betas are not centered within the 95% CI.
PlotCoefCi1(model=res.mod.tree.surv.dat36,yval = rev(c(0,2,4,6,8,10,12)))

# Plot y-axis labels  
text(y=rev(c(0,2,4,6,8,10,12)), x=-.6,
     labels=c("Height","Fungicide", "Insecticide","A.con", "A.het",
	"Fungicide:A.con","Insecticide:A.con"),
     adj=1,xpd = TRUE,las=0,cex=1.5)

# main axis name
mtext(text = "(b)", side = 3, line = 1.5, adj=0, at=-1.0, cex = 1.5, las = 0)

dev.off() # end of the Fig. 1

write.csv(res.mod.tree.rege.dat18,"results/res.mod.tree.rege.dat18.csv")
write.csv(res.mod.tree.surv.dat36,"results/res.mod.tree.surv.dat36.csv")


### 2. Prediction
## Recruitment
p.tree.recr.dat18 <- read.csv("results/p.tree.recr.dat18.csv")
summary(p.tree.recr.dat18)
# keep the original order appears in the data
p.tree.recr.dat18$pesticide <- factor(p.tree.recr.dat18$pesticide, levels=c('W', 'F', 'I'))

# plot
p.recr18 <- ggplot(p.tree.recr.dat18) +
  geom_smooth(method='lm', se = FALSE, 
              aes(x=lg_A.con, y=recruits, colour =pesticide, linetype=pesticide)) +
  scale_color_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  scale_linetype_manual(values = c("dashed","dashed","dashed"),
	labels=c("Control","Fungicide","Insecticide")) +
  geom_smooth(method='lm', 
              aes(x=lg_A.con, y=fit, colour=pesticide, fill=pesticide),alpha=0.4) +
  labs(x="Conspecific adult density (m2)", y="Recruitment (seedlings m-2)", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")+
  coord_cartesian(ylim=c(0.2,6.0))



## survival

p.recr18.2 <- p.recr18 +  theme(legend.position = "none")
p.surv36.2 <- p.surv36 + theme(legend.background = element_rect(fill="lightgray"),
					legend.position=c(0.8, 0.8))

# add labels
p.recr18.2 <- arrangeGrob(p.recr18.2, top=textGrob("(a)", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))
p.surv36.2 <- arrangeGrob(p.surv36.2, top=textGrob("(b)", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))

pdf("results/fig.cbs20180324/Fig.2-Pest-prediction-regeneration&survival-overall.pdf",
	height=7, width=10,useDingbats=FALSE)													   

grid.arrange(p.recr18.2, p.surv36.2, nrow=1, ncol=2,
	vp=viewport(width=0.9, height=0.6), heights =1.5, widths = c(2.5, 2.5))
                  
dev.off()