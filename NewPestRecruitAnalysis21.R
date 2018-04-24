### "Recruitment Analysis of the CBS Pesticide Experiment"
## By: "Shihong Jia"
## "Apr 23, 2018"

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

### Analysis 2: Recruitment analysis
The presentation of the relationship between conspecific adult density and seedling survival under three pesticide treatments (i.e. fungicide, insecticide and control). Both prediction and the confidence intervals and the observed values will be shown in one graph.


### Load packages
```{r loadlibs, warning=FALSE, message=FALSE}
library(ggplot2)
library(cowplot)
library(lme4)
library(pbkrtest)
library(tidyr)
library(dplyr)
library(optimx)
library(broom)
library(glmmTMB)
```

#### Data description
## Three years of spring censuses data (15-17) were used.
## Information of neigboring trees (area of conspecific, heterospecific and overall species) 20 radius around the focal quadrats are added to the dataset. 

## poisson distribution
## Species-level analysis: only those 
## Former version: NewPestRecruitmentAnalysis9
## Individual species only with conspecific neighboring adults occuring
## Species with more than 80 individuals
## Model: census or census + site
## abundant and less abundant groups[new]
## NewPestRecruitAnalysis16: prediction with simulate function
## NewPestRecruitAnalysis17: glmmTMB model and present with the original data (without summarize)
## All seedlings
## Redo the simulated regression figures (back transform, replace the observed line with points)



### Analysis of Tree Recruitment
pest.recr.dat13 <- read.csv("pest.recr.dat.csv")

pest.recr.dat13$pesticide <- factor(pest.recr.dat13$pesticide, levels=c('W', 'F', 'I'))
pest.recr.dat13$census <- factor(pest.recr.dat13$census, levels=c('16sp', '17sp'))
pest.recr.dat13$exclosure <- factor(pest.recr.dat13$exclosure)
pest.recr.dat13$site <- factor(pest.recr.dat13$site)
pest.recr.dat13 <-  as.data.frame(pest.recr.dat13)
summary(pest.recr.dat13)
write.csv(pest.recr.dat13,"results/pest.recr.dat.csv")

### initial plots
# between seasons
group_by(pest.recr.dat13, pesticide, exclosure, census) %>% 
  summarise(recrs = mean(recruits), 
            se.recrs=sd(recruits)/sqrt(length(recruits))) %>% 
  ggplot(aes(x=pesticide, y=recrs, ymin=recrs-se.recrs, ymax=recrs+se.recrs, colour=exclosure)) + 
  labs(x="Pesticide", y="Recruitment") + 
  geom_errorbar(position=position_dodge(width = .5), width = 0.5) + 
  geom_point(position=position_dodge(width = .5), size=4) +
  facet_grid(~census)



## A.con > 0?
hist(pest.recr.dat13$lg_A.con)
hist(subset(pest.recr.dat13,A.con >0)$lg_A.con)
dim(pest.recr.dat13); dim(subset(pest.recr.dat13,A.con >0))


## inteaction between fungicide and adult trees
mod.tree.recr15.1 <- glmer(recruits ~ pesticide*(lg_A.con + lg_A.het) + 
                                             (1|plot/quad.unique) + (1|sp.), 
                                           data=subset(pest.recr.dat13, A.con>0), family=poisson)
summary(mod.tree.recr15.1, correlation=FALSE)

# add census as fixed-effect and ramdom slope effect
mod.tree.recr15.2 <- glmer(recruits ~ pesticide*(lg_A.con + lg_A.het) + census +
                                             (1+census|plot/quad.unique) + (1|sp.), 
                                           data=subset(pest.recr.dat13, A.con>0), family=poisson)
summary(mod.tree.recr15.2, correlation=FALSE)

mod.tree.recr15.3 <- glmer(recruits ~ pesticide*(lg_A.con + lg_A.het) + census +
                                             (1+census|quad.unique) + (1|sp.), 
                                           data=subset(pest.recr.dat13, A.con>0), family=poisson)
summary(mod.tree.recr15.3, correlation=FALSE)

# census and site as fixed-effect
mod.tree.recr15.4 <- glmer(recruits ~ pesticide*(lg_A.con + lg_A.het) + census + site  +
                                             (1+census|quad.unique) + (1|sp.), 
                                           data=subset(pest.recr.dat13, A.con>0), family=poisson)
summary(mod.tree.recr15.4, correlation=FALSE)

anova(mod.tree.recr15.1, mod.tree.recr15.2, mod.tree.recr15.3, mod.tree.recr15.4)

# get results of the best model
summary(mod.tree.recr15.4, correlation=FALSE)

## plot the coefficients
require(coefplot2)
coefplot2(mod.tree.recr15.4, top.axis=F)


## glmmTMB model
modTMB.tree.recr18 <- glmmTMB(recruits ~ pesticide*(lg_A.con + lg_A.het) + census + site  +
                                             (1+census|quad.unique) + (1|sp.), 
                                           data=pest.recr.dat13, family=poisson)
summary(modTMB.tree.recr18, correlation=FALSE)


# include exclosure
modTMB.tree.recr18.5 <- glmmTMB(recruits ~ pesticide*(lg_A.con + lg_A.het) + exclosure + census + site  +
                                             (1+census|quad.unique) + (1|sp.), 
                                           data=pest.recr.dat13, family=poisson)
summary(modTMB.tree.recr18.5, correlation=FALSE)



###-------------------- predictions -------------------------###
### remove some outliers
range(p.tree.recr.dat18$A.con)
hist(p.tree.recr.dat18$A.con)
pest.recr.dat20 <- subset(pest.recr.dat13, A.con < 3000)
dim(pest.recr.dat13); dim(pest.recr.dat20)

# not consider the census and site as fixed-effect and remove census as a random slope in the model
mod.tree.recr20 <- glmer(recruits ~ pesticide*(lg_A.con + lg_A.het) + (1|quad.unique) + (1|sp.),
			 data=pest.recr.dat20, family=poisson(link=log))
summary(mod.tree.recr20, correlation=FALSE)


# link=log?
mod.tree.recr20.1 <- glmer(recruits ~ pesticide*(lg_A.con + lg_A.het) + (1|quad.unique) + (1|sp.),
			 data=pest.recr.dat20, family=poisson)
summary(mod.tree.recr20.1, correlation=FALSE)


# glmmTMB model
modTMB.tree.recr20 <- glmmTMB(recruits ~ pesticide*(lg_A.con + lg_A.het) + (1|quad.unique) + (1|sp.),
			 data=pest.recr.dat20, family=poisson(link=log))
summary(modTMB.tree.recr20, correlation=FALSE)

# same result




## prediction
p.tree.recr.dat20 <- pest.recr.dat20

recr.newdat20 <-  with(p.tree.recr.dat20, 
		expand.grid(lg_A.con=seq(min(lg_A.con), max(lg_A.con), length=100),
		pesticide=c("W", "F", "I"), lg_A.het=mean(lg_A.het)))



## 95% CI of the predicted values
## bootstrap
boot.recr20 <- bootMer(mod.tree.recr20, FUN=function(x)
		predict(x, recr.newdat20, re.form=~0, type='response'), nsim=99)


surv.newdat41$lci.pred <- apply(boot.surv41$t, 2, quantile, probs = 0.025) 
surv.newdat41$uci.pred <- apply(boot.surv41$t, 2, quantile, probs = 0.975)
# surv.newdat41$mean.pred <- apply(boot.surv41$t, 2, quantile, probs = 0.5)
recr.newdat20$mean.pred <- predict(mod.tree.recr20, recr.newdat20, re.form=~0, 
					type='response',allow.new.levels=TRUE)
recr.newdat20$mean.pred <- exp(recr.newdat20$mean.pred)


# back to the original data of x
p.tree.recr.dat20$logx_s <- scale(log(p.tree.recr.dat20$A.con + 1))
mean_logx <- mean(log(p.tree.recr.dat20$A.con + 1))
sd_logx <- sd(log(p.tree.recr.dat20$A.con + 1))
recr.newdat20$x <- exp((recr.newdat20$lg_A.con*sd_logx) + mean_logx)-1
summary(recr.newdat20)


ggplot(recr.newdat20) +
  geom_line(aes(x=x, y=mean.pred, colour=pesticide),size=0.8) +
#  geom_ribbon(aes(x=x, ymin = lci.pred, ymax = uci.pred, fill = pesticide), alpha = 0.3) +
  scale_color_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  labs(x="Conspecific adult density (m2)", y="Recruitment", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")


ggplot(recr.newdat20) +
  geom_line(aes(x=lg_A.con, y=mean.pred, colour=pesticide),size=0.8) +
#  geom_ribbon(aes(x=lg_A.con, ymin =lci.pred, ymax = uci.pred, fill = pesticide), alpha = 0.3) +
  scale_color_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  labs(x="Conspecific adult density (m2)", y="Recruitment", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")






## why the predicted values are so small
## another approach

# temporary data 
recr.tmpdat <- p.tree.recr.dat20[, c("pesticide", "lg_A.con", "lg_A.het",
		 "census", "site", "quad.unique", "sp.")]
dim(recr.tmpdat)


recr.jvalues <- with(p.tree.recr.dat20, seq(from = min(lg_A.con), to = max(lg_A.con), length.out = 100))


# calculate predicted probabilities and store in a list 
recr.biprobs <- lapply(levels(p.tree.recr.dat20$pesticide), function(treat) 
	{recr.tmpdat$pesticide[] <- treat 
	lapply(recr.jvalues, function(j) {
	recr.tmpdat$lg_A.con <- j 
	predict(mod.tree.recr20, newdata = recr.tmpdat, type = "response") 
  }) 
})


# get means and quartiles for all jvalues for each level of Pesticide 
recr.plotdat2 <- lapply(recr.biprobs, function(X) {
	 recr.temp <- t(sapply(X, function(x) {
	 c(M=mean(x), quantile(x, c(.25, .75))) 
	}))
	recr.temp <- as.data.frame(cbind(recr.temp, recr.jvalues))
	colnames(recr.temp) <- c("Predicted", "Lower", "Upper", "Lg_A.con_s") 
	return(recr.temp)
})

# collapse to one data frame 
recr.plotdat2 <- do.call(rbind, recr.plotdat2) 

# add pesticide treatment
recr.plotdat2$pesticide <- factor(rep(levels(p.tree.recr.dat20$pesticide), 
					each = length(recr.jvalues)))

recr.plotdat2$x <- exp((recr.plotdat2$Lg_A.con_s*sd_logx) + mean_logx)-1
recr.plotdat2$pesticide <- factor(recr.plotdat2$pesticide, levels=c("W","F","I"))
summary(recr.plotdat2)

ggplot(recr.plotdat2, aes(x = x, y = Predicted, color=pesticide)) + 
	geom_line(size=1) +
	geom_ribbon(aes(x=x, ymin = Lower, ymax = Upper, fill = pesticide), alpha = 0.2) +
	scale_color_manual(values = c('gray40','orangered','darkblue'),
		labels=c("Control","Fungicide","Insecticide")) + 
  	scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  		labs(x="Conspecific adult density (m2)", y="Recruitment", 
  	colour="Pesticide",linetype="Pesticide", fill="Pesticide")






#### add the observed values
## generate k sets of simulated y at these new points.
ob.recr.newdat <- with(p.tree.recr.dat20, 
		expand.grid(lg_A.con=lg_A.con, pesticide=levels(pesticide), 
		lg_A.het=mean(lg_A.het)))

dim(ob.recr.newdat)


ob.recr.newdat$pred <- predict(mod.tree.recr20, newdata=ob.recr.newdat, re.form=~0, type='link')
p.tree.recr.dat20$resd <- resid(mod.tree.recr20)


# select the lg_A.con in the original data
ob.recr.newdat$resd <- p.tree.recr.dat20[p.tree.recr.dat20$lg_A.con %in% ob.recr.newdat$lg_A.con,]$resd


ob.recr.newdat$obse_lk <- ob.recr.newdat$pred + ob.recr.newdat$resd
ob.recr.newdat$obse <- exp(ob.recr.newdat$obse_lk)
summary(ob.recr.newdat)




# replace the pesticide
ob.recr.newdat$pesticide <- p.tree.recr.dat20[p.tree.recr.dat20$lg_A.con %in% ob.recr.newdat$lg_A.con,]$pesticide


# back to the original data of x
ob.recr.newdat$x <- exp((ob.recr.newdat$lg_A.con*sd_logx) + mean_logx)-1
summary(ob.recr.newdat); dim(ob.recr.newdat)



## bin the x to 300 density
ob.recr.newdat$A.con_cat <- cut(ob.recr.newdat$x, seq(-1, 2999, 300))

obse.recr_cat <- summarise(group_by(ob.recr.newdat, A.con_cat, pesticide), 
			mean = mean(obse), se=sd(obse)/sqrt(length(obse)))
summary(obse.recr_cat); length(obse.recr_cat$A.con_cat)


# change the A.con_cat to continuous data
obse.recr_cat$x <- rep(seq(150,2850,300), each=3)


ggplot(recr.newdat20) +
  geom_line(aes(x=x, y=mean.pred, colour=pesticide),size=0.8) +
#  geom_ribbon(aes(x=x, ymin = lci.pred, ymax = uci.pred, fill = pesticide), alpha = 0.3) +
  scale_color_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  labs(x="Conspecific adult density (m2)", y="Recruitment", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")


ggplot() +
  geom_line(data=recr.newdat20, aes(x=x, y=mean.pred, colour=pesticide),size=0.8) +
#  geom_ribbon(recr.newdat20, aes(x=x, ymin = lci.pred, ymax = uci.pred, fill = pesticide), alpha = 0.2) +
  geom_point(data=obse.recr_cat,aes(x=x, y = mean, color=pesticide),
		size=2,position=position_dodge(width = 115)) +
  geom_errorbar(data=obse.recr_cat,aes(x=x, ymin=mean-1.96*se, ymax=mean+1.96*se,color=pesticide),
		width=0.1, position=position_dodge(width = 115)) +
  scale_color_manual(values = c('gray40','orangered','limegreen'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','limegreen'),
	labels=c("Control","Fungicide","Insecticide")) +
  labs(x="Conspecific adult density (m2)", y="Survival", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide") +
  theme(legend.title=element_blank(), legend.background = element_rect(fill="lightgray"),
		legend.position ="right") + theme_bw()


ggsave(plot=p.ob.pr.recr, 
	filename="results/fig.cbs20180418/Fig.1-recruitment-prediction-observation(glmmTMB)-overall.pdf", 
	dpi=300,width=10,height=6,useDingbats=FALSE)




## include random effects
## observed means and 1se
ggplot() +
  geom_line(data=recr.plotdat2, aes(x=x, y=Predicted, colour=pesticide),size=0.8) +
  geom_ribbon(data=recr.plotdat2, aes(x=x, ymin = Lower, ymax = Upper, fill = pesticide), alpha = 0.2) +
  geom_point(data=obse.recr_cat,aes(x=x, y = mean, color=pesticide),
		size=2,position=position_dodge(width = 115)) +
  geom_errorbar(data=obse.recr_cat,aes(x=x, ymin=mean-se, ymax=mean+se,color=pesticide),
		width=0.1, position=position_dodge(width = 115)) +
  scale_color_manual(values = c('gray40','orangered','limegreen'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','limegreen'),
	labels=c("Control","Fungicide","Insecticide")) +
  labs(x="Conspecific adult density (m2)", y="Recruitment", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide") +
  theme(legend.title=element_blank(), legend.background = element_rect(fill="lightgray"),
		legend.position ="right") + theme_bw()





###---------------- individual species analysis ----------------------------###
# For individual species, to fit the model better, we exclude those without conspecific neighboring adults
pest.recr.dat13.1 <- subset(pest.recr.dat13, A.con >0)
dim(pest.recr.dat13); dim(pest.recr.dat13.1)
# Exclude 89 data in total


## how many individuals in each species
ind.sp.pest.recr <- table(pest.recr.dat13.1$sp.)
write.csv(ind.sp.pest.recr, "results/recruitment-individual in each species.csv")


## criteria1 
sp.recr.counts <- summarise(group_by(pest.recr.dat13.1, sp.), 
                            n.quadrat=length(recruits), 
	max.recr=max(recruits), min.recr=min(recruits))

sp.recr.counts$range <- sp.recr.counts$max.recr - sp.recr.counts$min.recr
sp.recr.sel <- sp.recr.counts$sp.[sp.recr.counts$n.quadrat > 10 & sp.recr.counts$range > 5]
sp.recr.counts1 <- sp.recr.counts[sp.recr.counts$sp. %in% sp.recr.sel,]

write.csv(sp.recr.counts1, "results/recruitment-species-counts1(by recruit range).csv")

## criteria1
sp.recr.counts2 <- as.data.frame(table(pest.recr.dat13.1$sp., pest.recr.dat13.1$pesticide))
colnames(sp.recr.counts2) <- c('sp.', 'pesticide', 'no.quad')
sp.recr.counts2 <- sp.recr.counts2 %>% spread(pesticide, no.quad) 

write.csv(sp.recr.counts2, "results/recruitment-species-counts.csv")

sel.sp.pest.recr <- unique(sp.recr.counts2$sp.[sp.recr.counts2$W >= 10 & 
                                                          (sp.recr.counts2$F >= 10 |
                                                             sp.recr.counts2$I >= 10)])


# sp.recr.sel2 <- c('ABNE', 'ACBA',  'ACPS', 'ACTE', 'FRMA', 'PIKO', 'QUMO', 'TIAM')
# sp.recr.sel12 <- c('ACBA',  'ACPS', 'FRMA','TIAM')
sp.recr.sel2 <- c('ABNE', 'ACBA', 'ACPS', 'ACTE', 'FRMA', 'TIAM')


## glmmTMB model

recr_TMBmods18 <- sapply(sp.recr.sel2, function(sp){
  print(sp)
  spdat <- filter(pest.recr.dat13, sp.==sp)
  spdat <- droplevels(spdat)
  mod <- glmmTMB(recruits ~ pesticide*(lg_A.con + lg_A.het) + (1|quad.unique), 
               data=spdat, family=poisson)
}, simplify=FALSE)

lapply(recr_TMBmods18, function(x) summary(x, correlation=FALSE))


##################################
# predict with the glmmTMB model results

## prediction
recr.TMBmods18 <- sapply(sp.recr.sel2, function(sp){
  print(sp)
  spdat <- filter(pest.recr.dat13,sp.==sp)
  spdat <- droplevels(spdat)
  
  # model
  mod <- glmmTMB(recruits ~ pesticide*(lg_A.con + lg_A.het) + (1|quad.unique), 
               data=spdat, family=poisson)

  p.recr.dat <- spdat
  
  # calculate the confidence intervals
  p.recr.dat <- data.frame(p.recr.dat, 
	predict(mod, interval='confidence', type='link', se.fit=TRUE))
  
  # results
  p.recr.dat
  
}, simplify=FALSE)


## plot species with more than 80 individuals
sel.pred.recr.dat18 <- do.call(rbind.data.frame, recr.TMBmods18)
str(sel.pred.recr.dat18)


## plot one species by one species
## ABNE
ABNE.recr.dat18 <- subset(sel.pred.recr.dat18, sp.=='ABNE')

p.ABNE.recr18 <- ggplot(data=ABNE.recr.dat18) +
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
  labs(x="", y="", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")+
  coord_cartesian(ylim=c(0.2,8.5))


## ACBA
ACBA.recr.dat18 <- subset(sel.pred.recr.dat18, sp.=='ACBA')

p.ACBA.recr18 <- ggplot(data=ACBA.recr.dat18) +
  geom_smooth(method='lm', se = FALSE, 
              aes(x=lg_A.con, y=recruits, colour =pesticide, linetype=pesticide)) +
  scale_color_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values =c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  scale_linetype_manual(values = c("dashed","dashed","dashed"),
	labels=c("Control","Fungicide","Insecticide")) +
  geom_smooth(method='lm', 
              aes(x=lg_A.con, y=fit, colour=pesticide, fill=pesticide),alpha=0.4) +
  labs(x="", y="", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")+
  coord_cartesian(ylim=c(0.2,8.5))


## ACPS
ACPS.recr.dat18 <- subset(sel.pred.recr.dat18, sp.=='ACPS')

p.ACPS.recr18 <- ggplot(data=ACPS.recr.dat18) +
  geom_smooth(method='lm', se = FALSE, 
              aes(x=lg_A.con, y=recruits, colour =pesticide, linetype=pesticide)) +
  scale_color_manual(values =c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) + 
  scale_fill_manual(values = c('gray40','orangered','darkblue'),
	labels=c("Control","Fungicide","Insecticide")) +
  scale_linetype_manual(values = c("dashed","dashed","dashed"),
	labels=c("Control","Fungicide","Insecticide")) +
  geom_smooth(method='lm', 
              aes(x=lg_A.con, y=fit, colour=pesticide, fill=pesticide),alpha=0.4) +
  labs(x="", y="", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")+
  coord_cartesian(ylim=c(0.2,8.5))


## ACTE
ACTE.recr.dat18 <- subset(sel.pred.recr.dat18, sp.=='ACTE')

p.ACTE.recr18 <- ggplot(data=ACTE.recr.dat18) +
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
  labs(x="", y="", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")+
  coord_cartesian(ylim=c(0.2,8.5))


## FRMA
FRMA.recr.dat18 <- subset(sel.pred.recr.dat18, sp.=='FRMA')

p.FRMA.recr18 <- ggplot(data=FRMA.recr.dat18) +
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
  labs(x="", y="", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")+
  coord_cartesian(ylim=c(0.2,8.5))


## TIAM
TIAM.recr.dat18 <- subset(sel.pred.recr.dat18, sp.=='TIAM')

p.TIAM.recr18 <- ggplot(data=TIAM.recr.dat18) +
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
  labs(x="", y="", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")+
  coord_cartesian(ylim=c(0.2,8.5))


## plot the results of individual species together

require(gridExtra)
require(grid)

p.ABNE.recr18.2 <- p.ABNE.recr18 +   theme(legend.position = "none")
p.ACBA.recr18.2 <- p.ACBA.recr18 +   theme(legend.position = "none")
p.ACPS.recr18.2 <- p.ACPS.recr18 +   theme(legend.background = element_rect(fill="lightgray"),
							legend.position = c(0.8,0.85))
p.ACTE.recr18.2 <- p.ACTE.recr18 +   theme(legend.position = "none")
p.FRMA.recr18.2 <- p.FRMA.recr18 +   theme(legend.position = "none")
p.TIAM.recr18.2 <- p.TIAM.recr18 +   theme(legend.position = "none")


# add labels
p.ABNE.recr18.2 <- arrangeGrob(p.ABNE.recr18.2, top=textGrob("(a) A. nephrolepis", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))
p.ACBA.recr18.2 <- arrangeGrob(p.ACBA.recr18.2, top=textGrob("(b) A. barbinerve", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))
p.ACPS.recr18.2 <- arrangeGrob(p.ACPS.recr18.2, top=textGrob("(c) A. pseudo-sieboldianum", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))
p.ACTE.recr18.2 <- arrangeGrob(p.ACTE.recr18.2, top=textGrob("(d) A. tegmentosum", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))
p.FRMA.recr18.2 <- arrangeGrob(p.FRMA.recr18.2, top=textGrob("(e) F. mandshurica", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))
p.TIAM.recr18.2 <- arrangeGrob(p.TIAM.recr18.2, top=textGrob("(f) T. amurensis", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))

pdf("results/fig.cbs20180324/Fig.7-Pest-regeneration-prediction-individualspecies.pdf", 
	height=9, width=12, useDingbats=FALSE)

grid.arrange(p.ABNE.recr18.2, p.ACBA.recr18.2, p.ACPS.recr18.2, p.ACTE.recr18.2, 
	p.FRMA.recr18.2, p.TIAM.recr18.2, 
             ncol=3, nrow = 2, vp=viewport(width=0.9, height=0.7),
             left = textGrob("Recruitment (seedlings m-2)", 
				     rot = 90, vjust =0.2, hjust=0.3, 
                             gp=gpar(fontsize=14)),
             bottom = textGrob("Conspecific adult density (m2)", vjust=0.2, 
                               gp=gpar(fontsize=14)),
             layout_matrix = rbind(c(1,2,3), c(4,5,6)),
             widths = c(2.5, 2.5, 2.5), heights = c(2.5, 2.5))
dev.off()



#####################

## Most and less common species

#####################

# We consider TIAM and FRMA as most common species,
# the rest of species as less common species
pest.moco.recr.dat18 <- subset(pest.recr.dat13, sp.=="TIAM"|sp.=="FRMA")

modTMB.pest.moco.recr18 <- glmmTMB(recruits ~ pesticide*(lg_A.con + lg_A.het) + census + site +
	(1+census|quad.unique),data=pest.moco.recr.dat18, family=poisson)
summary(modTMB.pest.moco.recr18, correlation=FALSE)

## less common
pest.leco.recr.dat18 <- subset(pest.recr.dat13, sp.!="TIAM"|sp.!="FRMA")

modTMB.pest.leco.recr18 <- glmmTMB(recruits ~ pesticide*(lg_A.con + lg_A.het) + census + site +
	(1+census|quad.unique)+(1|sp.),data=pest.leco.recr.dat18, family=poisson)
summary(modTMB.pest.leco.recr18, correlation=FALSE)


##--------------------##
## prediction

pest.moco.recr.dat18 <- data.frame(pest.moco.recr.dat18, 
	predict(modTMB.pest.moco.recr18, interval='confidence', type='link', se.fit=TRUE))

pest.leco.recr.dat18 <- data.frame(pest.leco.recr.dat18, 
	predict(modTMB.pest.leco.recr18, interval='confidence', type='link', se.fit=TRUE))

# plot
p.pest.moco.recr18 <- ggplot(data=pest.moco.recr.dat18) +
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
  labs(x="", y="", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")+
  coord_cartesian(ylim=c(0.2,8.0))


p.pest.leco.recr18 <- ggplot(data=pest.leco.recr.dat18) +
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
  labs(x="", y="", 
  colour="Pesticide",linetype="Pesticide", fill="Pesticide")+
  coord_cartesian(ylim=c(0.2,8.0))


## plot results of abundant and less abundant species together
p.pest.moco.recr18.2 <- p.pest.moco.recr18 +   theme(legend.position = "none")
p.pest.leco.recr18.2 <- p.pest.leco.recr18 +   theme(legend.background = element_rect(fill="lightgray"),
								     legend.position = c(0.8,0.85))

# add labels
p.pest.moco.recr18.2 <- arrangeGrob(p.pest.moco.recr18.2, top=textGrob("(a) Most common", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))
p.pest.leco.recr18.2 <- arrangeGrob(p.pest.leco.recr18.2, top=textGrob("(b) Less common", x = unit(0, "npc"), 
                                                       y=unit(1, "npc"), just=c("left","top")))

pdf("results/fig.cbs20180324/Fig.5-Pest-regeneration-prediction-commonness.pdf",
	height=7, width=10,useDingbats=FALSE)													   

grid.arrange(p.pest.moco.recr18.2, p.pest.leco.recr18.2, 
	nrow=1, ncol =2,vp=viewport(width=0.9, height=0.7),
                  left = textGrob("Recruitment (seedlings m-2)", 
				     rot = 90, vjust =0.2, hjust=0.3,gp=gpar(fontsize=14)),
                  bottom = textGrob("Conspecific adult density (m2)", vjust=0.2,gp=gpar(fontsize=14)),
                  layout_matrix = rbind(c(1,2)),
                  widths = c(2.5, 2.5), heights = 2.5)
dev.off()





############################
# Conditional R square
############################
## old results from glmer

## change to a new appraoch
## https://ecologyforacrowdedplanet.wordpress.com/2013/08/27/r-squared-in-mixed-models-the-easy-way/
require(MuMIn)
R2.mod.tree.recr13.2 <- r.squaredGLMM(mod.tree.recr13.2)




############################
# Extract the model results
############################
# lower returns to the lower 95% CI of the coefficient
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
res.mod.tree.rege18 <- res.TMBmod(modTMB.tree.recr18)
res.mod.tree.rege18 <- res.mod.tree.rege18[-c(5:7),]

## format the data frame
res.rege.dat <- function(data){
  names(data) <- c("betas", "se", "zval", "pval", "ci.lb", "ci.ub")
  rownames(data) <- c("Fungicide", "Insecticide","A.con", "A.het",
	"Fungicide:A.con","Insecticide:A.con","Fungicide:A.het","Insecticide:A.het")
  return(data)
}

res.mod.tree.rege.dat18 <- res.rege.dat(res.mod.tree.rege18)
write.csv(res.mod.tree.rege.dat18, "results/res.mod.tree.rege.dat18.csv")




##----------------- individual species ---------------------##
res.mod.sp.recr18 <- lapply(recr_TMBmods18, function(x) res.TMBmod(x))
res.mod.sp.recr18 <- do.call(rbind.data.frame, res.mod.sp.recr18)
str(res.mod.sp.recr18)

## renames the rownames
res.mod.sp.recr18$Parameter <- rep(c("Fungicide", "Insecticide","A.con", "A.het",
	"Fungicide:A.con", "Insecticide:A.con", 
	"Fungicide:A.het","Insecticide:A.het"), times=6)
res.mod.sp.recr18$sp. <- rep(c("A. nephrolepis","A. barbinerve ","A. pseudo-sieboldianum",
				"A. tegmentosum","F. mandshurica","T. amurensis"), each=8)
names(res.mod.sp.recr18) <- c("betas", "se", "zval", "pval", "ci.lb", "ci.ub", "parameter", "sp.")
str(res.mod.sp.recr18)

# get the idea of pval
# res.mod.sp.recr18$tf <- ifelse(res.mod.sp.recr18$pval <.05, TRUE, FALSE)

res.mod.sp.recr18$sig <- cut(res.mod.sp.recr18$pval, breaks = c(0, 0.05, 0.1, 1), labels = c('Sig', 'Marg.sig', 'Non-sig'))

# reverse the order of origianl data
res.mod.sp.recr18 <- transform(res.mod.sp.recr18, 
				parameter=factor(parameter, levels=rev(unique(parameter))))

res.mod.sp.recr18$sp_f = factor(res.mod.sp.recr18$sp., 
					levels=c("A. nephrolepis","A. barbinerve ","A. pseudo-sieboldianum",
					"A. tegmentosum","F. mandshurica","T. amurensis"))

p.mod.sp.recr18 <- ggplot(data=res.mod.sp.recr18, 
	aes(y=parameter, x=betas, xmin=ci.lb, xmax=ci.ub)) + 
	labs(y="", x="Estimated Coefficients") +
	geom_errorbarh(height = 0,size=0.3) +  
	geom_point(aes(shape=sig,fill=sig),size=2) + 
	scale_shape_manual(values=c(21,21,21))+
	scale_fill_manual(values = c("black", "gray", "white"))+
	geom_vline(xintercept = 0, linetype=2, color="darkgrey")+ 
	facet_wrap(~sp_f) +
	theme(legend.position = "none")+
  	coord_cartesian(xlim=c(-8,8))


ggsave(plot=p.mod.sp.recr18, 
	filename="results/fig.cbs20180324/Fig.S5-Pest-regeneration-coefficients-individualspecies.pdf", 
	dpi=300,width=10,height=6, useDingbats=FALSE)



#####################

## Most and less common species

#####################

## Coefficients

## most common
res.mod.pest.moco.recr18 <- res.TMBmod(modTMB.pest.moco.recr18)

## do not present the census and site
res.mod.pest.moco.recr18 <- res.mod.pest.moco.recr18[-c(5:7),]
res.mod.pest.moco.recr.dat18 <- res.rege.dat(res.mod.pest.moco.recr18)

## less common
res.mod.pest.leco.recr18 <- res.TMBmod(modTMB.pest.leco.recr18)

## do not present the census and site
res.mod.pest.leco.recr18 <- res.mod.pest.leco.recr18[-c(5:7),]
res.mod.pest.leco.recr.dat18 <- res.rege.dat(res.mod.pest.leco.recr18)



## Plot the coefficients
# Figure parameters
pdf('results/fig.cbs20180324/Fig.S3-coefficients-pest-recruitment-commonness.pdf',
	width=15,height=6, useDingbats=FALSE)

par(mfrow=c(1,2),lwd = 2, las = 1, mgp = c(2.5, .5, 0), tcl = -0.3,bty='n', mar = c(4, 12, 3.5, 1))

# Abundant
plot(x = 1:4, y = 1:4, type = 'n',yaxt = "n",bty="o",
     xlim = c(-0.5, 1.0),ylim = c(-.3, 14.3), 
     xlab = "",ylab = "", 
     cex.lab = 1.2, cex.axis = 1.2)
# xlab
title(xlab = "Estimated coefficients", cex.lab=1.5)

# horizontal line at 0
abline(v = 0, lty = 2, col = "darkgrey")

# axis at bottom
axis(side = 2, c(0,2,4,6,8,10,12,14),tck=-0.02, labels = FALSE)


## The betas are not centered within the 95% CI.
PlotCoefCi1(model=res.mod.pest.moco.recr.dat18,yval = c(14,12,10,8,6,4,2,0))

# Plot x-axis labels  
text(y=c(14,12,10,8,6,4,2,0), x=-.635,
     labels=c("Fungicide", "Insecticide","A.con", "A.het",
	"Fungicide:A.con","Insecticide:A.con","Fungicide:A.het","Insecticide:A.het"),
     adj=1,xpd = TRUE,las=0,cex=1.5)

# main axis name
mtext(text = "(a) Most common", side = 3, line = 1.7, adj=0, at=-1.25, cex = 1.5, las = 0)


# Less abundant
plot(x = 1:4, y = 1:4, type = 'n',yaxt = "n",bty="o",
     xlim = c(-.5, 1.0),ylim = c(-.3, 14.3), 
     xlab = "",ylab = "", 
     cex.lab = 1.2, cex.axis = 1.2)
# xlab
title(xlab = "Estimated coefficients", cex.lab=1.5)

# horizontal line at 0
abline(v = 0, lty = 2, col = "darkgrey")

# axis at bottom
axis(side = 2, c(0,2,4,6,8,10,12,14),tck=-0.02, labels = FALSE)


## The betas are not centered within the 95% CI.
PlotCoefCi1(model=res.mod.pest.leco.recr.dat18,yval = c(14,12,10,8,6,4,2,0))

# Plot x-axis labels  
text(y=c(14,12,10,8,6,4,2,0), x=-.65,
     labels=c("Fungicide", "Insecticide","A.con", "A.het",
	"Fungicide:A.con","Insecticide:A.con","Fungicide:A.het","Insecticide:A.het"),
     adj=1,xpd = TRUE,las=0,cex=1.5)

# main axis name
mtext(text = "(b) Less common", side = 3, line = 1.7, adj=0, at=-1.25, cex = 1.5, las = 0)

dev.off() # end of the Fig. S3