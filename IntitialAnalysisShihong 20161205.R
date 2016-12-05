## Initial analysis

library(ggplot2)

#dat <- read.csv("../data/Woodyseedlings 20161205.csv")
dat <- read.csv("data/Woodyseedlings 20161205.csv")
summary(dat)

# remove the NAs
# dat15b <- height15b dataset
dat.15b <- subset(dat,dat$height15b!="NA")
summary(dat.15b)

# survival for the 16a
dat.15b$surv1 <-as.numeric(dat.15b$height16a !=0 & !is.na(dat.15b$height16a))
# survival for the 16b
dat.15b$surv2 <- as.numeric(dat.15b$height16b !=0 & !is.na(dat.15b$height16b))
quad <-  as.factor(dat$code):dat$sub.code


dat.ag <-  aggregate(cbind(surv1, surv2) ~ treatment + exclosure + species + quad, data=dat.15b, FUN=mean)
summary(dat.ag)
ggplot(data=dat.ag, aes(x=treatment, y=surv1, colour=as.factor(exclosure))) + geom_point(position='jitter')
ggplot(data=dat.ag, aes(x=treatment, y=surv1, colour=as.factor(exclosure))) + geom_violin()
ggplot(data=dat.ag, aes(x=treatment, y=surv2, colour=as.factor(exclosure))) + geom_violin()


## play with data
# quadrats code
quad <-  as.factor(dat.15b$code):dat.15b$sub.code

# aggregte the data for the survival
dat.ag <-  aggregate(cbind(surv1, surv2, site, code, sub.code) ~ treatment + exclosure + quad, data=dat.15b, FUN=mean)

summary(dat.ag1)
View(dat.ag1)

ggplot(data=dat.ag, aes(x=treatment, y=surv1, colour=as.factor(exclosure))) + geom_violin()
ggplot(data=dat.ag, aes(x=treatment, y=surv2, colour=as.factor(exclosure))) + geom_violin()


p <- ggplot(data=dat.ag, aes(x=treatment, y=surv1,colour=as.factor(exclosure)))
p <- p + geom_violin()
p + facet_wrap( "site" )


# mixed model text
library(Matrix)
library(lme4)

# simple model
mod <- lmer(surv1~treatment+(1|code),data=dat.15b)
summary(mod)

mod <- lmer(surv1~treatment+(1|code/sub.code),data=dat.15b)
summary(mod)

mod <- lmer(surv1~treatment+(1|code/sub.code)+(1|species),data=dat.15b)
summary(mod)

summary(mod)$coef
summary(mod)$varcor
ranef(mod)

fixef(mod)

## gives a plot of residuals vs fitted values
plot(mod)
## homoscedascity assumption
plot(mod, sqrt(abs(resid(.))) ~ fitted(.))
plot(mod, sqrt(abs(resid(.))) ~ fitted(.), type=c('p', 'smooth'))

## test the normality of the residuals
library(car)
qqPlot(resid(mod))
qqPlot(ranef(mod)$code[,1])


# bootstrapping
sims <- bootMer(mod3, FUN=fixef, nsim=100)
summary(sims)
dim(sims$t)
apply(sims$t, 2, quantile, c(0.025, 0.975))

anova(mod)

library(pbkrtest)
## first build a null model
mod0 <- update(mod, ~.-treatment)
## Then do a bootstrap comparison
sims.anova <- PBmodcomp(mod, mod0, nsim=100)
summary(sims.anova)



