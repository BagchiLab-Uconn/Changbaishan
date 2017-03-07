########################################################################
#########>>>>>>>>>>> Other vegetation metric analyses ##################
########################################################################

#### Abundance, shannon diversity, simpson diversity, evenness (not all) and richness
#### Analyzing pesticide and snow removal treatment seperately




####################
## no. of species ##
####################


## think about how to select with or without NA
## below is not correct
gathered.w.pest.spec <- subset(gathered.w.pest,gathered.w.pest$height !=0 & !is.na(gathered.w.pest$height))
gathered.w.pest.spec$no.spec <-as.numeric(gathered.w.pest.spec$species)
gathered.w.pest.spri  <- aggregate(no.spec ~ pesticide + exclosure + site + quadrat + census, 
                                   data=gathered.w.pest.spec, FUN=function(x)length(unique(x)))

# par(mfrow=c(2,2))
hist(gathered.w.pest.spri$no.spec)
#hist(log(gathered.w.pest.spri$no.spec))
#hist(log10(gathered.w.pest.spri$no.spec))
#hist(sqrt(gathered.w.pest.spri$no.spec))

###### Transform no.spec
# model1
m.w.pest.spri1 <- lmer(sqrt(no.spec) ~ pesticide*exclosure + census + (1|site),
                       data=gathered.w.pest.spri)
summary(m.w.pest.spri1)

plot(m.w.pest.spri1)

# model2
m.w.pest.spri2 <- lmer(log(no.spec) ~ pesticide*exclosure + census + (1|site),
                       data=gathered.w.pest.spri)
summary(m.w.pest.spri2)

plot(m.w.pest.spri2)

# give the result of sqrt-transformation
# m.w.pest.spri <- m.w.pest.spri1

####### Without transforming no.spec
m.w.pest.spri3 <- lmer(no.spec ~ pesticide*exclosure + census + (1|site),
                       data=gathered.w.pest.spri)
summary(m.w.pest.spri3)

plot(m.w.pest.spri3)

# model2
m.w.pest.spri4 <- lmer(no.spec ~ pesticide + exclosure + census + (1|site),
                       data=gathered.w.pest.spri)
summary(m.w.pest.spri4)
plot(m.w.pest.spri4)

anova(m.w.pest.spri3, m.w.pest.spri4)
m.w.pest.spri <- m.w.pest.spri3

## residuals vs. fittd values
test.w.pest.spri <- data.frame(res = resid(m.w.pest.spri), fit = fitted(m.w.pest.spri), 
                               trt = gathered.w.pest.spri$pesticide, excl = gathered.w.pest.spri$exclosure,
                               site = gathered.w.pest.spri$site, cens = gathered.w.pest.spri$census)

###
qplot(x=fit, y=res, geom='smooth', data=test.w.pest.spri) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("no. of species")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.w.pest.spri, facets=~trt) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("no. of species")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.w.pest.spri, colour=trt, facets=~site) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("no. of species in each site")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()

qplot(x=fit, y=res, geom='smooth', data=test.w.pest.spri, colour=trt, facets=~cens) + 
  xlab("fitted")+
  ylab("resi16a")+
  ggtitle("no. of species in each census")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point()


####################
### bootstrap
####################
new.w.pest.spri <- expand.grid(pesticide = c("C", "F","I","W"), exclosure = c("0","1"),
                               census= c("15fa","16sp","16fa"))

bootfit.w.pest.spri <- bootMer(m.w.pest.spri, FUN=function(x)predict(x, new.w.pest.spri, re.form=~0), nsim=999)

new.w.pest.spri$lci <- apply(bootfit.w.pest.spri$t, 2, quantile, 0.025) 
new.w.pest.spri$uci <- apply(bootfit.w.pest.spri$t, 2, quantile, 0.975) 
new.w.pest.spri$shanpred <- predict(m.w.pest.spri, newdata= new.w.pest.spri, re.form=~0)
new.w.pest.spri$shanpred1 <- predict(m.w.pest.spri, newdata= new.w.pest.spri, type = "link", re.form=~0)

#pred.pest.spri <- within(new.pest.spri, {
#  Predprob <- plogis(shanpred)
#  LL <- plogis(lci)
#  UL <- plogis(uci)
#})

#pred.pest.spri.dat <- cbind(new.pest.spri, pred.pest.spri)


### plot the original results
ggplot(new.w.pest.spri, aes(x = pesticide, y = shanpred, ymin = lci, ymax = uci,
                            color=as.factor(exclosure))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("census") + theme_bw() +
  ylab("No.of species")


ggplot(new.w.pest.spri, aes(x = census, y = shanpred, ymin = lci, ymax = uci,
                            color=as.factor(pesticide))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("exclosure") + theme_bw() +
  ylab("No.of species")



#################################################################
###------------------ Data: snow -----------------------------###
gathered.w.snow.spec <- subset(gathered.w.snow, gathered.w.snow$height !=0 & !is.na(gathered.w.snow$height))
gathered.w.snow.spec$no.spec <-as.numeric(gathered.w.snow.spec$species)
gathered.w.snow.spri  <- aggregate(no.spec ~ treatment + site + quadrat + census, 
                                 data=gathered.w.snow.spec, FUN=function(x)length(unique(x)))

#par(mfrow=c(2,2))
hist(gathered.w.snow.spri$no.spec)
#hist(log(gathered.snow.spri$no.spec))
#hist(sqrt(gathered.snow.spri$no.spec))

###### model
### Transforming 
m.w.snow.spri <- lmer(sqrt(no.spec) ~ treatment + census + (1|site),
                    data=gathered.w.snow.spri)
summary(m.w.snow.spri)
plot(m.w.snow.spri)

### Without Transforming 
m.w.snow.spri2 <- lmer(no.spec ~ treatment + census + (1|site),
                      data=gathered.w.snow.spri)
m.w.snow.spri3 <- lmer(no.spec ~ treatment*census + (1|site),
                       data=gathered.w.snow.spri)

anova(m.w.snow.spri1, m.w.snow.spri2, m.w.snow.spri3)

m.w.snow.spri <- m.w.snow.spri3

summary(m.w.snow.spri3)
plot(m.w.snow.spri3)


summary(m.w.snow.spri2)
plot(m.w.snow.spri2)


### bootstrap

new.w.snow.spri <- expand.grid(treatment = c("C","S"), census= c("15sp","15fa","16sp","16fa"))

bootfit.w.snow.spri <- bootMer(m.w.snow.spri, FUN=function(x)predict(x, new.w.snow.spri, re.form=~0), nsim=999)

new.w.snow.spri$lci <- apply(bootfit.w.snow.spri$t, 2, quantile, 0.025) 
new.w.snow.spri$uci <- apply(bootfit.w.snow.spri$t, 2, quantile, 0.975) 
new.w.snow.spri$spripred <- predict(m.w.snow.spri, newdata= new.w.snow.spri, re.form=~0)

### plot the original results

ggplot(new.w.snow.spri, aes(x = census, y = spripred, ymin = lci, ymax = uci,
                          color=as.factor(treatment))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  theme_bw() +
  ylab("No.of species(sqrt)")




############################################################################
################################ Abundance #################################
############################################################################

###-------------------------- Pest ----------------------------###
gathered.w.pest$abun <- as.numeric(gathered.w.pest$height !=0 & !is.na(gathered.w.pest$height))

pest.w.abun  <- aggregate(abun ~ pesticide + exclosure + site + quadrat + census, 
                        data=gathered.w.pest, FUN=sum)
par(mfrow = c(2,2))
hist(pest.w.abun$abun)
hist(log(pest.w.abun$abun))
hist(log10(pest.w.abun$abun))
hist(sqrt(pest.w.abun$abun))

###### model
m.w.pest.abun <- lmer(sqrt(abun) ~ pesticide*exclosure + census + (1|site), 
                      data=pest.w.abun)
summary(m.w.pest.abun)
plot(m.w.pest.abun)

### Without transforming abun
m.w.pest.abun2 <- lmer(abun ~ pesticide*exclosure + census + (1|site), 
                      data=pest.w.abun)
summary(m.w.pest.abun2)
plot(m.w.pest.abun2)

m.w.pest.abun2 <- lmer(abun ~ pesticide*exclosure + census + (1|site), 
                       data=pest.w.abun)
m.w.pest.abun3 <- lmer(abun ~ pesticide + exclosure + census + (1|site), 
                       data=pest.w.abun)
m.w.pest.abun4 <- lmer(abun ~ pesticide + census + (1|site), 
                       data=pest.w.abun)
anova(m.w.pest.abun2, m.w.pest.abun3, m.w.pest.abun4)
anova(m.w.pest.abun3, m.w.pest.abun4)

summary(m.w.pest.abun4)
plot(m.w.pest.abun4)



############# bootstrap ##################
new.w.pest.abun <- expand.grid(pesticide = c("C", "F","I","W"), exclosure = c("0","1"),
                             census= c("16sp","16fa"))

bootfit.w.pest.abun <- bootMer(m.w.pest.abun, FUN=function(x)predict(x, new.w.pest.abun, re.form=~0), nsim=999)

new.w.pest.abun$lci <- apply(bootfit.w.pest.abun$t, 2, quantile, 0.025) 
new.w.pest.abun$uci <- apply(bootfit.w.pest.abun$t, 2, quantile, 0.975) 
new.w.pest.abun$abunpred <- predict(m.w.pest.abun, newdata= new.w.pest.abun, re.form=~0)

### plot the original results
ggplot(new.w.pest.abun, aes(x = pesticide, y = abunpred, ymin = lci, ymax = uci,
                          color=as.factor(exclosure))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("census") + theme_bw() +
  ylab("Abundance(sqrt)")


ggplot(new.w.pest.abun, aes(x = census, y = abunpred, ymin = lci, ymax = uci,
                          color=as.factor(pesticide))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("exclosure") + theme_bw() +
  ylab("Abundance(sqrt)")




###-------------------------- Snow ----------------------------###
gathered.w.snow$abun <- as.numeric(gathered.w.snow$height !=0 & !is.na(gathered.w.snow$height))

w.snow.abun  <- aggregate(abun ~ treatment +  site + quadrat + census, 
                        data=gathered.w.snow, FUN=sum)

dev.new()
par(mfrow = c(2,2))
hist(w.snow.abun$abun)
hist(log(w.snow.abun$abun))
hist(log10(w.snow.abun$abun))
hist(sqrt(w.snow.abun$abun))

###### model
m.w.snow.abun <- lmer(sqrt(abun) ~ treatment + census + (1|site), data=w.snow.abun)
summary(m.w.snow.abun)
plot(m.w.snow.abun)

### Without transforming abun
m.w.snow.abun1 <- lmer(abun ~ treatment + census + (1|site), data=w.snow.abun)
summary(m.w.snow.abun1)
plot(m.w.snow.abun1)



############# bootstrap ##################
new.w.snow.abun <- expand.grid(treatment = c("C","S"), census= c("15sp","15fa","16sp","16fa"))

bootfit.w.snow.abun <- bootMer(m.w.snow.abun, FUN=function(x)predict(x, new.w.snow.abun, re.form=~0), nsim=999)

new.w.snow.abun$lci <- apply(bootfit.w.snow.abun$t, 2, quantile, 0.025) 
new.w.snow.abun$uci <- apply(bootfit.w.snow.abun$t, 2, quantile, 0.975) 
new.w.snow.abun$abunpred <- predict(m.w.snow.abun, newdata= new.w.snow.abun, re.form=~0)

### plot the original results

ggplot(new.snow.abun, aes(x = census, y = abunpred, ymin = lci, ymax = uci,
                          color=as.factor(treatment))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  theme_bw() +
  ylab("Abundance(sqrt)")







############################################################################
######################## community analysis ################################
############################################################################
### Abundance for each species in each quadrat in each census
gathered.w.pest$abun <- as.numeric(gathered.w.pest$height !=0 & !is.na(gathered.w.pest$height))
gathered.w.pest.spec.abun  <- aggregate(abun ~ pesticide + exclosure + quad1 + census + species, 
                                        data=gathered.w.pest, FUN=sum)
# View(gathered.pest.spec.abun)

unique(gathered.w.pest.spec.abun$species)

### Make the data available for community analysis
gathered.w.pest.spec.abun1 <- unite(gathered.w.pest.spec.abun, "new", pesticide:census, sep = "_")

# CommData1 <- as.data.frame.matrix(with(gathered.pest.spec.abun1,table(new1,species)))
# View(CommData1)

gathered.w.pest.spec.abun2 <- spread(gathered.w.pest.spec.abun1, species, abun)
#View(gathered.w.pest.spec.abun2)
gathered.w.pest.spec.abun2[is.na(gathered.w.pest.spec.abun2)] <- 0
gathered.w.pest.spec.abun3 <- separate(gathered.w.pest.spec.abun2,new,into=c("pesticide","exclosure","quad1","census"),"_")
gathered.w.pest.spec.abun4 <- separate(gathered.w.pest.spec.abun3,quad1,into=c("site","quadrat"),":")
gathered.w.pest.spec.abun5 <- unite(gathered.w.pest.spec.abun4, "quad1", site:quadrat, sep = "")


# gathered.w.pest.spec.abun6 <- gathered.w.pest.spec.abun5
# CommPest1 <- subset(CommPest,,-c("pesticide", "exclosure", "census"))
# CommPest1 <- CommPest[,-c("pesticide", "exclosure", "census")]
Comm.w.pest <- gathered.w.pest.spec.abun5[,-c(1, 2, 3, 4)]


###
library(vegan)
library(permute)
library(lattice)


## Shannon's index
w.pest.shan <- diversity(Comm.w.pest)

hist(w.pest.shan)

Comm.w.pest.dat <- data.frame(gathered.w.pest.spec.abun5$pesticide, 
                              gathered.w.pest.spec.abun5$exclosure, 
                              gathered.w.pest.spec.abun5$quad1, 
                              gathered.w.pest.spec.abun5$census,
                              w.pest.shan)
### rename columns
colnames(Comm.w.pest.dat) <- c("pesticide","exclosure","quad","census","shan")
### fit in the model
m.w.pest.shan1 <- lmer(shan ~ pesticide*exclosure + census + (1|quad), data=Comm.w.pest.dat)
m.w.pest.shan2 <- lmer(shan ~ pesticide + exclosure + census + (1|quad), data=Comm.w.pest.dat)

anova(m.w.pest.shan1,m.w.pest.shan2)

Comm.w.pest.dat$pesticide <- relevel(Comm.w.pest.dat$pesticide, ref="W")
m.w.pest.shan <- lmer(shan ~ pesticide + exclosure + census + (1|quad), data=Comm.w.pest.dat)

summary(m.w.pest.shan)

plot(m.w.pest.shan)
# this model is wrong, for the random effect should be site, rather than quad


#######################
### Site and quadrat
### It's site not code
gathered.w.pest.spec.abun  <- aggregate(abun ~ pesticide + exclosure + site + quadrat + census + species, 
                                        data=gathered.w.pest, FUN=sum)
# View(gathered.pest.spec.abun)

unique(gathered.pest.spec.abun$species)

### Make the data available for community analysis
gathered.w.pest.spec.abun1 <- unite(gathered.w.pest.spec.abun, "new", pesticide:census, sep = "_")

# CommData1 <- as.data.frame.matrix(with(gathered.pest.spec.abun1,table(new1,species)))
# View(CommData1)

gathered.w.pest.spec.abun2 <- spread(gathered.w.pest.spec.abun1, species, abun)
# View(gathered.pest.spec.abun2)
gathered.w.pest.spec.abun2[is.na(gathered.w.pest.spec.abun2)] <- 0

Comm.w.pest <- gathered.w.pest.spec.abun2[,-c(1)]
# CommPest1 <- subset(CommPest,,-c("pesticide", "exclosure", "census"))
# CommPest1 <- CommPest[,-c("pesticide", "exclosure", "census")]

###
w.pest.shan1 <- diversity(Comm.w.pest)
hist(w.pest.shan1)

### dataset
### This data can be used later
gathered.w.pest.spec.abun3 <- separate(gathered.w.pest.spec.abun2,new,into=c("pesticide","exclosure","site","quadrat","census"),"_")

###

Comm.w.pest.dat1 <- data.frame(gathered.w.pest.spec.abun3$pesticide, 
                               gathered.w.pest.spec.abun3$exclosure, 
                               gathered.w.pest.spec.abun3$site, 
                               gathered.w.pest.spec.abun3$quadrat, 
                               gathered.w.pest.spec.abun3$census, 
                               w.pest.shan1)

### rename columns
colnames(Comm.w.pest.dat1) <- c("pesticide","exclosure","site","quad","census","shan")
View(Comm.w.pest.dat1)
### Comm.w.pest.dat1 include shannon diversity for each quadrat in each census 

### fit in the model
m.w.pest.shan1 <- lmer(shan ~ pesticide*exclosure + census + (1|site), data=Comm.w.pest.dat1)
m.w.pest.shan2 <- lmer(shan ~ pesticide + exclosure + census + (1|site), data=Comm.w.pest.dat1)
m.w.pest.shan3 <- lmer(shan ~ pesticide + exclosure + census + (1+pesticide|site), data=Comm.w.pest.dat1)
m.w.pest.shan4 <- lmer(shan ~ (pesticide + exclosure + census)^2 + (1|site), data=Comm.w.pest.dat1)
m.w.pest.shan5 <- lmer(shan ~ (pesticide + exclosure + census)^2 + (1+pesticide|site), data=Comm.w.pest.dat1)

anova(m.w.pest.shan1, m.w.pest.shan2,m.w.pest.shan3,m.w.pest.shan4,m.w.pest.shan5)
# first model is the best

Comm.w.pest.dat1$pesticide <- relevel(Comm.w.pest.dat1$pesticide, ref="W")
m.w.pest.shan <- lmer(shan ~ pesticide*exclosure + census + (1|site), data=Comm.w.pest.dat1)

summary(m.w.pest.shan)
plot(m.w.pest.shan)
qqnorm(resid(m.w.pest.shan))



#########################
### include the growth.form
### Abundance for each species,quadrat,census and growth form

####################
### Need to figure out what factors induce the wired residuals plot

### Considering growth form is not good





####################################################################
### census15f is the first census, and this census just conducted at the begining of the experiment
### remove the first census
####################################################################
## gathered.pest

gathered.w.pest.new <- gathered.w.pest

gathered.w.pest.new$census[gathered.w.pest.new$census=="15fa"] <- NA

gathered.w.pest.new$abun <- as.numeric(gathered.w.pest.new$height !=0 & !is.na(gathered.w.pest.new$height))
gathered.w.pest.dive.new  <- aggregate(abun ~ pesticide + exclosure + site + quadrat + census + species, 
                                     data=gathered.w.pest.new, FUN=sum)
summary(gathered.w.pest.dive.new)


## Make the data available for community analysis
gathered.w.pest.dive.new1 <- unite(gathered.w.pest.dive.new, "new", pesticide:census, sep = "_")
gathered.w.pest.dive.new2 <- spread(gathered.w.pest.dive.new, species, abun)
gathered.w.pest.dive.new2[is.na(gathered.w.pest.dive.new2)] <- 0
# View(gathered.w.pest.dive.new2)

Comm.w.pest.new <- gathered.w.pest.dive.new2
# CommPest1 <- subset(CommPest,,-c("pesticide", "exclosure", "census"))
# CommPest1 <- CommPest[,-c("pesticide", "exclosure", "census")]
Comm.w.pest.new1 <- Comm.w.pest.new[,-c(1:5)]
### Shannon index
w.pest.shan.new <- diversity(Comm.w.pest.new1)
### Simpson index
w.pest.simp.new <- diversity(Comm.w.pest.new1, "simpson")

hist(w.pest.shan.new)
hist(w.pest.simp.new)


### rename columns
w.pest.shan.new.data <- gathered.pest.dive.new2
w.pest.simp.new.data <- gathered.pest.dive.new2
### Shannon index
Comm.w.pest.shan.dat.new <- data.frame(w.pest.shan.new.data$pesticide, w.pest.shan.new.data$exclosure, 
                                    w.pest.shan.new.data$site, w.pest.shan.new.data$quadrat, 
                                    w.pest.shan.new.data$census, w.pest.shan.new)

### Simpson index
Comm.w.pest.simp.dat.new <- data.frame(w.pest.simp.new.data$pesticide, w.pest.simp.new.data$exclosure, 
                                    w.pest.simp.new.data$site, w.pest.simp.new.data$quadrat, 
                                    w.pest.simp.new.data$census, w.pest.simp.new)


### rename columns
colnames(Comm.w.pest.shan.dat.new) <- c("pesticide","exclosure","site","quadrat","census","shan")
colnames(Comm.w.pest.simp.dat.new) <- c("pesticide","exclosure","site","quadrat","census","simp")
### fit in the model
### Shannon
m.w.pest.shan.new1 <- lmer(shan ~ pesticide*exclosure + census + (1|site), 
                           data=Comm.w.pest.shan.dat.new)
m.w.pest.shan.new2 <- lmer(shan ~ pesticide + exclosure + census + (1|site), 
                           data=Comm.w.pest.shan.dat.new)
m.w.pest.shan.new3 <- lmer(shan ~ pesticide + census + (1|site), 
                           data=Comm.w.pest.shan.dat.new)
m.w.pest.shan.new4 <- lmer(shan ~ pesticide + (1|site), 
                           data=Comm.w.pest.shan.dat.new)

anova(m.w.pest.shan.new1, m.w.pest.shan.new2, m.w.pest.shan.new3, m.w.pest.shan.new4)

summary(m.w.pest.shan.new1)
plot(m.w.pest.shan.new1)
qqnorm(resid(m.w.pest.shan.new1))


### Simpson
m.w.pest.simp.new1 <- lmer(simp ~ pesticide*exclosure + census + (1|site), 
                          data=Comm.w.pest.simp.dat.new)
m.w.pest.simp.new2 <- lmer(simp ~ pesticide + exclosure + census + (1|site), 
                           data=Comm.w.pest.simp.dat.new)
anova(m.w.pest.simp.new1, m.w.pest.simp.new2)

drop1(m.w.pest.simp.new2)

m.w.pest.simp.new3 <- lmer(simp ~ pesticide + census + (1|site), 
                           data=Comm.w.pest.simp.dat.new)
anova(m.w.pest.simp.new2, m.w.pest.simp.new3)

m.w.pest.simp.new4 <- lmer(simp ~ pesticide*census + (1|site), 
                           data=Comm.w.pest.simp.dat.new)
anova(m.w.pest.simp.new3, m.w.pest.simp.new4)

summary(m.w.pest.simp.new3)
plot(m.w.pest.simp.new3)
qqnorm(resid(m.w.pest.simp.new3))



test.w.shan.pest1 <- data.frame(res = resid(m.w.pest.shan.new1), fit = fitted(m.w.pest.shan.new1), 
                              pest = Comm.w.pest.shan.dat.new$pesticide, excl = Comm.w.pest.shan.dat.new$exclosure,
                              cens = Comm.w.pest.shan.dat.new$census, quad = Comm.w.pest.shan.dat.new$quad)

qplot(x=fit, y=res, geom='smooth', data=test.w.shan.pest1, colour=pest, facets=~excl) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("between exclosure")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() 


qplot(x=fit, y=res, geom='smooth', data=test.w.shan.pest1, colour=pest, facets=~cens) + 
  xlab("fitted")+
  ylab("resi")+
  ggtitle("among census")+ theme_bw() +
  geom_hline(aes(yintercept=0)) + geom_point() 


####################
### bootstrap
####################


###################### shannon's index ######################
View(new.w.pest)
new.w.pest.shan <- expand.grid(pesticide = c("C", "F","I","W"), exclosure = c("0","1"),
                             census= c("16sp","16fa"))

bootfit.w.pest.shan <- bootMer(m.w.pest.shan, FUN=function(x)predict(x, new.w.pest.shan, re.form=~0), nsim=999)

new.w.pest.shan$lci <- apply(bootfit.w.pest.shan$t, 2, quantile, 0.025) 
new.w.pest.shan$uci <- apply(bootfit.w.pest.shan$t, 2, quantile, 0.975) 
new.w.pest.shan$shanpred <- predict(m.w.pest.shan, newdata= new.w.pest.shan, re.form=~0)

### plot
ggplot(new.w.pest.shan, aes(x = pesticide, y = shanpred, ymin = lci, ymax = uci,
                          color=as.factor(exclosure))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("census") + theme_bw() +
  ylab("Shannon's D")

ggplot(new.w.pest.shan, aes(x = census, y = shanpred, ymin = lci, ymax = uci,
                          color=as.factor(pesticide))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("exclosure") + theme_bw() +
  ylab("Shannon's D")


####
### 90% confidence intervals
new.w.pest.shan90 <- new.w.pest.shan
new.w.pest.shan90$lci <- apply(bootfit.w.pest.shan$t, 2, quantile, 0.01) 
new.w.pest.shan90$uci <- apply(bootfit.w.pest.shan$t, 2, quantile, 0.90) 
new.w.pest.shan90$shanpred <- predict(m.w.pest.shan, newdata= new.w.pest.shan, re.form=~0)

### plot
ggplot(new.w.pest.shan90, aes(x = census, y = shanpred, ymin = lci, ymax = uci,
                            color=as.factor(pesticide))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("exclosure") + theme_bw() +
  ylab("Shannon's D")






#########################
### Pielou's evenness
#########################
# Shannon index

J.w.pest <- w.pest.shan1/log(specnumber(Comm.w.pest))

hist(J.w.pest)

### evenness index data frame
Comm.w.pest.dat.even <- data.frame(gathered.w.pest.spec.abun3$pesticide, 
                                   gathered.w.pest.spec.abun3$exclosure, 
                                   gathered.w.pest.spec.abun3$site, 
                                   gathered.w.pest.spec.abun3$quadrat, 
                                   gathered.w.pest.spec.abun3$census, 
                                   J.w.pest)
# View(Comm.w.pest.dat.even)


### rename columns
colnames(Comm.w.pest.dat.even) <- c("pesticide","exclosure","site","quad","census","even")

### fit in the model
m.w.pest.even1 <- lmer(even ~ pesticide*exclosure + census + (1|site), 
                       data=Comm.w.pest.dat.even)
m.w.pest.even2 <- lmer(even ~ pesticide + exclosure + census + (1|site), 
                       data=Comm.w.pest.dat.even)
anova(m.w.pest.even1, m.w.pest.even2)
drop1(m.w.pest.even2)

summary(m.w.pest.even2)
plot(m.w.pest.even2)


####################
### bootstrap
####################
# View(new.w.pest.shan)
new.w.pest.even <- new.w.pest.shan

bootfit.w.pest.even <- bootMer(m.w.pest.even2, 
                               FUN=function(x)predict(x, new.w.pest.even, re.form=~0), 
                               nsim=999)

new.w.pest.even$lci <- apply(bootfit.pest.even$t, 2, quantile, 0.025) 
new.w.pest.even$uci <- apply(bootfit.pest.even$t, 2, quantile, 0.975) 
new.w.pest.even$evenpred <- predict(m.w.pest.even2, newdata= new.w.pest.even, re.form=~0)

### plot the original results
ggplot(new.w.pest.even, aes(x = pesticide, y = evenpred, ymin = lci, ymax = uci,
                          color=as.factor(exclosure))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("census") + theme_bw() +
  ylab("Evenness")

ggplot(new.w.pest.even, aes(x = census, y = evenpred, ymin = lci, ymax = uci,
                          color=as.factor(pesticide))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  facet_wrap("exclosure") + theme_bw() +
  ylab("Evenness")






##########################################################
######################## snow data set ###################
##########################################################
gathered.w.snow.new <- gathered.w.snow

### remove the first census
# gathered.snow.new$census[gathered.snow.new$census=="14f"] <- NA

gathered.w.snow.new$abun <- as.numeric(gathered.w.snow.new$height !=0 & !is.na(gathered.w.snow.new$height))
gathered.w.snow.dive.dat  <- aggregate(abun ~ treatment + site + quadrat + census + species, 
                                       data=gathered.w.snow.new, FUN=sum)
summary(gathered.w.snow.dive.dat)

## Make the data available for community analysis
gathered.w.snow.dive.dat1 <- unite(gathered.w.snow.dive.dat, "new", treatment:census, sep = "_")
gathered.w.snow.dive.dat2 <- spread(gathered.w.snow.dive.dat1, species, abun)
gathered.w.snow.dive.dat2[is.na(gathered.w.snow.dive.dat2)] <- 0
View(gathered.w.snow.dive.dat2)

Comm.w.snow.dat <- gathered.w.snow.dive.dat2[,-c(1)]
### Shannon index
w.snow.shan <- diversity(Comm.w.snow.dat)

### Simpson index
w.snow.simp <- diversity(Comm.w.snow.dat, "simpson")

hist(w.snow.shan)
hist(log(w.snow.shan))
hist(log10(w.snow.shan))


hist(w.snow.simp)
hist(log(w.snow.simp))
# both shann and simp are not normally distribution

### get data set of data frame for further models
gathered.w.snow.dive.dat3 <- separate(gathered.w.snow.dive.dat2,new,into=c("treatment","site","quadrat","census"),"_")
w.snow.shan.new.dat <- gathered.w.snow.dive.dat3


### Shannon index
Comm.w.snow.shan.dat <- data.frame(w.snow.shan.new.dat$treatment, w.snow.shan.new.dat$site, 
                                   w.snow.shan.new.dat$quadrat, w.snow.shan.new.dat$census, 
                                   w.snow.shan)

### rename columns
colnames(Comm.w.snow.shan.dat) <- c("treatment","site","quadrat","census","shan")

### fit in the model
# m.snow.shan <- lmer(log(shan) ~ treatment + census + (1|site), data=CommSnow.shan.dat)

# sum(CommSnow.shan.dat$shan != round(CommSnow.shan.dat$shan))  #How often do any non-integer values of shan occur in your data

m.w.snow.shan1 <- lmer(shan ~ treatment + census + (1|site), data=Comm.w.snow.shan.dat)
m.w.snow.shan2 <- lmer(shan ~ treatment + census + (1+treatment|site), data=Comm.w.snow.shan.dat)
m.w.snow.shan3 <- lmer(shan ~ treatment*census + (1|site), data=Comm.w.snow.shan.dat)
m.w.snow.shan4 <- lmer(shan ~ treatment*census + (1+treatment|site), data=Comm.w.snow.shan.dat)

anova(m.w.snow.shan1, m.w.snow.shan2, m.w.snow.shan3, m.w.snow.shan4)

m.w.snow.shan <- m.w.snow.shan4
summary(m.w.snow.shan)
plot(m.w.snow.shan1)
plot(m.w.snow.shan)



##################### Overdispersion? ######################
plot(m.w.snow.shan)
chisq.w.snow.shan <- sum(resid(m.w.snow.shan, type='pearson')^2)
chisq.w.snow.shan/df.residual(m.w.snow.shan)  ## 0.19, Less than 1
## significantly so?
1-pchisq(chisq.w.snow.shan, df.residual(m.w.snow.shan)) # 1,No

###### Normal distribution?
## qqplot
library(car)
qqPlot(resid(m.w.snow.shan))




###################### bootstrap: shannon's index ######################
new.w.snow.shan <- expand.grid(treatment = c("C", "S"), exclosure = c("0","1"),
                             census= c("15sp","15fa","16sp","16fa"))

bootfit.w.snow.shan <- bootMer(m.w.snow.shan, 
                               FUN=function(x)predict(x, new.w.snow.shan, re.form=~0), 
                               nsim=999)

new.w.snow.shan$lci <- apply(bootfit.w.snow.shan$t, 2, quantile, 0.025) 
new.w.snow.shan$uci <- apply(bootfit.w.snow.shan$t, 2, quantile, 0.975) 
new.w.snow.shan$shanpred <- predict(m.w.snow.shan, newdata= new.w.snow.shan, re.form=~0)

# plot
ggplot(new.w.snow.shan, aes(x = census, y = shanpred, ymin = lci, ymax = uci,
                          color=as.factor(treatment))) +
  geom_pointrange(position=pd) + 
  geom_point(alpha = 0.3, position=pd) +
  theme_bw() +
  ylab("Shannon's D")



####################
## no. of species ##
####################

gathered.w.snow.spec <- subset(gathered.w.snow,gathered.snow$height !=0 & !is.na(gathered.w.snow$height))
gathered.w.snow.spec$no.spec <-as.numeric(gathered.w.snow.spec$species)
gathered.w.snow.spri  <- aggregate(no.spec ~ treatment +  site + quadrat + census, 
                                 data=gathered.w.snow.spec, FUN=function(x)length(unique(x)))

hist(gathered.w.snow.spri$no.spec)

# model
m.w.snow.spri1 <- lmer(no.spec ~ treatment + census + (1|site),
                     data=gathered.w.snow.spri)
summary(m.w.snow.spri1)

plot(m.w.snow.spri1)
