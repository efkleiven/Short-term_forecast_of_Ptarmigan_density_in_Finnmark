# script to plot 

# load jagsoutput
#....
getwd()
dir("./output")


# load climate model
# --------------------
load("output/climatemodel_00-23.RData" )

# check convergence
hist(unlist(out23$Rhat))
max(unlist(out23$Rhat))
mean(unlist(out23$Rhat))
which(unlist(out23$Rhat)>1.1)

climmod <- out23$sims.list

# load foodweb model
#-------------------
load("output/foodwebmodel_00-23.RData" )

# check convergence
hist(unlist(out23$Rhat))
max(unlist(out23$Rhat))
mean(unlist(out23$Rhat))
which(unlist(out23$Rhat)>1.1)

# make traceplot where Rhat suggest issues
library(jagsUI)
traceplot(out23,parameters=c('btDoYft', 'btMOTHt', 'btrrt', 'bttgt', 'btR', 'btRD', 'btharvt', 'btcarc', 'bttrend') )
traceplot(out23,parameters=c('bt0') )

foodwebmod <- out23$sims.list

# load full model
#-----------------
load("output/fullmodel_00-23.RData" )

# check convergence
hist(unlist(out23$Rhat))
max(unlist(out23$Rhat))
mean(unlist(out23$Rhat))
which(unlist(out23$Rhat)>1.1)

fullmod <- out23$sims.list

#-----------------
# Prepare cov data
#-----------------

# load model data
load("data/RypeData_GBIF00-23_pred24_v1.rds")

# rename
BugsData <- input.data

# make reg and clust vector
Reg <- BugsData$Reg
Clust <- BugsData$Clust

# standardize predictor variables
library(AHMbook)
BugsData$DoYft <- standardize(colMeans(BugsData$OnsetF))
BugsData$DoYfs <- standardize(rowMeans(BugsData$OnsetF))
BugsData$DoYf <- standardize(BugsData$OnsetF)

# scale harv
BugsData$harvt <- standardize(colMeans(BugsData$harv))
BugsData$harvs <- standardize(rowMeans(BugsData$harv))
BugsData$harv <- standardize(BugsData$harv)

# scale rr
#BugsData$rrt <- standardize(colMeans(BugsData$rr))
#BugsData$rrs <- standardize(rowMeans(BugsData$rr))
#BugsData$rr <- standardize(BugsData$rr)

# scale tg
#BugsData$tgt <- standardize(colMeans(BugsData$tg))
#BugsData$tgs <- standardize(rowMeans(BugsData$tg))
#BugsData$tg <- standardize(BugsData$tg)

# scale anom
BugsData$motht <- standardize(colMeans(BugsData$anom))
BugsData$moths <- standardize(rowMeans(BugsData$anom))
BugsData$moth <- standardize(BugsData$anom)

#scale carc
BugsData$carc <- standardize(BugsData$carc$carc)

# plot covariates
years1 <- 2000:2024
years <- 2000:2023
years2 <- 2000:2022

par(mfrow=c(3,2))

# onsetfall
plot(colMeans(BugsData$OnsetF)~years, type='b', main = "Ankomst høst", ylab="Day Of The Year", xlab="")

#harvest
plot(colMeans(BugsData$harv)~years, type='b', ylab="Felte ryper per km^2", xlab="", main="Høsting")

#percipitation
#plot(colMeans(BugsData$rr)~years, type='b', ylim=c(0,55), ylab="Nedbør 1. uke i Juli", xlab="", main="Nedbør")

#temp
#plot(colMeans(input.data$tg)~years, type='b', ylim=c(6,18), ylab="Temperatur 1. uke i Juli", xlab="", main="Temperatur")

#kadaver
plot(BugsData$carc~years1, type='b', ylab="Antall kadaver", xlab="", main="Kadaver")

#Bjørkemåler
plot(colMeans(BugsData$moth)~years, type='b', ylim=c(-1,1), ylab="Anomaly", xlab="", main="Bjørkemåler")

# plot rodents, for rodent cov data from east are artificially low since not all areas where trapped
plot(BugsData$rod[,1] ~ c(2000:2024), type='l', ylab='rodents', xlab='year', lwd=2, main="Smågnagere")
lines(BugsData$rod[,2]~ c(2000:2024) , col=2, lwd=2)
lines(BugsData$rod[,3]~ c(2000:2024), col=3, lwd=2)

# Inits
SimT = BugsData$T

# make covariate
rod <- BugsData$rod
carc <- BugsData$carc

DoYf <- BugsData$DoYf
DoYft <- BugsData$DoYft
DoYfs <- BugsData$DoYfs

harv <- BugsData$harv
harvt <- BugsData$harvt
harvs <- BugsData$harvs

moth <- BugsData$moth
motht <- BugsData$motht
moths <- BugsData$moths

#rr <- BugsData$rr
#rrt <- BugsData$rrt
#rrs <- BugsData$rrs

#tg <- BugsData$tg
#tgt <- BugsData$tgt
#tgs <- BugsData$tgs

# write function to make predictions
pred <- function(pm, SimT){
  # make empty array
  logMU <- array(NA,dim=c(315,1))
  
  # extract last years mu
  MUt <- apply(pm$mu[,,SimT], 2, mean)
  
  # loop to make prediction
  for(s in 1:315){
    logMU[s] <- pm$bt0[Reg[s]] + 
      pm$rCl[Clust[s]]+ 
      mean(pm$btDD) * MUt[s] + mean(pm$btR) * rod[SimT+1, Reg[s]] + mean(pm$btRD) * rod[SimT, Reg[s]] + 
      mean(pm$btDoYf) * DoYf[Reg[s], SimT] + mean(pm$btDoYft) * DoYft[SimT] + mean(pm$btDoYfs) * DoYfs[Reg[s]] +
      mean(pm$btMOTH) * moth[Reg[s], SimT] + mean(pm$btMOTHt) * motht[SimT] + mean(pm$btMOTHs) * moths[Reg[s]] +
      mean(pm$btharv) * harv[s, SimT] + mean(pm$btharvt) * harvt[SimT] + mean(pm$btharvs) * harvs[s] +
      #mean(pm$btcarc) * carc[SimT]+
      mean(pm$bttrend) * (SimT)
  }
  
  mu <- exp(logMU)
  return(mu)
}

#*************************
# Climate model prediction
#-------------------------

# extract last years mu
MUt <- apply(climmod$mu[,,SimT], 2, mean)

logMU <- array(NA,dim=c(315,1))

# loop to make prediction
for(s in 1:315){
  logMU[s] <- mean(climmod$bt0[,Reg[s]]) + 
    mean(climmod$rCl[,Clust[s]])+ 
    #mean(climmod$btDD) * MUt[s] + mean(climmod$btR) * rod[SimT+1, Reg[s]] + mean(climmod$btRD) * rod[SimT, Reg[s]] + 
    mean(climmod$btDoYf) * DoYf[Clust[s], SimT] + mean(climmod$btDoYft) * DoYft[SimT] + mean(climmod$btDoYfs) * DoYfs[Clust[s]] +
    #mean(climmod$btMOTH) * moth[Clust[s], SimT] + mean(climmod$btMOTHt) * motht[SimT] + mean(climmod$btMOTHs) * moths[Clust[s]] +
    #mean(climmod$btharv) * harv[s, SimT] + mean(climmod$btharvt) * harvt[SimT] + mean(climmod$btharvs) * harvs[s] +
    #mean(climmod$btcarc) * carc[SimT]+
    mean(climmod$bttrend) * (SimT)
}

mu24 <- exp(mean(logMU))

# prepare data for ggplot

mean <- c(colMeans(apply(exp(climmod$mu), c(2,3), mean)),NA)
low <- c(apply(apply(exp(climmod$mu), c(2,3), mean), 2, quantile, probs=0.025),NA)
high <- c(apply(apply(exp(climmod$mu),c(2,3), mean), 2, quantile, probs=0.975),NA)
year <- 2000:2024

# store in tibble
dat <- tibble(mean=mean, low=low, high=high, year=year)

# add prediction
dat$pred <- c(rep(NA, times=24), mu24)
dat$pred_low <- c(rep(NA, times=24), mu24-sd(exp(logMU)))
dat$pred_high <- c(rep(NA, times=24), mu24+sd(exp(logMU)))

# make ggplot
ggplot(dat, aes(x = year))  + 
  geom_line(aes(y = mean), color = "#C2B391")+ 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, fill = "#C2B391") +
  geom_point(aes(y=pred, x=year))+
  geom_errorbar(aes(ymin = pred_low, ymax = pred_high, x=year), height = 0.5) +
  scale_x_continuous(
    breaks = c(2000:2024),
    limits = c(2000, 2024)) + 
  scale_y_continuous(
    breaks = seq(0, 40, by=5),
    limits = c(0,40)) + 
  ggtitle('Climate model')+
  ylab("Average population density (mu)") +  
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.75),
        title=element_text(size=16))

ggsave("plot/climatemod_predictions.png",
       width=30, height=18, units="cm")


#*************************
# Foodweb model prediction
#-------------------------

# extract last years mu
MUt <- apply(foodwebmod$mu[,,SimT], 2, mean)

logMU <- array(NA,dim=c(315,1))

# loop to make prediction
for(s in 1:315){
  logMU[s] <- mean(foodwebmod$bt0[,Reg[s]]) + 
    mean(foodwebmod$rCl[,Clust[s]])+ 
    mean(foodwebmod$btDD) * MUt[s] + mean(foodwebmod$btR) * rod[SimT+1, Reg[s]] + mean(foodwebmod$btRD) * rod[SimT, Reg[s]] + 
    #mean(foodwebmod$btDoYf) * DoYf[Clust[s], SimT] + mean(foodwebmod$btDoYft) * DoYft[SimT] + mean(foodwebmod$btDoYfs) * DoYfs[Clust[s]] +
    mean(foodwebmod$btMOTH) * moth[Clust[s], SimT] + mean(foodwebmod$btMOTHt) * motht[SimT] + mean(foodwebmod$btMOTHs) * moths[Clust[s]] +
    mean(foodwebmod$btharv) * harv[s, SimT] + mean(foodwebmod$btharvt) * harvt[SimT] + mean(foodwebmod$btharvs) * harvs[s] +
    mean(foodwebmod$btcarc) * carc[SimT]+
    mean(foodwebmod$bttrend) * (SimT)
}

mu24 <- mean(exp(logMU))

# prepare data for ggplot

mean <- c(colMeans(apply(exp(foodwebmod$mu), c(2,3), mean)),NA)
low <- c(apply(apply(exp(foodwebmod$mu), c(2,3), mean), 2, quantile, probs=0.025),NA)
high <- c(apply(apply(exp(foodwebmod$mu),c(2,3), mean), 2, quantile, probs=0.975),NA)
year <- 2000:2024

# store in tibble
dat <- tibble(mean=mean, low=low, high=high, year=year)

# add prediction
dat$pred <- c(rep(NA, times=24), mu24)
dat$pred_low <- c(rep(NA, times=24), mu24-sd(exp(logMU)))
dat$pred_high <- c(rep(NA, times=24), mu24+sd(exp(logMU)))

# make ggplot
ggplot(dat, aes(x = year))  + 
  geom_line(aes(y = mean), color = "#C2B391")+ 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, fill = "#C2B391") +
  geom_point(aes(y=pred, x=year))+
  geom_errorbar(aes(ymin = pred_low, ymax = pred_high, x=year), height = 0.5) +
  scale_x_continuous(
    breaks = c(2000:2024),
    limits = c(2000, 2024)) + 
  scale_y_continuous(
    breaks = seq(0, 50, by=5),
    limits = c(0,50)) + 
  ggtitle('Foodweb model')+
  ylab("Average population density (mu)") +  
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.75))

ggsave("plot/foowebmodel_predictions.png",
       width=30, height=18, units="cm")

#--------------------------------------------------------

#*************************
# Full model prediction
#-------------------------

# extract last years mu
MUt <- apply(fullmod$mu[,,SimT], 2, mean)

logMU <- array(NA,dim=c(315,1))

# loop to make prediction
for(s in 1:315){
  logMU[s] <- mean(fullmod$bt0[,Reg[s]]) + 
    mean(fullmod$rCl[,Clust[s]])+ 
    mean(fullmod$btDD) * MUt[s] + mean(fullmod$btR) * rod[SimT+1, Reg[s]] + mean(fullmod$btRD) * rod[SimT, Reg[s]] + 
    mean(fullmod$btDoYf) * DoYf[Clust[s], SimT] + mean(fullmod$btDoYft) * DoYft[SimT] + mean(fullmod$btDoYfs) * DoYfs[Clust[s]] +
    mean(fullmod$btMOTH) * moth[Clust[s], SimT] + mean(fullmod$btMOTHt) * motht[SimT] + mean(fullmod$btMOTHs) * moths[Clust[s]] +
    mean(fullmod$btharv) * harv[s, SimT] + mean(fullmod$btharvt) * harvt[SimT] + mean(fullmod$btharvs) * harvs[s] +
    mean(fullmod$btcarc) * carc[SimT]+
    mean(fullmod$bttrend) * (SimT)
}

mu24 <- mean(exp(logMU))

# prepare data for ggplot

mean <- c(colMeans(apply(exp(fullmod$mu), c(2,3), mean)),NA)
low <- c(apply(apply(exp(fullmod$mu), c(2,3), mean), 2, quantile, probs=0.025),NA)
high <- c(apply(apply(exp(fullmod$mu),c(2,3), mean), 2, quantile, probs=0.975),NA)
year <- 2000:2024

# store in tibble
dat <- tibble(mean=mean, low=low, high=high, year=year)
dat$pred_low <- c(rep(NA, times=24), mu24-sd(exp(logMU)))
dat$pred_high <- c(rep(NA, times=24), mu24+sd(exp(logMU)))

# add prediction
dat$pred <- c(rep(NA, times=24), mu24)
dat$pred_low <- c(rep(NA, times=24), mu24-sd(exp(logMU)))
dat$pred_high <- c(rep(NA, times=24), mu24+sd(exp(logMU)))

# make ggplot
ggplot(dat, aes(x = year))  + 
  geom_line(aes(y = mean), color = "#C2B391")+ 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, fill = "#C2B391") +
  geom_point(aes(y=pred, x=year))+
  geom_errorbar(aes(ymin = pred_low, ymax = pred_high, x=year), height = 0.5) +
  scale_x_continuous(
    breaks = c(2000:2024),
    limits = c(2000, 2024)) + 
  scale_y_continuous(
    breaks = seq(0, 55, by=5),
    limits = c(0,55)) + 
  ylab("Average population density (mu)") + 
  ggtitle('Full model')+
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.75),
        title=element_text(size=16))

ggsave("plot/fullmod_predictions.png",
       width=30, height=18, units="cm")

#--------------------------------------------------------

