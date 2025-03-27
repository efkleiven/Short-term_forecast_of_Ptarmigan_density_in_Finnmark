# script to plot 

# load jagsoutput
#....
getwd()
dir("./output")

#load 
load("output/HSD_reg&cluster_dd&rod2&DoYf&harv&rrtg&carc&anomaly&trend_00-16_v2.RData")
#out16 <- out16$sims.list

load("output/HSD_reg&cluster_dd&rod2&DoYf&harv&rrtg&carc&anomaly&trend_00-17_v2.RData")
out17 <- out17$sims.list

load("output/HSD_reg&cluster_dd&rod2&DoYf&harv&rrtg&carc&anomaly&trend_00-18_v2.RData")
out18 <- out18$sims.list

load("output/HSD_reg&cluster_dd&rod2&DoYf&harv&rrtg&carc&anomaly&trend_00-19_v2.RData")
out19 <- out19$sims.list

load("output/HSD_reg&cluster_dd&rod2&DoYf&harv&rrtg&carc&anomaly&trend_00-20_v2.RData")
out20 <- out20$sims.list

load("output/HSD_reg&cluster_dd&rod2&DoYf&harv&rrtg&carc&anomaly&trend_00-21_v2.RData")
out21 <- out21$sims.list

load("output/HSD_reg&cluster_dd&rod2&DoYf&harv&rrtg&carc&anomaly&trend_00-22_v2.RData")
out22 <- out22$sims.list

load("output/HSD_reg&cluster_dd&rod2&DoYf&harv&rrtg&carc&anomaly&trend_00-23_v2.RData")
out23 <- out23$sims.list

#--------------------------------------------------------------------------------------------
nyear <- 7
npar <- 10

dat_year <- tibble(value=c(out16$btDD-1, out16$btR, out16$btRD, out16$btDoYft, out16$btrrt, out16$bttgt, out16$btcarc, out16$btharvt, out16$btMOTHt, out16$bttrend,
                           out17$btDD-1, out17$btR, out17$btRD, out17$btDoYft, out17$btrrt, out17$bttgt, out17$btcarc, out17$btharvt, out17$btMOTHt, out17$bttrend,
                           out18$btDD-1, out18$btR, out18$btRD, out18$btDoYft, out18$btrrt, out18$bttgt, out18$btcarc, out18$btharvt, out18$btMOTHt, out18$bttrend,
                           out19$btDD-1, out19$btR, out19$btRD, out19$btDoYft, out19$btrrt, out19$bttgt, out19$btcarc, out19$btharvt, out19$btMOTHt, out19$bttrend,
                           out20$btDD-1, out20$btR, out20$btRD, out20$btDoYft, out20$btrrt, out20$bttgt, out20$btcarc, out20$btharvt, out20$btMOTHt, out20$bttrend,
                           out21$btDD-1, out21$btR, out21$btRD, out21$btDoYft, out21$btrrt, out21$bttgt, out21$btcarc, out21$btharvt, out21$btMOTHt, out21$bttrend,
                           out22$btDD-1, out22$btR, out22$btRD, out22$btDoYft, out22$btrrt, out22$bttgt, out22$btcarc, out22$btharvt, out22$btMOTHt, out22$bttrend,
                           out23$btDD-1, out23$btR, out23$btRD, out23$btDoYft, out23$btrrt, out23$bttgt, out23$btcarc, out23$btharvt, out23$btMOTHt, out23$bttrend),
                   par=rep(c(rep('DD', times=12000), rep('Rod t', times=12000), rep('Rod t-1', times=12000), rep('OnsetFall', times=12000), rep('Precipitation', times=12000), rep('Temp', times=12000), rep('Carcas', times=12000), rep('Harvest', times=12000), rep('Moth', times=12000), rep('Trend', times=12000)), times=nyear),
                   year=c( rep('2016', times=12000*npar), rep('2017', times=12000*npar), rep('2018', times=12000*npar), rep('2019', times=12000*npar),
                           rep('2020', times=12000*npar), rep('2021', times=12000*npar), rep('2022', times=12000*npar), rep('2023', times=12000*npar)) )

ggplot(dat_year, aes(y=value, x=year))+
  geom_violin()+
  facet_wrap(~par, ncol=5)+
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.75, size=16),
        axis.text.y = element_text(size=16),
        axis.title = element_blank(),
        strip.text = element_text(size = 16))

getwd()
ggsave("plot/cov_16&23.jpg",
       width=40, height = 20, units="cm", dpi=300)

#********************************
# make predictions
#********************************

# Prepare cov data
# load model data
load("BugsData.rds")

# make reg and clust vector
Reg <- BugsData$Reg
Clust <- BugsData$Clust

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

rr <- BugsData$rr
rrt <- BugsData$rrt
rrs <- BugsData$rrs

tg <- BugsData$tg
tgt <- BugsData$tgt
tgs <- BugsData$tgs

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
      mean(pm$btMOTH) * moth[Reg[s], SimT+1] + mean(pm$btMOTHt) * motht[SimT+1] + mean(pm$btMOTHs) * moths[Reg[s]] +
      mean(pm$btharv) * harv[s, SimT] + mean(pm$btharvt) * harvt[SimT] + mean(pm$btharvs) * harvs[s] +
      mean(pm$btrr) * rr[Reg[s], SimT+1] + mean(pm$btrrt) * rrt[SimT+1] + mean(pm$btrrs) * rrs[Reg[s]] + 
      mean(pm$bttg) * tg[Reg[s], SimT+1] + mean(pm$bttgt) * tgt[SimT+1] + mean(pm$bttgs) * tgs[Reg[s]] +
      mean(pm$btcarc) * carc[SimT] + mean(pm$bttrend) * (SimT)
  }
  
  mu <- exp(logMU)
  return(mu)
}

# make predictions
mu17 <- pred(pm=out16, SimT=17)
mu18 <- pred(pm=out17, SimT=18)
mu19 <- pred(pm=out18, SimT=19)
mu20 <- pred(pm=out19, SimT=20)
mu21 <- pred(pm=out20, SimT=21)
mu22 <- pred(pm=out21, SimT=22)
mu23 <- pred(pm=out22, SimT=23)
#mu24 <- pred(pm=out23, SimT=24)


#******************************
#* Plot predictions
#* ***************************


# prepare data for ggplot

mean <- colMeans(apply(exp(out23$mu), c(2,3), mean))
low <- apply(apply(exp(out23$mu), c(2,3), mean), 2, quantile, probs=0.025)
high <- apply(apply(exp(out23$mu),c(2,3), mean), 2, quantile, probs=0.975)
year <- 2000:2023

# store in tibble
dat <- tibble(mean=mean, low=low, high=high, year=year)

dat$pred <- c(rep(NA, times=17), mean(mu17), mean(mu18), mean(mu19), mean(mu20), mean(mu21), mean(mu22), mean(mu23) )

# make ggplot
ggplot(dat, aes(x = year))  + 
  geom_line(aes(y = mean), color = "#C2B391")+ 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, fill = "#C2B391") +
  geom_point(aes(y=pred, x=year))+
  scale_x_continuous(
    breaks = c(2000:2023),
    limits = c(2000, 2023)) + 
  scale_y_continuous(
    breaks = seq(0, 95, by=5),
    limits = c(0,95)) + 
  ylab("Average population density (mu)") +  
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.75))

ggsave("plot/Pt_predictions.png",
       width=30, height=18, units="cm")

#--------------------------------------------------------

pm <- out18
logMU <- array(NA,dim=c(315,1))
SimT<- 19
# extract last years mu
MUt <- apply(pm$mu[,,19], 2, mean)

# loop to make prediction
for(s in 1:315){
  logMU[s] <- #pm$bt0[Reg[s]] + 
    pm$rCl[Clust[s]]+ 
    mean(pm$btDD) * MUt[s] + mean(pm$btR) * rod[SimT+1, Reg[s]] + mean(pm$btRD) * rod[SimT, Reg[s]] + 
    mean(pm$btDoYf) * DoYf[Reg[s], SimT] + mean(pm$btDoYft) * DoYft[SimT] + mean(pm$btDoYfs) * DoYfs[Reg[s]] +
    mean(pm$btMOTH) * moth[Reg[s], SimT+1] + mean(pm$btMOTHt) * motht[SimT+1] + mean(pm$btMOTHs) * moths[Reg[s]] #+
    mean(pm$btharv) * harv[s, SimT] + mean(pm$btharvt) * harvt[SimT] + mean(pm$btharvs) * harvs[s] +
    mean(pm$btrr) * rr[Reg[s], SimT+1] + mean(pm$btrrt) * rrt[SimT+1] + mean(pm$btrrs) * rrs[Reg[s]] + 
    mean(pm$bttg) * tg[Reg[s], SimT+1] + mean(pm$bttgt) * tgt[SimT+1] + mean(pm$bttgs) * tgs[Reg[s]] +
    mean(pm$btcarc) * carc[SimT] + mean(pm$bttrend) * (SimT)
}

mu <- exp(logMU)
mu

hist(unlist(pm$bt0))
hist(pm$rCl)

mean(pm$btDD)*mean(MUt)

mean(pm$rCl)
