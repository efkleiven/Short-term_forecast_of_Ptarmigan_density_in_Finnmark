#******************************************************************************************
# In this script we make short-term forecasts for the willow ptarmigan density in Finnmark
#*****************************************************************************************

# Prepare cov data
# load data
load("Short-term_forecast_of_Ptarmigan_density_in_Finnmark/data/data_JAGS_2024_v4.rds")

# rename
BugsData <- input.data

# make reg and clust vector
Reg <- BugsData$Reg
Clust <- BugsData$Clust

# Standardize predictor variables
library(AHMbook)
BugsData$DoYft <- standardize(colMeans(BugsData$OnsetF))
BugsData$DoYfs <- standardize(rowMeans(BugsData$OnsetF))
BugsData$DoYf <- standardize(BugsData$OnsetF)

# scale harv
BugsData$harvt <- standardize(colMeans(BugsData$harv))
BugsData$harvs <- standardize(rowMeans(BugsData$harv))
BugsData$harv <- standardize(BugsData$harv)

# scale anom
BugsData$motht <- standardize(colMeans(BugsData$anom))
BugsData$moths <- standardize(rowMeans(BugsData$anom))
BugsData$moth <- standardize(BugsData$anom)

BugsData$carc <- standardize(BugsData$carc$carc)

# store covariates in vector
rod <- BugsData$rod

carc <- BugsData$carc

harv <- BugsData$harv
harvt <- BugsData$harvt
harvs <- BugsData$harvs

moth <- BugsData$moth
motht <- BugsData$motht
moths <- BugsData$moths

DoYft <- BugsData$DoYft 
DoYfs <- BugsData$DoYfs 
DoYf <- BugsData$DoYf 

# write function to make predictions
pred_early <- function(pm, SimT){
  # make empty array
  logMU <- array(NA,dim=c(176,1))
  
  # extract last years mu
  MUt <- apply(pm$mu[,,SimT], 2, mean)
  
  # loop to make prediction
  for(s in 1:176){
    logMU[s] <- pm$bt0[Reg[s]] + 
      pm$rCl[Clust[s]]+ 
      mean(pm$btDD) * MUt[s] + mean(pm$btR) * rod[SimT+1, Reg[s]] + mean(pm$btRD) * rod[SimT, Reg[s]] + 
      mean(pm$btDoYf) * DoYf[Clust[s], SimT] + mean(pm$btDoYft) * DoYft[SimT] + mean(pm$btDoYfs) * DoYfs[Clust[s]] +
      mean(pm$btharv) * harv[s, SimT] + mean(pm$btharvt) * harvt[SimT] + mean(pm$btharvs) * harvs[s] +
      mean(pm$btcarc) * carc[SimT] + mean(pm$bttrend) * (SimT)
  }
  
  mu <- exp(logMU)
  return(mu)
}

# load jags model

load("./output/early_00-21.RData")
load("./output/early_00-22.RData")
load("./output/early_00-23.RData")
load("./output/early_00-24.RData")


# make predictions
mu22 <- pred_early(pm=out_early_21, SimT=22)
mu23 <- pred_early(pm=out_early_22, SimT=23)
mu24 <- pred_early(pm=out_early_23, SimT=24)
mu25 <- pred_early(pm=out_early_24, SimT=25)

#******************************
#* Plot predictions
#* ***************************

# prepare data for ggplot
#mean <- colMeans(apply(exp(out_full_24$mu), c(2,3), mean))
#low <- apply(apply(exp(out_full_24$mu), c(2,3), mean), 2, quantile, probs=0.025)
#high <- apply(apply(exp(out_full_24$mu),c(2,3), mean), 2, quantile, probs=0.975)

year <- 2000:2024

mean24 <- apply(exp(out_early_24$mean_mu), 2, mean)
low24 <- apply(exp(out_early_24$mean_mu), 2, function(x) quantile(x, probs = 0.025))
high24 <- apply(exp(out_early_24$mean_mu), 2, function(x) quantile(x, probs = 0.975))

# store in tibble
dat <- tibble(mean=mean, low=low, high=high, year=year, pred=NA)

# add row for 2025
dat2 <- rbind(dat, c(NA,NA,NA,2025, mean(mu25)))

# fill in predictions for 21, 22, 23
dat2$pred[23] <- mean(mu22)
dat2$pred[24] <- mean(mu23)
dat2$pred[25] <- mean(mu24)

# look at table
print(dat2, n=100)

# standardfeil
se <- 2*sd(mu25)/sqrt(length(mu25))

# Error bar only for the second point (x = 2, y = 10)
error_df <- data.frame(
  year = c(2023,2024,2025),
  mean = c(mean(mu23),mean(mu24),mean(mu25)),
  ymin = c(mean(mu23)-sd(mu23),mean(mu24)-sd(mu24),mean(mu25)-sd(mu25)),
  ymax = c(mean(mu23)+sd(mu23),mean(mu24)+sd(mu24),mean(mu25)+sd(mu25))
)

print(dat2, n=100)

# make ggplot
ggplot(dat2, aes(x = year))  + 
  geom_line(aes(y = mean), color = "#C2B391")+ 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, fill = "#C2B391") +
  geom_point(aes(y=pred, x=year), size=4, color="brown")+
  #geom_errorbar(data = error_df, aes(ymin = ymin, ymax = ymax), width = 0.1)+
  #geom_errorbar(data = error_df, aes(ymin = ymin, ymax = ymax), width = 0.1)+
  scale_x_continuous(
    breaks = c(2007:2025),
    limits = c(2007, 2025)) + 
  scale_y_continuous(
    breaks = seq(0, 17, by=5),
    limits = c(0,17)) + 
  ylab("Bestandsindeks for lirype") +  
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=14, angle = 45, vjust = 0.75),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        legend.position = "none")

#ggsave("plot/Pt_predictions_2025.png",
#       width=30, height=18, units="cm")
