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

# scale rr
BugsData$rrt <- standardize(colMeans(BugsData$rr))
BugsData$rrs <- standardize(rowMeans(BugsData$rr))
BugsData$rr <- standardize(BugsData$rr)

# scale tg
BugsData$tgt <- standardize(colMeans(BugsData$tg))
BugsData$tgs <- standardize(rowMeans(BugsData$tg))
BugsData$tg <- standardize(BugsData$tg)

#scale carc
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

rr <- BugsData$rr
rrt <- BugsData$rrt
rrs <- BugsData$rrs

tg <- BugsData$tg
tgt <- BugsData$tgt
tgs <- BugsData$tgs

# write function to make predictions
pred_full <- function(pm, SimT){
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
      mean(pm$btMOTH) * moth[Clust[s], SimT] + mean(pm$btMOTHt) * motht[SimT] + mean(pm$btMOTHs) * moths[Clust[s]] +
      mean(pm$btrr) * rr[Reg[s], SimT+1] + mean(pm$btrrt) * rrt[SimT+1] + mean(pm$btrrs) * rrs[Reg[s]] + 
      mean(pm$bttg) * tg[Reg[s], SimT+1] + mean(pm$bttgt) * tgt[SimT+1] + mean(pm$bttgs) * tgs[Reg[s]] +
      mean(pm$btcarc) * carc[SimT] +  mean(pm$bttrend) * (SimT)
  }
  
  mu <- exp(logMU)
  return(mu)
}

# load jags model
load("./output/fullmod_00-20.RData")
load("./output/fullmod_00-21.RData")
load("./output/fullmod_00-22.RData")
load("./output/fullmod_00-23.RData")
load("./output/fullmod_00-24.RData")

# make predictions
mu21 <- pred_full(pm=out_full_20, SimT=21)
mu22 <- pred_full(pm=out_full_21, SimT=22)
mu23 <- pred_full(pm=out_full_22, SimT=23)
mu24 <- pred_full(pm=out_full_23, SimT=24)
mu25 <- pred_full(pm=out_full_24, SimT=25)

#******************************
#* Plot predictions
#* ***************************

# prepare data for ggplot
#mean <- colMeans(apply(exp(out_full_24$mu), c(2,3), mean))
#low <- apply(apply(exp(out_full_24$mu), c(2,3), mean), 2, quantile, probs=0.025)
#high <- apply(apply(exp(out_full_24$mu),c(2,3), mean), 2, quantile, probs=0.975)

year <- 2000:2024

mean <- apply(exp(out_full_24$mean_mu), 2, mean)
low <- apply(exp(out_full_24$mean_mu), 2, function(x) quantile(x, probs = 0.025))
high <- apply(exp(out_full_24$mean_mu), 2, function(x) quantile(x, probs = 0.975))

# store in tibble
dat <- tibble(mean=mean, low=low, high=high, year=year, pred=NA)

dat2 <- rbind(dat, c(NA,NA,NA,2025, mean(mu25)))

dat2$pred[22] <- mean(mu21)
dat2$pred[23] <- mean(mu22)
dat2$pred[24] <- mean(mu23)
dat2$pred[25] <- mean(mu24)

print(dat2, n=100)

# standardfeil
se <- 2*sd(mu25)/sqrt(length(mu25))

# Error bar only for the second point (x = 2, y = 10)
error_df <- data.frame(
  year = 2025,
  mean = mean(mu25),
  ymin = mean(mu25)-se,
  ymax = mean(mu25)+se
)

error_df <- data.frame(
  year = 2025,
  mean = mean(mu25),
  ymin = mean(mu25)-sd(mu25),
  ymax = mean(mu25)+sd(mu25)
)
print(dat2, n=100)

# make plot for full predictions
plot_1 <- ggplot(dat2, aes(x = year))  + 
  geom_line(aes(y = mean), color = "#C2B391")+ 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, fill = "#C2B391") +
  geom_point(aes(y=pred, x=year), size=5, color="#0072B2")+
  #geom_errorbar(data = error_df, aes(ymin = ymin, ymax = ymax), width = 0.1)+
  #geom_errorbar(data = error_df, aes(ymin = ymin, ymax = ymax), width = 0.1)+
  scale_x_continuous(
    breaks = c(2009:2023),
    limits = c(2009, 2023)) + 
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

ggsave("plot_full_23.png", plot = plot_1,  width=30, height=18, units="cm")

# make plot for full predictions
plot_1 <- ggplot(dat2, aes(x = year))  + 
  geom_line(aes(y = mean), color = "#C2B391")+ 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, fill = "#C2B391") +
  geom_point(aes(y=pred, x=year), size=5, color="#0072B2")+
  #geom_errorbar(data = error_df, aes(ymin = ymin, ymax = ymax), width = 0.1)+
  #geom_errorbar(data = error_df, aes(ymin = ymin, ymax = ymax), width = 0.1)+
  scale_x_continuous(
    breaks = c(2009:2024),
    limits = c(2009, 2024)) + 
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

ggsave("plot_full_24.png", plot = plot_1,  width=30, height=18, units="cm")

#-------------------------------------------------------------------
# without extream values

# write function to make predictions
pred_24 <- function(pm, SimT){
  # make empty array
  logMU <- array(NA,dim=c(176,1))
  
  # extract last years mu
  MUt <- apply(pm$mu[,,SimT], 2, mean)
  
  # loop to make prediction
  for(s in 1:176){
    logMU[s] <- pm$bt0[Reg[s]] + 
      pm$rCl[Clust[s]]+ 
      mean(pm$btDD) * MUt[s] + mean(pm$btR) * rod[SimT+1, Reg[s]] + mean(pm$btRD) * rod[SimT, Reg[s]] + 
      #mean(pm$btDoYf) * DoYf[Clust[s], SimT] + mean(pm$btDoYft) * DoYft[SimT] + mean(pm$btDoYfs) * DoYfs[Clust[s]] + # remove Doyf as it was extreamly early
      mean(pm$btharv) * harv[s, SimT] + mean(pm$btharvt) * harvt[SimT] + mean(pm$btharvs) * harvs[s] +
      mean(pm$btMOTH) * moth[Clust[s], SimT] + mean(pm$btMOTHt) * motht[SimT] + mean(pm$btMOTHs) * moths[Clust[s]] +
      mean(pm$btrr) * rr[Reg[s], SimT+1] + mean(pm$btrrt) * rrt[SimT+1] + mean(pm$btrrs) * rrs[Reg[s]] + 
      mean(pm$bttg) * tg[Reg[s], SimT+1] + mean(pm$bttgt) * tgt[SimT+1] + mean(pm$bttgs) * tgs[Reg[s]] +
      mean(pm$btcarc) * carc[SimT] +  mean(pm$bttrend) * (SimT)
  }
  
  mu <- exp(logMU)
  return(mu)
}

pred_25 <- function(pm, SimT){
  # make empty array
  logMU <- array(NA,dim=c(176,1))
  
  # extract last years mu
  MUt <- apply(pm$mu[,,SimT], 2, mean)
  
  # loop to make prediction
  for(s in 1:176){
    logMU[s] <- pm$bt0[Reg[s]] + 
      pm$rCl[Clust[s]]+ 
      mean(pm$btDD) * MUt[s] + mean(pm$btR) * rod[SimT+1, Reg[s]] + mean(pm$btRD) * rod[SimT, Reg[s]] + 
      mean(pm$btDoYf) * DoYf[Clust[s], SimT] + mean(pm$btDoYft) * DoYft[SimT] + mean(pm$btDoYfs) * DoYfs[Clust[s]] + # remove Doyf as it was extreamly early
      mean(pm$btharv) * harv[s, SimT] + mean(pm$btharvt) * harvt[SimT] + mean(pm$btharvs) * harvs[s] +
      #mean(pm$btMOTH) * moth[Clust[s], SimT] + mean(pm$btMOTHt) * motht[SimT] + mean(pm$btMOTHs) * moths[Clust[s]] +
      mean(pm$btrr) * rr[Reg[s], SimT+1] + mean(pm$btrrt) * rrt[SimT+1] + mean(pm$btrrs) * rrs[Reg[s]] + 
      mean(pm$bttg) * tg[Reg[s], SimT+1] + mean(pm$bttgt) * tgt[SimT+1] + mean(pm$bttgs) * tgs[Reg[s]] +
      mean(pm$btcarc) * carc[SimT] +  mean(pm$bttrend) * (SimT)
  }
  
  mu <- exp(logMU)
  return(mu)
}

# make new prediction
mu24_2 <- pred_24(pm=out_full_23, SimT=24)
mu25_2 <- pred_25(pm=out_full_24, SimT=25)

# change new values
dat3 <- cbind(dat2, pred2=(c(rep(NA, times=24), mean(mu24_2), mean(mu25_2))))

# make plot for full predictions
plot_2 <- ggplot(dat3, aes(x = year))  + 
  geom_line(aes(y = mean), color = "#C2B391", size=2)+ 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, fill = "#C2B391") +
  geom_point(aes(y=pred, x=year), size=5, color="#0072B2")+
  geom_point(aes(y=pred2, x=year), size=5, color="#E69F00")+
  geom_segment(data = dat3[dat3$year == 2024, ], 
               aes(x = 2024, xend = 2025, 
                   y = mean, yend = mean), 
               linetype = "dashed", color = "#C2B391", size = 1.2)+
  #geom_errorbar(data = error_df, aes(ymin = ymin, ymax = ymax), width = 0.1)+
  #geom_errorbar(data = error_df, aes(ymin = ymin, ymax = ymax), width = 0.1)+
  scale_x_continuous(
    breaks = c(2009:2024),
    limits = c(2009, 2024)) + 
  scale_y_continuous(
    breaks = seq(0, 17, by=5),
    limits = c(0,17)) + 
  ylab("Ptarmigan abundance index") +  
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=14, angle = 45, vjust = 0.75),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        legend.position = "none")

ggsave("plot_extrem_24.png", plot = plot_2,  width=30, height=18, units="cm")

# make plot for full predictions
plot_2 <- ggplot(dat3, aes(x = year))  + 
  geom_line(aes(y = mean), color = "#C2B391", size=2)+ 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, fill = "#C2B391") +
  geom_point(aes(y=pred, x=year), size=5, color="#0072B2")+
  geom_point(aes(y=pred2, x=year), size=5, color="#E69F00")+
  geom_segment(data = dat3[dat3$year == 2024, ], 
               aes(x = 2024, xend = 2025, 
                   y = mean, yend = mean), 
               linetype = "dashed", color = "#C2B391", size = 1.2)+
  #geom_errorbar(data = error_df, aes(ymin = ymin, ymax = ymax), width = 0.1)+
  #geom_errorbar(data = error_df, aes(ymin = ymin, ymax = ymax), width = 0.1)+
  scale_x_continuous(
    breaks = c(2009:2025),
    limits = c(2009, 2025)) + 
  scale_y_continuous(
    breaks = seq(0, 17, by=5),
    limits = c(0,17)) + 
  ylab("Ptarmigan abundance index") +  
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=14, angle = 45, vjust = 0.75),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        legend.position = "none")

ggsave("plot_extrem_25.png", plot = plot_2,  width=30, height=18, units="cm")

#------------------------------------------------------------------
# add early forecasts

# incorporate rodent prediction

#--------------------------------------
# make forecast for rodents
#--------------------------------------

# Function to fit AR(2) model to each column
fit_ar2 <- function(column) {
  arima(column, order = c(2, 0, 0))  # AR(2) model (p=2, d=0, q=0)
}

# Function to predict next value using the AR(2) model coefficients
predict_next_value <- function(model, column_data) {
  # Extract estimated coefficients (phi1, phi2)
  phi1 <- coef(model)[1]  # AR(1) coefficient
  phi2 <- coef(model)[2]  # AR(2) coefficient
  
  # Get the last two observed values in the column (y_t-1, y_t-2)
  y_t_1 <- tail(column_data, 1)      # Last value (y_t-1)
  y_t_2 <- tail(column_data, 2)[1]   # Second last value (y_t-2)
  
  # Predict the next value (y_t)
  y_t_pred <- phi1 * y_t_1 + phi2 * y_t_2
  
  return(y_t_pred)
}

#make empty array for predictions
rod_pred <- array(NA, dim=c(26,4))

# 2020
# Apply the AR(2) model to each column of the matrix
ar2_21 <- apply(rod[1:22,], 2, fit_ar2)
ar2_22 <- apply(rod[1:23,], 2, fit_ar2)
ar2_23 <- apply(rod[1:24,], 2, fit_ar2)
ar2_24 <- apply(rod[1:25,], 2, fit_ar2)
ar2_25 <- apply(rod[1:26,], 2, fit_ar2)

# Predict the next value for each column
rod_pred_21 <- sapply(1:ncol(rod[1:22,]), function(i) {
  predict_next_value(ar2_21[[i]], rod[1:22, i])
})

rod_pred_22 <- sapply(1:ncol(rod[1:23,]), function(i) {
  predict_next_value(ar2_22[[i]], rod[1:23, i])
})

rod_pred_23 <- sapply(1:ncol(rod[1:24,]), function(i) {
  predict_next_value(ar2_23[[i]], rod[1:24, i])
})

rod_pred_24 <- sapply(1:ncol(rod[1:25,]), function(i) {
  predict_next_value(ar2_24[[i]], rod[1:25, i])
})

rod_pred_25 <- sapply(1:ncol(rod[1:26,]), function(i) {
  predict_next_value(ar2_25[[i]], rod[1:26, i])
})

# fill in values
rod_pred[22,] <- rod_pred_21 
rod_pred[23,] <- rod_pred_22 
rod_pred[24,] <- rod_pred_23 
rod_pred[25,] <- rod_pred_24 
rod_pred[26,] <- rod_pred_25 

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
      mean(pm$btDD) * MUt[s] + mean(pm$btR) * rod_pred[SimT+1, Reg[s]] + mean(pm$btRD) * rod[SimT, Reg[s]] + 
      mean(pm$btDoYf) * DoYf[Clust[s], SimT] + mean(pm$btDoYft) * DoYft[SimT] + mean(pm$btDoYfs) * DoYfs[Clust[s]] +
      mean(pm$btharv) * harv[s, SimT] + mean(pm$btharvt) * harvt[SimT] + mean(pm$btharvs) * harvs[s] +
      mean(pm$btcarc) * carc[SimT] + mean(pm$bttrend) * (SimT)
  }
  
  mu <- exp(logMU)
  return(mu)
}

# load jags model
load("./output/early_00-20.RData")
load("./output/early_00-21.RData")
load("./output/early_00-22.RData")
load("./output/early_00-23.RData")
load("./output/early_00-24.RData")

# make predictions
mu_early_21 <- pred_early(pm=out_early_20, SimT=21)
mu_early_22 <- pred_early(pm=out_early_21, SimT=22)
mu_early_23 <- pred_early(pm=out_early_22, SimT=23)
mu_early_24 <- pred_early(pm=out_early_23, SimT=24)
mu_early_25 <- pred_early(pm=out_early_24, SimT=25)


# change new values
dat4 <- cbind(dat3, pred_early=(c(rep(NA, times=21), mean(mu_early_21), mean(mu_early_22), mean(mu_early_23), mean(mu_early_24), mean(mu_early_25))))

# plot with full model predictions
plot3 <- ggplot(dat4, aes(x = year))  + 
  geom_line(aes(y = mean), color = "#C2B391", size=2)+ 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, fill = "#C2B391") +
  geom_point(aes(y=pred, x=year), size=4, color="#0072B2", shape=16)+
  #geom_point(aes(y=pred_early, x=year), size=4, color="#D55E00", shape=15)+
  #geom_point(aes(y=pred2, x=year), size=4, color="#E69F00", shape=17)+
  geom_segment(data = dat3[dat3$year == 2024, ], 
               aes(x = 2024, xend = 2025, 
                   y = mean, yend = mean), 
               linetype = "dashed", color = "#C2B391", size = 1.2)+
  #geom_errorbar(data = error_df, aes(ymin = ymin, ymax = ymax), width = 0.1)+
  #geom_errorbar(data = error_df, aes(ymin = ymin, ymax = ymax), width = 0.1)+
  scale_x_continuous(
    breaks = c(2013:2025),
    limits = c(2013, 2025)) + 
  scale_y_continuous(
    breaks = seq(0, 17, by=5),
    limits = c(0,17)) + 
  ylab("Ptarmigan abundance index") +  
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=14, angle = 45, vjust = 0.75),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        legend.position = "none")

ggsave("plot/Pt_predictions__full_2025.png",
       width=30, height=18, units="cm")

# make plot for full predictions
ggplot(dat4, aes(x = year))  + 
  geom_line(aes(y = mean), color = "#C2B391", size=2)+ 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, fill = "#C2B391") +
  geom_point(aes(y=pred, x=year), size=4, color="#0072B2", shape=16)+
  #geom_point(aes(y=pred_early, x=year), size=4, color="#D55E00", shape=15)+
  geom_point(aes(y=pred2, x=year), size=4, color="#E69F00", shape=17)+
  geom_segment(data = dat3[dat3$year == 2024, ], 
               aes(x = 2024, xend = 2025, 
                   y = mean, yend = mean), 
               linetype = "dashed", color = "#C2B391", size = 1.2)+
  #geom_errorbar(data = error_df, aes(ymin = ymin, ymax = ymax), width = 0.1)+
  #geom_errorbar(data = error_df, aes(ymin = ymin, ymax = ymax), width = 0.1)+
  scale_x_continuous(
    breaks = c(2013:2025),
    limits = c(2013, 2025)) + 
  scale_y_continuous(
    breaks = seq(0, 17, by=5),
    limits = c(0,17)) + 
  ylab("Ptarmigan abundance index") +  
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=14, angle = 45, vjust = 0.75),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        legend.position = "none")

ggsave("plot/Pt_predictions_extvalues_2025.png",
       width=30, height=18, units="cm")

# make plot for full predictions
plot_3 <- ggplot(dat4, aes(x = year))  + 
  geom_line(aes(y = mean), color = "#C2B391", size=2)+ 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, fill = "#C2B391") +
  geom_point(aes(y=pred, x=year), size=5, color="#0072B2", shape=16)+
  geom_point(aes(y=pred_early, x=year), size=5, color="#D55E00", shape=15)+
#  geom_point(aes(y=pred2, x=year), size=4, color="#E69F00", shape=17)+
  geom_segment(data = dat3[dat3$year == 2024, ], 
               aes(x = 2024, xend = 2025, 
                   y = mean, yend = mean), 
               linetype = "dashed", color = "#C2B391", size = 1.2)+
  #geom_errorbar(data = error_df, aes(ymin = ymin, ymax = ymax), width = 0.1)+
  #geom_errorbar(data = error_df, aes(ymin = ymin, ymax = ymax), width = 0.1)+
  scale_x_continuous(
    breaks = c(2013:2025),
    limits = c(2013, 2025)) + 
  scale_y_continuous(
    breaks = seq(0, 17, by=5),
    limits = c(0,17)) + 
  ylab("Ptarmigan abundance index") +  
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=14, angle = 45, vjust = 0.75),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        legend.position = "none")

ggsave("plot_early_25.png", plot = plot_3,  width=30, height=18, units="cm")

ggsave("plot/Pt_predictions_early_2025.png",
       width=30, height=18, units="cm")
# the end