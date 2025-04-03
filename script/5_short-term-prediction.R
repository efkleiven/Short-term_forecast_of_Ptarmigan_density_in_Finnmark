#******************************************************************************************
# In this script we make short-term forecasts for the willow ptarmigan density in Finnmark
#*****************************************************************************************

# Prepare cov data
# load data
load("../data/data_JAGS_2024.rds")

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

# scale rr
BugsData$rrt <- standardize(colMeans(BugsData$rr))
BugsData$rrs <- standardize(rowMeans(BugsData$rr))
BugsData$rr <- standardize(BugsData$rr)

# scale tg
BugsData$tgt <- standardize(colMeans(BugsData$tg))
BugsData$tgs <- standardize(rowMeans(BugsData$tg))
BugsData$tg <- standardize(BugsData$tg)

# scale anom
BugsData$motht <- standardize(colMeans(BugsData$anom))
BugsData$moths <- standardize(rowMeans(BugsData$anom))
BugsData$moth <- standardize(BugsData$anom)

#scale carc
BugsData$carc <- standardize(BugsData$carc$carc)

# store covariates in vector
rod <- BugsData$rod

rod

# Function to fit AR(2) model to each column
fit_ar2 <- function(column) {
  arima(column, order = c(2, 0, 0))  # AR(2) model (p=2, d=0, q=0)
}

# Apply the AR(2) model to each column of the matrix
ar2_models <- apply(rod, 2, fit_ar2)

# Print model results for each column
ar2_models

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

# Predict the next value for each column
predictions <- sapply(1:ncol(rod), function(i) {
  predict_next_value(ar2_models[[i]], rod[, i])
})

# Print predictions
predictions

# Extend the original time series by appending the predictions
extended_matrix <-rbind(array(rod, dim=c(25,4)), predictions)

# Create a time axis for plotting (assuming time steps are sequential)
time_steps <- 2000:2025

# Plot the original time series with the predicted values
# Plot the original time series with the predicted values
matplot(time_steps, extended_matrix[, 1:ncol(rod)], type = "l", 
        lty = 1, col = 1:ncol(rod), xlab = "Time", 
        ylab = "Value", main = "Time Series with AR(2) Prediction")

# Add the predicted points with a different symbol (e.g., "pch = 19" for filled circles)
points(rep(max(time_steps), ncol(rod)), predictions, pch = 19, col = 1:ncol(rod))

# Add a legend to identify the lines
#legend("topright", legend = c(paste("Time Series", 1:ncol(rod)), "Predicted"), 
#       col = c(1:ncol(rod), "black"), lty = 1, cex = 0.8)

#-----------------------------------
# Function to predict two steps ahead using the AR(2) model coefficients
predict_two_steps_ahead <- function(model, column_data) {
  # Extract estimated coefficients (phi1, phi2)
  phi1 <- coef(model)[1]  # AR(1) coefficient
  phi2 <- coef(model)[2]  # AR(2) coefficient
  
  # Get the last two observed values in the column (y_t-1, y_t-2)
  y_t_1 <- tail(column_data, 1)      # Last value (y_t-1)
  y_t_2 <- tail(column_data, 2)[1]   # Second last value (y_t-2)
  
  # First prediction (y_t+1)
  y_t_1_pred <- phi1 * y_t_1 + phi2 * y_t_2
  
  # Second prediction (y_t+2), using y_t+1 and y_t-1
  y_t_2_pred <- phi1 * y_t_1_pred + phi2 * y_t_1
  
  return(c(y_t_1_pred, y_t_2_pred))  # Return both predictions
}

# Predict two steps ahead for each column
predictions_two_steps <- sapply(1:ncol(rod), function(i) {
  predict_two_steps_ahead(ar2_models[[i]], rod[, i])
})

# Extend the original time series by appending the predictions
extended_matrix <- rbind(rod, predictions_two_steps[1,], predictions_two_steps[2,])

# Create a time axis for plotting (assuming time steps are sequential)
time_steps <- 2000:2026

# Plot the original time series with the predicted values
matplot(time_steps, extended_matrix[, 1:ncol(rod)], type = "l", 
        lty = 1, col = 1:ncol(rod), xlab = "Time", 
        ylab = "Value", main = "Time Series with AR(2) Two-Step Prediction")

# Add the predicted points with a different symbol (e.g., "pch = 19" for filled circles)
points(rep(max(time_steps) - 1, ncol(rod)), predictions_two_steps[1,], pch = 19, col = 1:ncol(rod))
points(rep(max(time_steps), ncol(rod)), predictions_two_steps[2,], pch = 17, col = 1:ncol(rod))

#carc <- BugsData$carc

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
      mean(pm$btMOTH) * moth[Reg[s], SimT+1] + mean(pm$btMOTHt) * motht[SimT+1] + mean(pm$btMOTHs) * moths[Reg[s]] +
      mean(pm$btDoYf) * DoYf[Clust[s], SimT] + mean(pm$btDoYft) * DoYft[SimT] + mean(pm$btDoYfs) * DoYfs[Clust[s]] +
      mean(pm$btharv) * harv[s, SimT] + mean(pm$btharvt) * harvt[SimT] + mean(pm$btharvs) * harvs[s] +
      mean(pm$btcarc) * carc[SimT] + mean(pm$bttrend) * (SimT)
  }
  
  mu <- exp(logMU)
  return(mu)
}

# make predictions
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

dat$pred <- c(rep(NA, times=20), mean(mu20), mean(mu21), mean(mu22), mean(mu23) )

# make ggplot
ggplot(dat, aes(x = year))  + 
  geom_line(aes(y = mean), color = "#C2B391")+ 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, fill = "#C2B391") +
  geom_point(aes(y=pred, x=year))+
  scale_x_continuous(
    breaks = c(2000:2023),
    limits = c(2000, 2023)) + 
  scale_y_continuous(
    breaks = seq(0, 60, by=5),
    limits = c(0,60)) + 
  ylab("Average population density (mu)") +  
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.75))

ggsave("plot/Pt_predictions.png",
       width=30, height=18, units="cm")
