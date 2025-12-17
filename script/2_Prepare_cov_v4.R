# load covariates and prepare for JAGS
# v.4 - forecast made in early july
# updated 10.7.2025

# load info about transect and observations
load("data/d_trans_3.rds")
load("data/d_obs_3.rds")

d_obs <- d_obs3
d_trans <- d_trans3

dir("data")

# load spatial covariate (onsetS, max NDVI, OnsetF, Anomaly, rr, tg)
all_vars <- tibble(read.csv("data/all_vars_mean_2500m_buffer_July2025.csv"))
max(all_vars$year)

# manaually replacing names in d_trans to be in line with all_vars
d_trans$locality <- replace(d_trans$locality, d_trans$locality=="Karasjok Aitevarri","Aitevarri" )
d_trans$locality <- replace(d_trans$locality, d_trans$locality=="Karasjok Iskuras", "Iskuras" )
d_trans$locality <- replace(d_trans$locality, d_trans$locality=="Stilla - Joatka", "Stilla" )
d_trans$locality <- replace(d_trans$locality, d_trans$locality=="Laksfjordvidda/Kunes", "Kunes")
d_trans$locality <- replace(d_trans$locality, d_trans$locality=="Savngovann", "Sangovann")
d_trans$locality <- replace(d_trans$locality, d_trans$locality=="Nyborgmoen/Kroken", "Nyborgmoen")
d_trans$locality <- replace(d_trans$locality, d_trans$locality=="Ifjord-Lebesby", "Ifjord_Lebesby")
d_trans$locality <- replace(d_trans$locality, d_trans$locality=="Ifjord-Tana", "Ifjord_Tana")

# fix norwegian letters in all_vars
unique(all_vars$loc_coat)

all_vars$loc_coat <- replace(all_vars$loc_coat, all_vars$loc_coat=="Hjelms\xf8y","Hjelmsøy" )
all_vars$loc_coat <- replace(all_vars$loc_coat, all_vars$loc_coat=="Kj\xe6s", "Kjæs")
all_vars$loc_coat <- replace(all_vars$loc_coat, all_vars$loc_coat=="Rolvs\xf8y", "Rolvsøy")
all_vars$loc_coat <- replace(all_vars$loc_coat, all_vars$loc_coat=="S\xf8r\xf8ya", "Sørøya" )

all_sites <- unique(d_trans$locationID)
all_reg <- unique(d_trans$verbatimLocality)
all_clust <- unique(d_trans$locality)

dat <- d_trans %>%
  slice(match(all_sites, d_trans$locationID))%>%
  dplyr::select(locality, verbatimLocality, locationID, municipality, site)%>%
  arrange(locationID) %>%
  mutate(reg = match(verbatimLocality, all_reg), 
         clust = match(locality, all_clust) )

# add column with numeric cluster
all_vars$clust <- dat$clust[match(all_vars$loc_coat, dat$locality)]

# make array with mean OnsetF for each cluster
OnsetF <- all_vars %>%
  group_by(year, clust) %>%
  summarise(fall=mean(fall)) %>%
  pivot_wider(names_from = year, values_from = fall)%>%
  select(!clust)

# make empty array
OnsetF2 <- array(NA, dim = c(nrow(OnsetF), ncol(OnsetF)))

for(i in 1:nrow(OnsetF)){
  OnsetF2[i,] <- unlist(OnsetF[i,])
}

# plot onset of fall
plot(colMeans(OnsetF2)~c(2000:2024), type='b')

OnsetF2 <- OnsetF2[,-26]

# make array with mean rr for each cluster
rr <- all_vars %>%
  group_by(year, clust) %>%
  summarise(rr=mean(rr)) %>%
  pivot_wider(names_from = year, values_from = rr)%>%
  select(!clust)

# make empty array
rr2 <- array(NA, dim = c(nrow(rr), ncol(rr)))

for(i in 1:nrow(rr)){
  rr2[i,] <- unlist(rr[i,])
}


# make array with mean tg for each cluster
tg <- all_vars %>%
  group_by(year, clust) %>%
  summarise(tg=mean(tg)) %>%
  pivot_wider(names_from = year, values_from = tg)%>%
  select(!clust)

# make empty array
tg2 <- array(NA, dim = c(nrow(tg), ncol(tg)))

for(i in 1:nrow(tg)){
  tg2[i,] <- unlist(tg[i,])
}

# Assuming tg2 is your 27x26 matrix
# Each column is a year

# Step 1: Compute column means and standard deviations
year_means <- colMeans(tg2, na.rm = TRUE)
year_sd <- apply(tg2, 2, sd, na.rm = TRUE)

# Step 2: Create a basic plot
years <- 2000:2025  # Replace with actual years if available
plot(years, year_means, type = "b", pch = 19, ylim = range(c(year_means - year_sd, year_means + year_sd)),
     xlab = "Year", ylab = "Mean ± SD", main = "Yearly Mean tg")

# Step 3: Add error bars
arrows(years, year_means - year_sd, years, year_means + year_sd,
       angle = 90, code = 3, length = 0.05)

# plot rr

# Step 1: Compute column means and standard deviations
year_means <- colMeans(rr2, na.rm = TRUE)
year_sd <- apply(rr2, 2, sd, na.rm = TRUE)

# Step 2: Create a basic plot
years <- 2000:2025  # Replace with actual years if available
plot(years, year_means, type = "b", pch = 19, ylim = range(c(year_means - year_sd, year_means + year_sd)),
     xlab = "Year", ylab = "Mean ± SD", main = "Yearly Mean rr")

# Step 3: Add error bars
arrows(years, year_means - year_sd, years, year_means + year_sd,
       angle = 90, code = 3, length = 0.05)


# make array with mean anom for each cluster
anom <- all_vars %>%
  group_by(year, clust) %>%
  summarise(anom=mean(anom)) %>%
  pivot_wider(names_from = year, values_from = anom)%>%
  select(!clust)

# make empty array
anom2 <- array(NA, dim = c(nrow(anom), ncol(anom)))

for(i in 1:nrow(anom)){
  anom2[i,] <- unlist(anom[i,])
}

# replace NA with 0
anom2[is.na(anom2)] <- 0.01 # should I replace with 0 or mean(anom2, na.rm=T)?
 dim(anom2)

 # remove year with no data
anom2 <- anom2[,-26]
 
 scale(colMeans(anom2))
plot(colMeans(anom2))

# Step 1: Compute column means and standard deviations
year_means <- colMeans(anom2, na.rm = TRUE)
year_sd <- apply(anom2, 2, sd, na.rm = TRUE)

# Step 2: Create a basic plot
years <- 2000:2024  # Replace with actual years if available
plot(years, year_means, type = "b", pch = 19, ylim = range(c(year_means - year_sd, year_means + year_sd)),
     xlab = "Year", ylab = "Mean ± SD", main = "Yearly Mean anomaly")

# Step 3: Add error bars
arrows(years, year_means - year_sd, years, year_means + year_sd,
       angle = 90, code = 3, length = 0.05)

 # harvest 
harv <- tibble(read.delim("./data/HarvestTotal2000-2024.txt")) %>%
  mutate(harv = FelteWpt/Areal) %>%
  dplyr::select(!c(Region, Jaktdager, Areal, FelteWpt))

# 4 last years in Kvalsund is missing, replace with hammerfest
print(harv, n=1000)

# manually replace missing values
harv$harv[221:225] <- harv$harv[121:125] # replaced missing harv for kvalsund with hammerfest

# standardize names
unique(harv$Kommune)
unique(dat$municipality)

harv$Kommune[harv$Kommune=="Vadso"]  <- "Vadsø"
harv$Kommune[harv$Kommune=="Vardo"]  <- "Vardø"
harv$Kommune[harv$Kommune=="Sorvar"] <- "Sør-Varanger"
harv$Kommune[harv$Kommune=="Masoy"]  <- "Måsøy"

# make harvest wide format
harv2 <- harv %>%
  pivot_wider(names_from=Aar, values_from = harv) 

# make vector of positions 
pos_hr <-  match(dat$municipality, harv2$Kommune)

# remove kommune from harv
harv3 <- harv2 %>%
  dplyr::select(!(Kommune)) 

# make empty array
harv4 <- array(NA, dim = c(length(pos_hr), ncol(harv3)))

for(i in 1:length(pos_hr)){
  harv4[i,] <- unlist(harv3[pos_hr[i],])
}

#----------------------------------------------
carc <- tibble(read.csv2("data/TOTV30062024215135791.csv")) %>%
  filter(Skadetype=="Rein") %>% 
  mutate(year = year(as_datetime(Funnetdato, format="%d.%m.%Y")),
         month = month(as_datetime(Funnetdato, format="%d.%m.%Y")) )%>%
  filter(month %in% c(1:6) & year %in% c(2000:2024))%>%
  dplyr::select(year, Antall)%>%
  group_by(year) %>%
  summarise(carc = sum(Antall))

print(carc, n=50)

# checking that the formatting works
# compared to henden et al, and its similar
plot(carc, type='b')
dir("data")

library(readxl)
library(dplyr)
library(lubridate)
library(tibble)

carc25 <- read_excel("data/kadaver_2025.xlsx") %>%
  as_tibble() %>%
  filter(Skadetype == "Rein") %>% 
  mutate(
    year = year(Funnetdato),
    month = month(Funnetdato)
  )%>%
  filter(month %in% c(1:6) & year == 2025)%>%
  dplyr::select(year, Antall)%>%
  group_by(year) %>%
  summarise(carc = sum(Antall))

carc <- rbind(carc,carc25)
plot(carc, type='b')

#----------------------------
# Rodent data

# East
rodOst <- tibble(read.delim("data/storskala_04-25_spring.txt")) %>% 
  filter(season=="spring") %>%
  mutate(rod=Mruf)%>%
  dplyr::select(c(year, rod))%>%
  group_by(year)%>%
  summarise(rodO = sum(rod))

# replace the missing first years with the mean for all years in region east
rodOst2 <- tibble(year=2000:2003, rodO=rep(mean(rodOst$rodO), times=4))

# bind replace missing data to the real dataset
rodOst3 <- rbind(rodOst2, rodOst)
print(rodOst3, n=1000)

# make year a character
rodOst3$year <- as.character(rodOst3$year)

#Inner
rodInnerVest <- tibble(read.table("data/porsanger_mus_reg_2025.csv", sep=",", header=T)) %>%
  filter(seas=="SPRING")%>%
  dplyr::select(c(year,  rufoCold, rufoSold))

# make year a character
rodInnerVest$year <- as.character(rodInnerVest$year) 

# West  
rod1 <- rodOst3 %>% 
  right_join( rodInnerVest,  by=join_by(year))%>%
  left_join( rodOst3,   by=join_by(year))%>%
  dplyr::select(!year)%>% as.matrix()

# scale and put in array
# should scale columns seperatly
rod3 <- scale(array(as.numeric(rod1), dim=c(26,4)))

# Sample year labels (replace with your actual years if needed)
years <- 2000:2025  # Assuming 26 years

# Create the plot
matplot(years, rod1, type = "l", lty = 1, lwd = 2,
        col = 1:4, xlab = "Year", ylab = "Value",
        main = "Trends by Area")

# Add a legend
legend("topright", legend = paste("Area", 1:4), col = 1:4, lty = 1, lwd = 2)



#---------------
W=200
nD <- 8; delta=W/nD; midpt = seq(delta/2, W, delta)
N_obs <- dim(d_obs)[1]

#------------------------------------------
## Assembling all data in a list for JAGS
load("data/RypeData_GBIF00-24_v3.rds")
getwd()

input.data$OnsetF = OnsetF2
input.data$rr = rr2
input.data$tg =tg2
input.data$rr = rr2
input.data$anom=anom2
input.data$carc=carc
input.data$rod = rod3
input.data$harv = harv4

# save inputdata    
save(input.data, file = "data/data_JAGS_2024_v4.rds")
#- End of Script