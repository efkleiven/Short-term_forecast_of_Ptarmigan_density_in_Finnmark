# load covariates and prepare for JAGS

# load info about transect and observations
load("./data/d_trans_2.rds")
load("./data/d_obs_2.rds")

d_obs <- d_obs2
d_trans <- d_trans2

# load spatial covariate (onsetS, max NDVI, OnsetF, Anomaly, rr, tg)
all_vars <- tibble(read.csv("./data/all_vars.csv"))
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

all_vars$loc_coat <- replace(all_vars$loc_coat, all_vars$loc_coat=="Hjelms\xf8y", "HjelmsÃ¸y")
all_vars$loc_coat <- replace(all_vars$loc_coat, all_vars$loc_coat=="Kj\xe6s", "KjÃ¦s")
all_vars$loc_coat <- replace(all_vars$loc_coat, all_vars$loc_coat=="Rolvs\xf8y", "RolvsÃ¸y")
all_vars$loc_coat <- replace(all_vars$loc_coat, all_vars$loc_coat=="S\xf8r\xf8ya", "SÃ¸rÃ¸ya")

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

# checking that the formatting works
# compared to henden et al, and its similar
plot(carc, type='b')

#----------------------------

# Rodent data
rodOst <- tibble(read.delim("data/storskala_04-24_spring.txt")) %>% 
  filter(season=="spring") %>%
  mutate(rod=Mruf)%>%
  dplyr::select(c(year, rod))%>%
  group_by(year)%>%
  summarise(rodO = sum(rod))

# replace the missing first years with the mean for all years in region east
rodOst2 <- tibble(year=2000:2003, rodO=rep(mean(rodOst$rodO), times=4))

# bind replace missing data to the real dataset
rodOst3 <- rbind(rodOst2, rodOst)

# make year a character
rodOst3$year <- as.character(rodOst3$year)

#Inner
rodInner <- tibble(read.table("data/karmassumS.csv", sep=",", header=T)) %>%
  filter(seas=="SPRING")%>%
  dplyr::select(c(year,  rufocanus))

# make year a character
rodInner$year <- as.character(rodInner$year) 

# West  
rodwest <- tibble(read.table("data/karmassumC.csv", sep=",", header=T)) %>%
  filter(seas=="SPRING")%>%
  dplyr::select(c(year,  rufocanus))

# make year a character
rodwest$year <- as.character(rodwest$year) 

rod1 <- rodInner %>% 
  right_join( rodwest,  by=join_by(year))%>%
  right_join( rodOst3,  by=join_by(year))%>%
  left_join( rodOst3,   by=join_by(year))%>%
  dplyr::select(!year)%>% as.matrix()

# scale and put in array
rod3 <- scale(array(rod1, dim=c(25,4)))
#-----------------------------------------

W=200
nD <- 8; delta=W/nD; midpt = seq(delta/2, W, delta)
N_obs <- dim(d_obs)[1]

#------------------------------------------
## Assembling all data in a list for JAGS
load("data/RypeData_GBIF00-24.rds")

input.data$OnsetF = OnsetF2
input.data$rr = rr2
input.data$tg =tg2
input.data$rr = rr2
input.data$anom=anom2
input.data$carc=carc
input.data$rod = rod3
input.data$harv = harv4

# save inputdata    
save(input.data, file = "./data/data_JAGS_2024_v2.rds")
#- End of Script