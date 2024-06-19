# script to download data from GBIF/Living Norway
library(tidyverse)

# keys to line transect datasets
dataset_Keys <- c("b49a2978-0e30-4748-a99f-9301d17ae119", # Fjellstyrene
                  "6a948a1c-7e23-4d99-b1c1-ec578d0d3159", # Statskog
                  "c47f13c1-7427-45a0-9f12-237aad351040") # FeFo

# load library
library(LivingNorwayR)

# Download fefo data
Rype_arkiv <- getLNportalData("c47f13c1-7427-45a0-9f12-237aad351040",1.13)

# Core table
Core <- Rype_arkiv$getCoreTable()

# Complete event table
Eve_all <- tibble::as_tibble(Core$exportAsDataFrame()) %>%
  dplyr::mutate(sampleSizeValue = as.numeric(sampleSizeValue))

# Complete occurrence table
Occ <- tibble::as_tibble(Rype_arkiv$getExtensionTables()[[2]]$exportAsDataFrame())

Eve <- Eve_all %>% 
  dplyr::mutate(eventDate = as.Date(eventDate)) %>%
  dplyr::mutate(Year = lubridate::year(eventDate))

## Assemble transect level info
d_trans <- Eve %>% 
  dplyr::select(locationID, eventDate, eventID, modified, 
                samplingProtocol, eventRemarks, sampleSizeValue, 
                stateProvince, municipality, locality, 
                verbatimLocality, locationRemarks,footprintWKT, Year) %>%
  dplyr::filter(eventRemarks == "Line transect") %>%
  dplyr::mutate(locationRemarks = gsub("In the original data this is known as lineID ", '', locationRemarks)) %>%
  tidyr::separate(., col = locationRemarks, sep = ",", into = c("LineName", "locationRemarks")) %>%
  dplyr::select(-locationRemarks) 

## Identify and remove transects with suspiciously short length (< 200 m) and duplicate transects
duplTransects <- c(
  "3EAD98CE-C7C9-472D-9C79-48E735C4525D",
  "54F312E7-9E21-44C5-ADE3-90F47EC95F62",
  "B82C1256-D359-4075-B65B-4304CB086A4B",
  "443FF8B8-071C-4AEB-B834-3C43748D194F",
  "D7F9B1FA-75F9-43A7-9AD4-9479BB9B7C34",
  "20C98FA8-AB29-43D9-BC18-90973367B080",
  "05E5C283-74AB-4192-8E5F-E57370DDD527",
  "DA651143-0B57-44E1-8266-C9CC169391F2"
)

bad_transects <- c(d_trans$eventID[which(d_trans$sampleSizeValue < 200)], duplTransects)

d_trans <- d_trans %>%
  dplyr::filter(!(eventID %in% bad_transects))

## Double-check no duplicate transects remain
transect_duplicates <- d_trans %>%
  dplyr::mutate(year = lubridate::year(eventDate)) %>%
  dplyr::group_by(locationID, locality, verbatimLocality, year) %>%
  dplyr::summarise(transectCount = dplyr::n(), .groups = 'keep') %>%
  dplyr::filter(transectCount > 1)

transect_duplicates

## Assemble observation level info
# Observations: distance to transect lines
d_obsTemp <- Eve %>% 
  dplyr::select(locationID, locality, verbatimLocality, parentEventID, eventID, eventRemarks, 
                dynamicProperties, eventDate) %>%
  dplyr::filter(eventRemarks == "Human observation" & !is.na(dynamicProperties)) %>%
  dplyr::mutate(dynamicProperties = purrr::map(dynamicProperties, ~ jsonlite::fromJSON(.) %>% as.data.frame())) %>%
  tidyr::unnest(dynamicProperties) %>% 
  dplyr::rename(DistanceToTransectLine = "perpendicular.distance.in.meters.from.transect.line.as.reported.by.the.field.worker") %>%
  dplyr::mutate(DistanceToTransectLine = as.numeric(DistanceToTransectLine), 
                Year = lubridate::year(eventDate)) %>%
  dplyr::select(locationID, locality, verbatimLocality, parentEventID, eventID, DistanceToTransectLine, Year)

# Observations: remaining information (willow ptarmigan only)
d_obs <- Occ %>% 
  dplyr::filter(!is.na(eventID)) %>%
  dplyr::select(eventID, scientificName, individualCount, sex, lifeStage) %>%
#  dplyr::select(-sex, -lifeStage) %>%
  dplyr::right_join(., d_obsTemp, by = "eventID") %>%
  dplyr::filter(scientificName == "Lagopus lagopus")

# check which year are in (if not newest year is in it, check that the version number is correct)
unique(d_obs$Year)

## Remove observations from transects with suspiciously short length (< 200) and duplicate transects
d_obs <- d_obs %>%
  filter(!(parentEventID %in% bad_transects))

## merge adults and chicks from same observation event
d_obs

d_obs <- d_obs %>%
  group_by(eventID, Year, locality, DistanceToTransectLine, locationID, verbatimLocality) %>%
  summarise(indtot = sum(individualCount), chicks= sum(individualCount[lifeStage=="Juvenile"]))

# Remove info about transects that are not in observation dataset
which(d_trans$locationID %in% d_obs$locationID )

d_trans <- d_trans[which(d_trans$locationID %in% d_obs$locationID),]

#--------------------------------------------------------------
# manual merging of transect lines and correction of region

# merge Blåsenberg lines with Jarfjord
d_trans$locationID <- replace(d_trans$locationID,d_trans$locationID=="32CD5DB7-A97D-4591-9EC5-1F421E611A5C", "1C64965F-D590-4201-9C4B-79CFE732D30C") 
d_trans$locationID <- replace(d_trans$locationID,d_trans$locationID=="BFEE67A0-156E-4099-BCD7-B9C606DA8C63", "C9E68D8E-B7BE-49EF-95D7-6D265BCDAA46")
d_trans$locality   <- replace(d_trans$locality,d_trans$locality=="Blåsenborg", "Jarfjord")

d_obs$locationID <- replace(d_obs$locationID,d_obs$locationID=="32CD5DB7-A97D-4591-9EC5-1F421E611A5C", "1C64965F-D590-4201-9C4B-79CFE732D30C") 
d_obs$locationID <- replace(d_obs$locationID,d_obs$locationID=="BFEE67A0-156E-4099-BCD7-B9C606DA8C63", "C9E68D8E-B7BE-49EF-95D7-6D265BCDAA46")
d_obs$locality   <- replace(d_obs$locality,d_obs$locality=="Blåsenborg", "Jarfjord")

# change locality name for Komag lines near Byvann
d_trans$locality <- replace(d_trans$locality, d_trans$locationID=="73E642A8-13EC-485A-97E7-565695814A09", "Byvann")
d_obs$locality <- replace(d_obs$locality, d_obs$locationID=="73E642A8-13EC-485A-97E7-565695814A09", "Byvann")

# add column for site
sitename <- unique(d_obs$locationID)

d_obs <- d_obs %>%
  mutate(site=match(locationID, sitename))

max(d_obs$site)
#------------------------------------------------------------------------------------------------------------
#prepare data for distance analysis

# Truncation distance
W <- 600 # Chloe uses 200, John seems to use 600

sort(unique(d_trans$LineName))

## Count observations per area
obs_count <- d_obs %>% 
  dplyr::filter(between(DistanceToTransectLine, -0.1, W)) %>% 
  dplyr::group_by(locality) %>%
  dplyr::summarise(count = n())
#---------------------------------

unique(d_obs$locality)
unique(d_obs$locationID)

## Set up arrays for area-specific data
N_yearsTot <- (max(d_obs$Year) - min(d_obs$Year))+1
N_sUnits <- length(unique(d_obs$locationID))
N_sites <- N_years <- min_years <- max_years <- N_obs <- N_sumR_obs <- rep(NA, N_sUnits)

dist <- c(); Year_obs <- c(); zeros_dist <- c()
A <- matrix(NA, nrow=N_sUnits, ncol=N_yearsTot) 
#y <- Year_obs <- zeros_dist <- matrix(NA, nrow = N_sUnits, ncol = max(obs_count$count))
sumR_obs <- sumR_obs_year <- sumAd_obs<- matrix(NA, nrow = N_sUnits, ncol = max(obs_count$count))

L <- N_line_year <- N_J_line_year <- N_A_line_year <- array(0, dim = c(n_distinct(d_trans$locationID), N_yearsTot))

#N_a_line_year <- array(0, dim = c(N_sUnits, N_ageC, max(site_count$count), N_yearsTot))

  # Constants #
  #-----------#
  
  ## Numbers of sites
  N_sites <- n_distinct(d_obs$locationID)
  
  ## (Number of) years with data
  years_mon <- sort(unique(d_trans$Year)) - min(d_obs$Year) + 1 
  min_years <- min(years_mon)
  max_years <- max(years_mon)
  
  # Transect characteristics #
  #--------------------------#
  
  ## Transect lengths
  # extract lengths from d_trans
  TransLen <- d_trans %>% 
    dplyr::select(locationID, Year, sampleSizeValue) %>%
    reshape2::dcast(locationID~Year, value.var = "sampleSizeValue", sum) %>%
    arrange(locationID)
  
  # save translen for use in observation exploration
  #save(TransLen, file="data/trans_length.rds")
  
  # make seperate matrix from transect lengths
  TransLen_mat <- TransLen %>% dplyr::select(-locationID) %>% as.matrix()
  colnames(TransLen_mat) <- as.numeric(colnames(TransLen_mat)) - min(2000) + 1 # rename cols from year to numbers from 1
  
  # fill length array
  for(t in 1:N_yearsTot){
    
    if(t %in% years_mon){
      L[1:N_sites, t] <- TransLen_mat[1:N_sites, which(years_mon == t)]
    }
  }
  
  ## Total covered area 
  A[1:N_sUnits, 1:N_yearsTot] <- L[,]*W*2/1000000
  
  # Observation distance from transect #
  #------------------------------------#
  temp_dist <- d_obs %>% 
    dplyr::filter(between(DistanceToTransectLine, -0.1, W)) %>% 
    dplyr::select(Year, DistanceToTransectLine, site) %>%
    dplyr::mutate(Year2 = Year - (min_years) + 1)
  
  ## Number of observations
  N_obs <- length(temp_dist$DistanceToTransectLine)
  #N_obs <- nrow(d_obs)
  
  ## Distance to transect line
  dist[1:N_obs] <- temp_dist$DistanceToTransectLine
  
  ## Observation year
  Year_obs[1:N_obs] <- temp_dist$Year2
  
  ## Vector of 0's of same length as y
  zeros_dist[1:N_obs] <- rep(0, N_obs)
  
  # Number of birds/line #
  #-------------------------------------------#
  temp <- TransLen %>% dplyr::select(locationID)
  
  TaksObs <- d_obs %>% 
    filter(between(DistanceToTransectLine, -0.1, W)) %>% # 
    dplyr::mutate(cs = indtot) %>%
    reshape2::dcast(locationID~Year, value.var = "cs", sum) %>%
    dplyr::right_join(., temp, by = c("locationID" = "locationID")) %>%
    replace(., is.na(.), 0) %>%
    dplyr::arrange(locationID)
  
  TaksObs_mat <- TaksObs %>% dplyr::select(-locationID) %>% as.matrix() 
  colnames(TaksObs_mat) <- as.numeric(colnames(TaksObs_mat)) - min_years + 1
  
  years_mon <- as.numeric(colnames(TaksObs_mat)) # filtering observations on DistanceToTransectLine will in some cases remove all observations for a transect line within 1 year. This recalculates number of years with samples
  
  for(t in 1:N_yearsTot){
      N_line_year[1:N_sites, t] <- TaksObs_mat[1:N_sites, t]}
  
# converting to names from Henden model
  
  # define coefisintes for bind distance
  nD <- 24; delta=W/nD; midpt = seq(delta/2, W, delta)
  
  # define distance class
  dclass <- ifelse(temp_dist$DistanceToTransectLine<25, 1,
                   ifelse(temp_dist$DistanceToTransectLine<50, 2,
                          ifelse(temp_dist$DistanceToTransectLine<75, 3,
                                 ifelse(temp_dist$DistanceToTransectLine<100, 4,
                                        ifelse(temp_dist$DistanceToTransectLine<125, 5,
                                               ifelse(temp_dist$DistanceToTransectLine<150, 6,
                                                      ifelse(temp_dist$DistanceToTransectLine<175, 7,
                                                             ifelse(temp_dist$DistanceToTransectLine<200, 8,
                                                                    ifelse(temp_dist$DistanceToTransectLine<225, 9,
                                                                           ifelse(temp_dist$DistanceToTransectLine<250, 10,
                                                                                  ifelse(temp_dist$DistanceToTransectLine<275, 11,
                                                                                         ifelse(temp_dist$DistanceToTransectLine<300, 12,
                                                                                                ifelse(temp_dist$DistanceToTransectLine<325, 13,
                                                                                                       ifelse(temp_dist$DistanceToTransectLine<350, 14,
                                                                                                              ifelse(temp_dist$DistanceToTransectLine<375, 15,
                                                                                                                     ifelse(temp_dist$DistanceToTransectLine<400, 16,
                                                                                                                            ifelse(temp_dist$DistanceToTransectLine<425, 17,
                                                                                                                                   ifelse(temp_dist$DistanceToTransectLine<450, 18,
                                                                                                                                          ifelse(temp_dist$DistanceToTransectLine<475, 19,
                                                                                                                                                 ifelse(temp_dist$DistanceToTransectLine<500, 20,
                                                                                                                                                        ifelse(temp_dist$DistanceToTransectLine<525, 21,
                                                                                                                                                               ifelse(temp_dist$DistanceToTransectLine<550, 22,
                                                                                                                                                                      ifelse(temp_dist$DistanceToTransectLine<575, 23, 24
                                                                                                                                                                      )))))))))))))))))))))))
  # decide site for each observation
  all_sites <- unique(d_trans$locationID)

d_trans <- d_trans %>%
  mutate(site = match(locationID, all_sites))

all_reg <- unique(d_trans$verbatimLocality)
all_clust <- unique(d_trans$locality)

#----------------------------------------------
dat <- d_trans %>%
  slice(match(all_sites, d_trans$locationID))%>%
  dplyr::select(locality, verbatimLocality, locationID, municipality, site)%>%
  arrange(locationID) %>%
  mutate(reg = match(verbatimLocality, all_reg), 
         clust = match(locality, all_clust) )

# check number of regions and clusters
length(unique(dat$locality))
length(unique(dat$verbatimLocality))

# check number of lines in each region
table(dat$verbatimLocality)

# check observations per region
table(d_obs$verbatimLocality)

# Save d_trans for use in extracting weather variables 
save(d_trans, file="data/d_trans2.rds")

# Save d_obs for summary of observations
save(d_obs, file="data/d_obs2.rds")

#******** End download and processing GBIF data **********************************

#******** Start Load and prepare Covariates **************************************

# load spatial covariate (onsetS, max NDVI, OnsetF, Anomaly, rr, tg)
all_vars<-tibble(read.csv("data/all_vars.csv"))
#load('data/d_trans.rds')

# manaually replacing names in d_trans
d_trans$locality <- replace(d_trans$locality, d_trans$locality=="Karasjok Aitevarri","Aitevarri" )
d_trans$locality <- replace(d_trans$locality, d_trans$locality=="Karasjok Iskuras", "Iskuras" )
d_trans$locality <- replace(d_trans$locality, d_trans$locality=="Stilla – Joatka", "Stilla" )
d_trans$locality <- replace(d_trans$locality, d_trans$locality=="Laksfjordvidda/Kunes", "Kunes")
d_trans$locality <- replace(d_trans$locality, d_trans$locality=="Savngovann", "Sangovann")
d_trans$locality <- replace(d_trans$locality, d_trans$locality=="Nyborgmoen/Kroken", "Nyborgmoen")
d_trans$locality <- replace(d_trans$locality, d_trans$locality=="Ifjord-Lebesby", "Ifjord_Lebesby")
d_trans$locality <- replace(d_trans$locality, d_trans$locality=="Ifjord-Tana", "Ifjord_Tana")

# fix norwegian letters in all_vars
unique(all_vars$loc_coat)

all_vars$loc_coat <- replace(all_vars$loc_coat, all_vars$loc_coat=="Hjelms\xf8y", "Hjelmsøy")
all_vars$loc_coat <- replace(all_vars$loc_coat, all_vars$loc_coat=="Kj\xe6s", "Kjæs")
all_vars$loc_coat <- replace(all_vars$loc_coat, all_vars$loc_coat=="Rolvs\xf8y", "Rolvsøy")
all_vars$loc_coat <- replace(all_vars$loc_coat, all_vars$loc_coat=="S\xf8r\xf8ya", "Sørøya")

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
all_vars$reg <- dat$clust[match(all_vars$loc_coat, dat$locality)]

# make array with mean OnsetF for each cluster
OnsetF <- all_vars %>%
  group_by(year, reg) %>%
  summarise(fall=mean(fall)) %>%
  pivot_wider(names_from = year, values_from = fall)%>%
  select(!reg)

# make empty array
OnsetF2 <- array(NA, dim = c(nrow(OnsetF), ncol(OnsetF)))

for(i in 1:nrow(OnsetF)){
  OnsetF2[i,] <- unlist(OnsetF[i,])
}

# make array with mean rr for each cluster
rr <- all_vars %>%
  group_by(year, reg) %>%
  summarise(rr=mean(rr)) %>%
  pivot_wider(names_from = year, values_from = rr)%>%
  select(!reg)

# make empty array
rr2 <- array(NA, dim = c(nrow(rr), ncol(rr)))

for(i in 1:nrow(rr)){
  rr2[i,] <- unlist(rr[i,])
}


# make array with mean tg for each cluster
tg <- all_vars %>%
  group_by(year, reg) %>%
  summarise(tg=mean(tg)) %>%
  pivot_wider(names_from = year, values_from = tg)%>%
  select(!reg)

# make empty array
tg2 <- array(NA, dim = c(nrow(tg), ncol(tg)))

for(i in 1:nrow(tg)){
  tg2[i,] <- unlist(tg[i,])
}

# make array with mean anom for each cluster
anom <- all_vars %>%
  group_by(year, reg) %>%
  summarise(anom=mean(anom)) %>%
  pivot_wider(names_from = year, values_from = anom)%>%
  select(!reg)

# make empty array
anom2 <- array(NA, dim = c(nrow(anom), ncol(anom)))

for(i in 1:nrow(anom)){
  anom2[i,] <- unlist(anom[i,])
}

# replace NA with 0
anom2[is.na(anom2)] <- 0

# harvest 
harv <- tibble(read.delim("data/HarvestTotal2000-2022.txt")) %>%
  mutate(harv = FelteWpt/Areal) %>%
  dplyr::select(!c(Region, Jaktdager, Areal, FelteWpt))

# 3 last years in Kvalsund is missing, replace with hammerfest
print(harv, n=1000)

# manually replace missing values
harv$harv[205:207] <- harv$harv[113:115] # replaced based on location after manual check
  
# standardize names
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
harv4 <- array(NA, dim = c(315, 23))

for(i in 1:length(pos_hr)){
  harv4[i,] <- unlist(harv3[pos_hr[i],])
}

#----------------------------------------------
# Carcass data
carc <- tibble(read.csv2("data/TOTV13062023130101337.csv")) %>%
  filter(Skadetype=="Rein") %>% 
  mutate(year = year(as_datetime(Funnetdato, format="%d.%m.%Y")),
         month = month(as_datetime(Funnetdato, format="%d.%m.%Y")) )%>%
  filter(month %in% c(1:6) & year %in% c(2000:2023))%>%
  dplyr::select(year, Antall)%>%
  group_by(year) %>%
    summarise(carc = sum(Antall))

# checking that the formatting works
# compared to henden et al, and its similar
plot(carc, type='b')

#------------------------------------------
  ## Assembling all data in a list for JAGS
  input.data <- list(
    B=W, # truncation distance
    nD=nD, # number of distance binds
    delta = delta, #width of distance binds
    midpt = midpt,
    dist = dist, # Distance to transect line for each individual observation
    #zeros_dist = zeros_dist, # Vector of 0's of same length as y
    nobs = N_obs, # Total number of observations
    y = N_line_year, # Number of birds observed per site per year
    #area = L, # Transect length per site and year
    T = N_yearsTot, # Max number of years with data
    nsites = N_sites, # Total number of monitored sites per area
    area = A, # Total covered area per year per site
    dclass = dclass, # distance class for all observations
    site = temp_dist$site, # site for each observation
    Reg = dat$reg,
    Clust = dat$clust,
    NClust = length(unique(dat$clust)),
#    VerbClust = dat$locality,
    OnsetF = OnsetF2,
    harv = harv4,
    rr = rr2,
    tg =tg2,
    carc=carc,
    anom=anom2
  )

# save inputdata    
save(input.data, file = "data/RypeData_GBIF00-23_v3.rds")

#--- End of script
