#
Sys.setlocale("LC_ALL", "no_NO.UTF-8")

# script to download data from GBIF/Living Norway
library(tidyverse)

# keys to line transect datasets
dataset_Keys <- c("b49a2978-0e30-4748-a99f-9301d17ae119", # Fjellstyrene
                  "6a948a1c-7e23-4d99-b1c1-ec578d0d3159", # Statskog
                  "c47f13c1-7427-45a0-9f12-237aad351040") # FeFo

# load library
library(LivingNorwayR)

# Download fefo data
Rype_arkiv <- getLNportalData("c47f13c1-7427-45a0-9f12-237aad351040",1.14)

# Core table
Core <- Rype_arkiv$getCoreTable()

# Complete event table
Eve_all <- tibble::as_tibble(Core$exportAsDataFrame()) %>%
  dplyr::mutate(sampleSizeValue = as.numeric(sampleSizeValue))

# Complete occurrence table and sort by locationID
Occ <- tibble::as_tibble(Rype_arkiv$getExtensionTables()[[1]]$exportAsDataFrame()) 

Eve <-OccEve <- Eve_all %>% 
  dplyr::mutate(eventDate = as.Date(eventDate)) %>%
  dplyr::mutate(Year = lubridate::year(eventDate)) %>%
  arrange(locationID, eventDate)

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
  dplyr::filter(scientificName == "Lagopus lagopus")%>%
  arrange(locationID, Year)

# check which year are in (if not newest year is in it, check that the version number is correct)
unique(d_obs$Year)

## Remove observations from transects with suspiciously short length (< 200) and duplicate transects
d_obs <- d_obs %>%
  filter(!(parentEventID %in% bad_transects))

## merge adults and chicks from same observation event
d_obs <- d_obs %>%
  group_by(eventID, Year, locality, DistanceToTransectLine, locationID, verbatimLocality) %>%
  summarise(indtot = sum(individualCount), chicks= sum(individualCount[lifeStage=="Juvenile"])) %>%
  arrange(locationID, Year)

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
clustname <- unique(d_obs$locality)
regname <- unique(d_obs$verbatimLocality)

d_obs <- d_obs %>%
  mutate(site=match(locationID, sitename),
         clust=match(locality, clustname),
         reg=match(verbatimLocality, regname))

d_trans <- d_trans %>%
  mutate(site=match(locationID, sitename),
         clust=match(locality, clustname),
         reg=match(verbatimLocality, regname))

# how many unique site?
max(d_obs$site)

# remove transectlines that has been counted less than 5 years
d_trans2 <- d_trans %>%
  group_by(locationID) %>%
  filter(n_distinct(Year) >= 5) %>%
  ungroup()

length(unique(d_trans2$locationID))

d_obs2 <- d_obs %>%
  semi_join(d_trans2, by = "locationID")
# this filtering reduced the number of transect lines from 339 to 190 and number of total observations from 5978 to 5423

#---------------------------------------------------------------
# Save d_trans for use in extracting weather variables 
save(d_trans2, file="./data/d_trans_2.rds")

# Save d_obs for summary of observations
save(d_obs2, file="./data/d_obs_2.rds")
#------------------------------------------------------------------------------------------------------------
#******** End download and processing GBIF data **********************************

