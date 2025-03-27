# load libraries
library(tidyverse)

# look at files in the directory
dir() 

setwd("./data/met/processed")
filenames_rr <- dir()[1:24]
filenames_tg <- dir()[25:48]

# make empty lists
rrdata <- list()
tgdata <- list()

# make year vector
year <- 2000:2023

# import rr data
for(i in 1:length(filenames_rr)){
  rrdata[[i]] <- read.delim(filenames_rr[i], header = T, sep="")
  rrdata[[i]]$year <- year[i]
  rrdata[[i]]$site <- 1:315}

# merge list objects to one df
df1 <- do.call(rbind, rrdata) 

# format rr data 
df_rr <- df1 %>%
  tibble %>%
  dplyr::select(rr, year, site)%>%
  pivot_wider(names_from = year, values_from = rr)%>%
  dplyr::select(!site)

# import tg data
for(i in 1:length(filenames_tg)){
  tgdata[[i]] <- read.delim(filenames_tg[i], header = T, sep="")
  tgdata[[i]]$year <- year[i]
  tgdata[[i]]$site <- 1:315}

# merge list objects to one df
df2 <- do.call(rbind, tgdata) 

# format tg data 
df_tg <- df2 %>%
  tibble %>%
  dplyr::select(tg, year, site)%>%
  pivot_wider(names_from = year, values_from = tg)%>%
  dplyr::select(!site)

# set wd to save location
setwd("../../")

# save data
save(df_rr, file="rr.rds")
save(df_tg, file="tg.rds")

#-- end of script