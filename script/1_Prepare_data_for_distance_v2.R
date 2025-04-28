#******** Start Load and prepare Covariates **************************************
# load libraries
library(tidyverse)

# make R read norwegian letters
Sys.setlocale("LC_ALL", "no_NO.UTF-8")

# load ptarmigan data
load("./data/d_trans_2.rds")
load("./data/d_obs_2.rds")

d_obs <- d_obs2
d_trans <- d_trans2

#prepare data for distance analysis

# Truncation distance
W <- 200 # Chloe uses 200, John seems to use 600

sort(unique(d_trans$LineName))

## Count observations per area
obs_count <- d_obs %>% 
  dplyr::filter(between(DistanceToTransectLine, -0.1, W)) %>% 
  dplyr::group_by(locality) %>%
  dplyr::summarise(count = n())

print(obs_count, n=30)
sum(obs_count$count)

#--------------------------------
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

# make maxtrix for number of birds per line
TaksObs_mat <- TaksObs %>% dplyr::select(-locationID) %>% as.matrix() 
colnames(TaksObs_mat) <- as.numeric(colnames(TaksObs_mat)) - min_years + 1

years_mon <- as.numeric(colnames(TaksObs_mat)) # filtering observations on DistanceToTransectLine will in some cases remove all observations for a transect line within 1 year. This recalculates number of years with samples

for(t in 1:N_yearsTot){
  N_line_year[1:N_sites, t] <- TaksObs_mat[1:N_sites, t]}

# converting to names from Henden model

# define coefisintes for bind distance
nD <- 8; delta=W/nD; midpt = seq(delta/2, W, delta)

# define distance class
dclass <- ifelse(temp_dist$DistanceToTransectLine<25, 1,
                 ifelse(temp_dist$DistanceToTransectLine<50, 2,
                        ifelse(temp_dist$DistanceToTransectLine<75, 3,
                               ifelse(temp_dist$DistanceToTransectLine<100, 4,
                                      ifelse(temp_dist$DistanceToTransectLine<125, 5,
                                             ifelse(temp_dist$DistanceToTransectLine<150, 6,
                                                    ifelse(temp_dist$DistanceToTransectLine<175, 7, 8)))))))
# make vector with unique sites
all_sites <- unique(d_trans$locationID)

# make tibble with transect info, one row per line
dat <- d_trans %>%
  slice(match(all_sites, d_trans$locationID))%>%
  dplyr::select(locality, verbatimLocality, locationID, municipality, site, clust, reg, footprintWKT)%>%
  arrange(locationID) 

dat

# check number of regions and clusters
length(unique(d_trans$locality))
length(unique(d_trans$verbatimLocality))

# check number of lines in each region
table(dat$verbatimLocality)

# check observations per region
table(d_obs$verbatimLocality)

# Filter rows where 'indtot' is less than or equal to 20
d_obs_filtered <- d_obs %>%
  filter(indtot <= 20)

# Histogram of number of adults seen in each observation
ggplot(d_obs_filtered, aes(x = indtot-chicks)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Number of adults per observation", x = "n. adults", y = "Frequency")

# Create a histogram of the remaining values in 'indtot'
ggplot(d_obs_filtered, aes(x = indtot-chicks)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  facet_wrap(~Year)+
  labs(title = "Number of adults per observation", x = "n. adults", y = "Frequency")

#-------------------------------------------
# checking observations with many individuals
# observations of a lot of individuals, overrepresented in Hjemls?y and Rolvs?y in 2023 and 2024.
d_obs[d_obs$indtot==59,]
d_obs[d_obs$indtot==50,]
d_obs[d_obs$indtot==46,]
d_obs[d_obs$indtot==40,]
d_obs[d_obs$indtot==30,]
d_obs[d_obs$indtot==29,]
d_obs[d_obs$indtot==28,]
d_obs[d_obs$indtot==25,]
d_obs[d_obs$indtot==24,]
d_obs[d_obs$indtot==23,]
d_obs[d_obs$indtot==22,]
d_obs[d_obs$indtot==21,]
d_obs[d_obs$indtot==20,]

#--------------------------------------------------
# calcualte effort (kilometers per year)
TransLen$site <- dat$site
TransLen$region <- dat$verbatimLocality
TransLen$locality <- dat$locality

effort <- tibble(TransLen) %>%
  select(!c(locationID,site)) %>%
  group_by(locality) %>%
  summarise(across(where(is.numeric), ~sum(.x)/1000, .names = "Sum_{.col}")) %>%
  mutate()

print(effort, n=1000)

# check yearly total
colSums(effort[,-1])

# write table with effort
#write.table(effort, file="data/effort_2024.txt", sep="\t", row.names=F, quote=F)
# list of sites in each cluster

# Make list of sites in different cluster
result <- dat %>%
  group_by(clust) %>%
  summarise(values_in_a = list(site), .groups = "drop")

print(result$values_in_a)

# number of individuals seen per cluster per year
N_line_year_2 <- as_tibble(N_line_year)
dim(N_line_year_2)

N_line_year_2 

N_line_year_2$site <- dat$site
N_line_year_2$clust <- dat$clust
N_line_year_2$locality <- dat$locality

colnames(N_line_year_2)[1:25] <- as.character(2000:2024)

# number of birds seen in each cluster
ind_cluster <- tibble(N_line_year_2) %>%
  select(!c(site,clust)) %>%
  group_by(locality) %>%
  summarise(across(where(is.numeric), sum, .names = "Sum_{.col}")) %>%
  mutate()

print(ind_cluster, n=1000)

# number of observations per cluster
# Summarize count of rows for each combination of 'Year' and 'locality'
summary_counts <- d_obs %>%
  group_by(Year, locality) %>%
  summarize(row_count = n(), .groups = 'drop')

# Pivot the summary data to get 'Year' as columns and 'locality' as rows
wide_summary <- summary_counts %>%
  pivot_wider(names_from = Year, values_from = row_count, values_fill = list(row_count = 0))

# View the result
print(wide_summary)

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
  T = N_yearsTot, # Max number of years with data
  nsites = N_sites, # Total number of monitored sites per area
  area = A, # Total covered area per year per site
  dclass = dclass, # distance class for all observations
  site = temp_dist$site, # site for each observation
  Reg = dat$reg,
  Clust = dat$clust,
  NClust = length(unique(dat$clust)),
  VerbClust = dat$locality
)

# save inputdata    
save(input.data, file = "data/RypeData_GBIF00-24_v2.rds")
#- end of script
