## Script to run the model from Henden et al 2020
# the model is simplified; it contains no regional effects or covariates
# the data used is downloaded form GBIF and prepared using LivingNorway-r-package

library(tidyverse)
#rm(list = ls())

setwd("/data/P-Prosjekter2/182411_20_rypemodul_finnmark")

cat(file= "HDS_reg&clust_DD&rod2&DoYf&harv&rrtg&carc&moth&trend.txt", 
"model{

  # Priors
    # detecion par
  alpha0 ~ dunif(-5,5)      # intercept detection probability
  
    # initial year par
  btfR  ~ dunif(-5,5)      # rodent effect in initial year
    
    # dynamic model par
  btDD ~ dunif(-5,5)     # density dependence effect
  btR ~ dunif(-5,5)      # rodent effect from year t=2
  btRD ~ dunif(-5,5)     # delayed rodent effect from year t=2
  btcarc ~ dunif(-5,5)   # carcas effect from year t=2
  bttrend ~ dunif(-5,5)  # trend effect from year t=2
  
  btDoYf ~ dunif(-5,5)   # onset of fall residual effect
  btDoYft ~ dunif(-5,5)  # onset of fall temporal effect
  btDoYfs ~ dunif(-5,5)  # onset of fall spatial effect
  
  btMOTH ~ dunif(-5,5)    # moth residual effect
  btMOTHt ~ dunif(-5,5)   # moth temporal effect 
  btMOTHs ~ dunif(-5,5)   # moth spatial effect 
  
  btharv ~ dunif(-5,5)   # havest residual effect
  btharvt ~ dunif(-5,5)  # havest temporal effect
  btharvs ~ dunif(-5,5)  # harvest spatial effect
  
  btrr ~ dunif(-5,5)   # havest residual effect
  btrrt ~ dunif(-5,5)  # havest temporal effect
  btrrs ~ dunif(-5,5)  # harvest spatial effect
  
  bttg ~ dunif(-5,5)   # havest residual effect
  bttgt ~ dunif(-5,5)  # havest temporal effect
  bttgs ~ dunif(-5,5)  # harvest spatial effect
  
    #regional fixed effects
  for(r in 1:4){
    bt0[r] ~ dunif(-5,5)  # regional fixed effect (T>1)
    btf0[r] ~ dunif(-5,5)   # regional fixed effect (T=1)
  }

  ## Precision (T=1) and process variation (T>1)  
  PrEc <- 1 / (sdprectau * sdprectau) 
  sdprectau ~ dunif(0, 1) 

  PrOc <- 1 /(sdproctau * sdproctau) 
  sdproctau ~ dunif(0,1)  
  
  ## Random transect cluster effect 
  for (j in 1:NClust){ # number of clusters of lines
    rCl1[j]~ dnorm(0,rtau)    # initial year
    rCl[j]~ dnorm(0,rtau2)   # consequtive year
    }# end cluster loop
    
  rtau <- 1 / (sdrtau * sdrtau)     # precision
  rtau2 <- 1 /(sdrtau2 * sdrtau2)   # precision
  sdrtau ~ dunif(0,1)
  sdrtau2 ~ dunif(0,1)

  ## Distance sampling observation model for observed (binned) 
  # Distance data:
  for(i in 1:nobs){ 
    dclass[i] ~ dcat(fc[site[i],]+ 0.01)                       
  } 
  
  ## Likelihood
  for(s in 1:nsites){  # sites/lines 
  
    # Construct cell probabilities for nD multinomial cells
    for(g in 1:nD){    # distance classes
      # Half-normal detection function     
      log(p[s,g]) <- -midpt[g] * midpt[g] / (2 *sigma[s]*sigma[s]) 
      pi[s,g] <- delta / B          
      f[s,g] <- p[s,g] * pi[s,g]    
      fc[s,g] <- f[s,g] / pcap[s]   
    } 
    pcap[s] <- sum(f[s,]) # Pr(capture): sum of rectangular areas
    
    ## Process model:  
    # Initial density model for Year 1 (cf. Sillett et al. 2012):
    y[s,1] ~ dbin(pcap[s], N[s,1])                          
    N[s,1] ~ dpois( lambda[s,1] )                                  	 
    lambda[s,1] <-  D[s,1] * area[s,1]
    logD[s,1] ~ dnorm(mu[s,1], PrEc)                                         
    
    mu[s,1] <- btf0[Reg[s]] + rCl1[Clust[s]] + btfR * rod[1, Reg[s]]
    D[s,1] <- exp(logD[s,1])      
    
    log(sigma[s]) <- alpha0     
    
    # Population dynamics model for subsequent years: 
    for (t in 2:T){
      y[s,t] ~ dbin(pcap[s], N[s,t])        
      N[s,t] ~ dpois(lambda[s,t])  
      lambda[s,t] <- D[s,t] * area[s,t]
      logD[s,t] ~ dnorm( mu[s,t] , PrOc )
      
      mu[s,t] <- bt0[Reg[s]] + rCl[Clust[s]]+ btDD * mu[s,t-1] + btR * rod[t, Reg[s]] + btRD * rod[t-1, Reg[s]] + 
      btDoYf * DoYf[Reg[s], t-1] + btDoYft * DoYft[t-1] + btDoYfs * DoYfs[Reg[s]] +
      btMOTH * moth[Reg[s], t] + btMOTHt * motht[t] + btMOTHs * moths[Reg[s]] +
      btharv * harv[s, t-1] + btharvt * harvt[t-1] + btharvs * harvs[s] +
      btrr * rr[Reg[s], t] + btrrt * rrt[t] + btrrs * rrs[Reg[s]] + 
      bttg * tg[Reg[s], t] + bttgt * tgt[t] + bttgs * tgs[Reg[s]] +
      btcarc * carc[t] + bttrend * (t-1)
      
      D[s,t] <- exp(logD[s,t])
      
    }  # End pop dynamics subsequent years
  } # End likelihood
  
  ## Derived quantities
  # Total population size and density:
  
  for (t in 1:T){
    Ntotal[t] <- sum(N[,t])
    areaTot[t] <- sum(area[,t])
    Dtot[t] <- Ntotal[t]/areaTot[t]
  } # end time loop
  
}" #end model
)

#read data
#load("BugsData_Gbif_2022.rda")
load("data/RypeData_GBIF00-23_v3.rds")

# load pre-processed data
BugsData<-input.data

# year 2000 to 2023
BugsData$T <- 24

# load rodent data
# Øst
rodOst <- tibble(read.delim("data/storskala_04-23spring.txt")) %>% 
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

#Indre
rodInnerVest <- tibble(read.table("data/porsanger_mus_reg_2023.csv", sep=",", header=T)) %>%
  filter(seas=="SPRING")%>%
  dplyr::select(c(year,  rufoCold, rufoSold))

# make year a character
rodInnerVest$year <- as.character(rodInnerVest$year) 

# Vest  
rod1 <- rodOst3 %>% 
        right_join( rodInnerVest,  by=join_by(year))%>%
        left_join( rodOst3,   by=join_by(year))%>%
          dplyr::select(!year)%>% as.matrix()

# scale and make array
rod3 <- scale(array(rod1, dim=c(24,4)))

BugsData$rod <- rod3

# scale DoYf
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

# Inits
SimT = BugsData$T

#Nst <- BugsData$y[,1:SimT] + 5 # for no of individuals
Nst <- BugsData$y[,1:SimT] # for no of individuals

inits <- function(){list(N=Nst,sdproctau=runif(1,0,1),sdprectau=runif(1,0,1),alpha0=0, 
                         bt0=array(0, dim=c(4)), btf0=rep(0,times=4), btDD=0, 
                         btfR=0, btR=0, btRD=0, btDoYf=0, btDoYft=0, btDoYfs=0,
                         btharv=0, btharvt=0, btharvs=0, btcar=0, bttrend=0)}      

# Params to save
params <- c("sdproctau","sdprectau","PrOc","PrEc", "alpha0","btf0", "bt0", "bt1", "rCl",
            "btDD", "btR", "btRD", "btDoYf", "btDoYft","btDoYfs", "btMOTH", "btMOTHt", "btMOTHs", 
            "btharv", "btharvt","btharvs","btrr", "btrrt","btrrs","bttg", "bttgt","bttgs","btcarc","bttrend",
            "Ntotal", "mu", "D", "Dtot")   

ni <- 40000;   nb <- 20000;   nt <- 10;   nc <- 6;  na=5000 #

library(jagsUI)

out23 <- jags(BugsData, inits, params, "HDS_reg&clust_DD&rod2&DoYf&harv&rrtg&carc&moth&trend.txt", n.thin=nt,
                 n.chains=nc, n.burnin=nb, n.adapt=na, n.iter=ni, parallel=T)   # 

setwd("./output")
save(out23, file="HSD_reg&cluster_dd&rod2&DoYf&harv&rrtg&carc&anomaly&trend_00-23_v2.RData")
#- end of script
