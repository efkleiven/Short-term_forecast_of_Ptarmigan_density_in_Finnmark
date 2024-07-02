## Script to run the model from Henden et al 2020
# the model is simplified; it contains no regional effects or covariates
# the data used is downloaded form GBIF and prepared using LivingNorway-r-package

library(tidyverse)
#rm(list = ls())

setwd("/data/P-Prosjekter2/182411_20_rypemodul_finnmark")

cat(file= "foodwebmodel.txt", 
"model{

  # Priors
    # detecion par
  alpha0 ~ dunif(-10,10)      # intercept detection probability
  
    # initial year par
  btfR  ~ dunif(-10,10)      # rodent effect in initial year
    
    # dynamic model par
  btDD ~ dunif(-10,10)     # density dependence effect
  btR ~ dunif(-10,10)      # rodent effect from year t=2
  btRD ~ dunif(-10,10)     # delayed rodent effect from year t=2
  btcarc ~ dunif(-10,10)   # carcas effect from year t=2
  bttrend ~ dunif(-10,10)  # trend effect from year t=2
  
  #btDoYf ~ dunif(-10,10)   # onset of fall residual effect
  #btDoYft ~ dunif(-10,10)  # onset of fall temporal effect
  #btDoYfs ~ dunif(-10,10)  # onset of fall spatial effect
  
  btMOTH ~ dunif(-10,10)    # moth residual effect
  btMOTHt ~ dunif(-10,10)   # moth temporal effect 
  btMOTHs ~ dunif(-10,10)   # moth spatial effect 
  
  btharv ~ dunif(-10,10)   # havest residual effect
  btharvt ~ dunif(-10,10)  # havest temporal effect
  btharvs ~ dunif(-10,10)  # harvest spatial effect
  
    #regional fixed effects
  for(r in 1:4){
    bt0[r] ~ dunif(-10,10)  # regional fixed effect (T>1)
    btf0[r] ~ dunif(-10,10)   # regional fixed effect (T=1)
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
    #  btDoYf * DoYf[Clust[s], t-1] + btDoYft * DoYft[t-1] + btDoYfs * DoYfs[Clust[s]] +
      btMOTH * moth[Clust[s], t-1] + btMOTHt * motht[t-1] + btMOTHs * moths[Clust[s]] +
      btharv * harv[s, t-1] + btharvt * harvt[t-1] + btharvs * harvs[s] +
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
load("data/RypeData_GBIF00-23_pred24_v1.rds")

# load pre-processed data
BugsData<-input.data

# year 2000 to 2023
#BugsData$T <- 24

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
            "btDD", "btR", "btRD", "btMOTH", "btMOTHt", "btMOTHs", 
            "btharv", "btharvt","btharvs","btcarc","bttrend",
            "Ntotal", "mu", "D", "Dtot")   

ni <- 60000;   nb <- 40000;   nt <- 10;   nc <- 6;  na=10000 #

library(jagsUI)

out23 <- jags(BugsData, inits, params, "foodwebmodel.txt", n.thin=nt,
                 n.chains=nc, n.burnin=nb, n.adapt=na, n.iter=ni, parallel=T)   # 

setwd("./output")
save(out23, file="foodwebmodel_00-23.RData")
#- end of script
