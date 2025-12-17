## Script to run the model from Henden et al 2020
# the model is simplified; it contains no regional effects or covariates
# the data used is downloaded form GBIF and prepared using LivingNorway-r-package

library(tidyverse)
#rm(list = ls())

setwd("/data/P-Prosjekter2/182411_20_rypemodul_finnmark")

cat(file= "earlymodel.txt", 
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
  bttrend ~ dunif(-10,10)  # trend effect from year t=2
  btcarc ~ dunif(-10,10)
  
  btDoYf ~ dunif(-10,10)   # onset of fall residual effect
  btDoYft ~ dunif(-10,10)  # onset of fall temporal effect
  btDoYfs ~ dunif(-10,10)  # onset of fall spatial effect
  
  btMOTH ~ dunif(-10,10)    # moth residual effect
  btMOTHt ~ dunif(-10,10)   # moth temporal effect 
  btMOTHs ~ dunif(-10,10)   # moth spatial effect 
  
  btharv ~ dunif(-10,10)   # havest residual effect
  btharvt ~ dunif(-10,10)  # havest temporal effect
  btharvs ~ dunif(-10,10)  # harvest spatial effect
  
  btrr ~ dunif(-10,10)   # havest residual effect
  btrrt ~ dunif(-10,10)  # havest temporal effect
  btrrs ~ dunif(-10,10)  # harvest spatial effect
  
  bttg ~ dunif(-10,10)   # havest residual effect
  bttgt ~ dunif(-10,10)  # havest temporal effect
  bttgs ~ dunif(-10,10)  # harvest spatial effect

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
      
      mu[s,t] <- bt0[Reg[s]] + rCl[Clust[s]]+ btDD * mu[s,t-1] + btR * rod[t, Reg[s]] + btRD * rod[t-1, Reg[s]] + 
                 btDoYf * DoYf[Clust[s], t-1] + btDoYft * DoYft[t-1] + btDoYfs * DoYfs[Clust[s]] +
                 btMOTH * moth[Clust[s], t-1] + btMOTHt * motht[t-1] + btMOTHs * moths[Clust[s]] +
                 btharv * harv[s, t-1] + btharvt * harvt[t-1] + btharvs * harvs[s] +
                 btrr * rr[Clust[s], t] + btrrt * rrt[t] + btrrs * rrs[Clust[s]] +
                 bttg * tg[Clust[s], t] + bttgt * tgt[t] + bttgs * tgs[Clust[s]] +
                 btcarc * carc[t] + bttrend * (t-1)
                 
      logD[s,t] ~ dnorm( mu[s,t] , PrOc )
      D[s,t] <- exp(logD[s,t])
      
    }  # End pop dynamics subsequent years
  } # End likelihood
  
  # Derived parameters
  
  # Mean cluster density
    for (t in 1:T){
    Dclust_1[t] <- mean(D[clust_id_1, t])
    Dclust_2[t] <- mean(D[clust_id_2, t])
    Dclust_3[t] <- mean(D[clust_id_3, t])
    Dclust_4[t] <- mean(D[clust_id_4, t])
    Dclust_5[t] <- mean(D[clust_id_5, t])
    Dclust_6[t] <- mean(D[clust_id_6, t])
    Dclust_7[t] <- mean(D[clust_id_7, t])
    Dclust_8[t] <- mean(D[clust_id_8, t])
    Dclust_9[t] <- mean(D[clust_id_9, t])
    Dclust_10[t] <- mean(D[clust_id_10, t])

    Dclust_11[t] <- mean(D[clust_id_11, t])
    Dclust_12[t] <- mean(D[clust_id_12, t])
    Dclust_13[t] <- mean(D[clust_id_13, t])
    Dclust_14[t] <- mean(D[clust_id_14, t])
    Dclust_15[t] <- mean(D[clust_id_15, t])
    Dclust_16[t] <- mean(D[clust_id_16, t])
    Dclust_17[t] <- mean(D[clust_id_17, t])
    Dclust_18[t] <- mean(D[clust_id_18, t])
    Dclust_19[t] <- mean(D[clust_id_19, t])
    Dclust_20[t] <- mean(D[clust_id_20, t])
    
    Dclust_21[t] <- mean(D[clust_id_21, t])
    Dclust_22[t] <- mean(D[clust_id_22, t])
    Dclust_23[t] <- mean(D[clust_id_23, t])
    Dclust_24[t] <- mean(D[clust_id_24, t])
    Dclust_25[t] <- mean(D[clust_id_25, t])
    Dclust_26[t] <- mean(D[clust_id_26, t])
    Dclust_27[t] <- mean(D[clust_id_27, t])

    
# mean region density
Dreg_1[t] <- mean(c(Dclust_1[t], Dclust_5[t], Dclust_8[t], Dclust_9[t], Dclust_10[t], Dclust_13[t], Dclust_15[t], Dclust_21[t], Dclust_22[t], Dclust_24[t], Dclust_25[t])) # indre finnmark
Dreg_2[t] <- mean(c(Dclust_2[t], Dclust_4[t], Dclust_7[t], Dclust_12[t], Dclust_23[t]))                                            # vest finnmark
Dreg_3[t] <- mean(c(Dclust_3[t]))                                                                                                                                           # pasvik
Dreg_4[t] <- mean(c(Dclust_6[t], Dclust_11[t], Dclust_14[t], Dclust_16[t], Dclust_17[t], Dclust_18[t], Dclust_19[t], Dclust_20[t], Dclust_26[t], Dclust_27[t]))             # ?st finnmark

# total density
Dtot[t] <- mean(c(Dreg_1[t],Dreg_2[t],Dreg_3[t],Dreg_4[t])) # mean density for all 4 regions  
Dtot2[t] <- mean(c(Dreg_1[t],Dreg_2[t],Dreg_4[t])) # mean density without pasvik

#mean mu
mean_mu[t] <- mean(mu[,t])
    } # end time loop
    
 }" #end model
)

#read data
load("Short-term_forecast_of_Ptarmigan_density_in_Finnmark/data/data_JAGS_2024_v4.rds")

# load pre-processed data
BugsData<-input.data

max(BugsData$site)
str(BugsData)

# remove objects not needed
BugsData$VerbClust <- NULL

# vector giving which lines belong to which cluster
BugsData$clust_id_1 <- which(input.data$Clust==1)
BugsData$clust_id_2 <- which(input.data$Clust==2)
BugsData$clust_id_3 <- which(input.data$Clust==3)
BugsData$clust_id_4 <- which(input.data$Clust==4)
BugsData$clust_id_5 <- which(input.data$Clust==5)
BugsData$clust_id_6 <- which(input.data$Clust==6)
BugsData$clust_id_7 <- which(input.data$Clust==7)
BugsData$clust_id_8 <- which(input.data$Clust==8)
BugsData$clust_id_9 <- which(input.data$Clust==9)
BugsData$clust_id_10 <- which(input.data$Clust==10)

BugsData$clust_id_11 <- which(input.data$Clust==11)
BugsData$clust_id_12 <- which(input.data$Clust==12)
BugsData$clust_id_13 <- which(input.data$Clust==13)
BugsData$clust_id_14 <- which(input.data$Clust==14)
BugsData$clust_id_15 <- which(input.data$Clust==15)
BugsData$clust_id_16 <- which(input.data$Clust==16)
BugsData$clust_id_17 <- which(input.data$Clust==17)
BugsData$clust_id_18 <- which(input.data$Clust==18)
BugsData$clust_id_19 <- which(input.data$Clust==19)
BugsData$clust_id_20 <- which(input.data$Clust==20)

BugsData$clust_id_21 <- which(input.data$Clust==21)
BugsData$clust_id_22 <- which(input.data$Clust==22)
BugsData$clust_id_23 <- which(input.data$Clust==23)
BugsData$clust_id_24 <- which(input.data$Clust==24)
BugsData$clust_id_25 <- which(input.data$Clust==25)
BugsData$clust_id_26 <- which(input.data$Clust==26)
BugsData$clust_id_27 <- which(input.data$Clust==27)


# which reg does the clusters correspond too
# Get unique Clust values
unique_clust <- unique(BugsData$Clust)

# Match unique Clust values with corresponding Reg values
matched_reg <- sapply(unique_clust, function(clust_value) {
  unique(BugsData$Reg[BugsData$Clust == clust_value])
})

# Print result
df <- data.frame(Clust = unique_clust, Reg = matched_reg)

which(df$Reg==1) # inner
which(df$Reg==2) # west
which(df$Reg==3) # pasvik
which(df$Reg==4) # east

# reorganize rodent data to match regions
BugsData$rod <- BugsData$rod[, c(3, 2, 1, 4)]

### Standardize predictor variables
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

# define number of years in the dataset
BugsData$T <- dim(BugsData$y)[2]-4  # to get data until 2020

# Spesify initial values for MCMC
SimT = BugsData$T 

Nst <- BugsData$y[,1:SimT] # for no of individuals

# list all initial values
inits <- function(){list(N=Nst,sdproctau=runif(1,0,1),sdprectau=runif(1,0,1),alpha0=0, 
                         bt0=array(0, dim=c(4)), btf0=rep(0,times=4), btDD=0, 
                         btfR=0, btR=0, btRD=0, btDoYf=0, btDoYft=0, btDoYfs=0,
                         btMOTH=0, btMOTHt=0, btMOTHs=0,
                         btcarc=0,
                         btrr=0, btrrt=0, btrrs=0,
                         bttg=0, bttgt=0, bttgs=0,
                         btharv=0, btharvt=0, btharvs=0, bttrend=0)}      

# Params to save
params <- c( "btDD", "btR", "btRD", "btDoYf", "btDoYft","btDoYfs", "btharv", "btharvt","btharvs","bttrend",
             "btMOTH", "btMOTHt", "btMOTHs",
             "btcarc","btrr", "btrrt","btrrs","bttg", "bttgt","bttgs",
            "Ntotal", "mu", "D", "Dtot", "mean_mu", "bt0", "rCl", 
            "Dclust_1", "Dclust_2", "Dclust_3","Dclust_4","Dclust_5","Dclust_6","Dclust_7","Dclust_8","Dclust_9","Dclust_10",
            "Dclust_11", "Dclust_12", "Dclust_13","Dclust_14","Dclust_15","Dclust_16","Dclust_17","Dclust_18","Dclust_19","Dclust_20",
            "Dclust_21", "Dclust_22", "Dclust_23","Dclust_24","Dclust_25","Dclust_26","Dclust_27","Dclust_28","Dclust_29","Dclust_30",
            "Dreg_1","Dreg_2","Dreg_3","Dreg_4", "Dtot", "Dtot2")   

#MCMC settings
ni <- 60000;   nb <- 30000;   nt <- 10;   nc <- 6;  na=10000 #
#ni <- 60;   nb <- 30;   nt <- 1;   nc <- 1 #

library(jagsUI)

# run model in JAGS
out_full_20 <- jags(BugsData, inits, params, "earlymodel.txt", n.thin=nt,
                 n.chains=nc, n.burnin=nb, n.adapt=na, n.iter=ni, parallel=T)   # 

# keep only MCMC values
out_full_20 <- out_full_20$sims.list

# save output
save(out_full_20, file="output/fullmod_00-20.RData")

#- the end of script