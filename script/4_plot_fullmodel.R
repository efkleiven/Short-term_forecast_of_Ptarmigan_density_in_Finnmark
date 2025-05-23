# script to plot 

# read libraries
library(tidyverse)
library(jagsUI)

# load jagsoutput
#....
getwd()

# full model with all years
load("./output/earlymodel_00-24.RData")

#---------------------------
# plot overall mean
#---------------------------
# prepare data for ggplot
mean <- apply(out24$Dtot2, 2, mean)
low <- apply(out24$Dtot2, 2, function(x) quantile(x, probs = 0.025))
high <- apply(out24$Dtot2, 2, function(x) quantile(x, probs = 0.975))
year <- 2000:2024

# store in tibble
dat <- tibble(mean=mean, low=low, high=high, year=year)
print(dat, n=25)

# make ggplot
ggplot(dat, aes(x = year))  + 
  geom_line(aes(y = mean), color = "#C2B391")+ 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, fill = "#C2B391") + 
  scale_x_continuous(
    breaks = c(2000:2024),
    limits = c(2000, 2024)) + 
  scale_y_continuous(
    breaks = c(0,5,10,15,20,25,30,35,40),
    limits = c(0,42)) + 
  ylab(bquote("Average population density " (birds/km^2))) +  
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.75))


ggsave("plot/Pt_Finnmark_24.png",
       width=30, height=18, units="cm")
# for mean mu

#---------------------------
# plot overall mean
#---------------------------
# prepare data for ggplot
mean <- apply(out24$mean_mu, 2, mean)
low <- apply(out24$mean_mu, 2, function(x) quantile(x, probs = 0.025))
high <- apply(out24$mean_mu, 2, function(x) quantile(x, probs = 0.975))
year <- 2000:2024

# store in tibble
dat <- tibble(mean=mean, low=low, high=high, year=year)
print(dat, n=25)

# make ggplot
ggplot(dat, aes(x = year))  + 
  geom_line(aes(y = mean), color = "#C2B391")+ 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, fill = "#C2B391") + 
  scale_x_continuous(
    breaks = c(2000:2024),
    limits = c(2000, 2024)) + 
  scale_y_continuous(
    breaks = c(0:3),
    limits = c(0,3)) + 
  ylab(bquote("Average population density " (birds/km^2))) +  
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.75))

#-----------------------------------------------
# For MU
# prepare data for ggplot
mean <- colMeans(apply(exp(out24$mu), c(2,3), mean))
low <- apply(apply(exp(out24$mu), c(2,3), mean), 2, quantile, probs=0.025)
high <- apply(apply(exp(out24$mu),c(2,3), mean), 2, quantile, probs=0.975)
year <- 2000:2024

colMeans(apply(exp(out24$mu), c(2,3), mean))

# store in tibble
dat <- tibble(mean=mean, low=low, high=high, year=year)

# make ggplot
ggplot(dat, aes(x = year))  + 
  geom_line(aes(y = mean), color = "#C2B391")+ 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, fill = "#C2B391") + 
  scale_x_continuous(
    breaks = c(2000:2024),
    limits = c(2000, 2024)) + 
  scale_y_continuous(
    breaks = c(0,5,10,15,20,25,30,35,40,45),
    limits = c(0,35)) + 
  ylab(bquote("Average population density " (birds/km^2))) +  
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.75))
#----------------------------------------------------------------------------------------------------------

# Plot region mean
load("data/data_JAGS_2024.rds")
load("d_chicks.rds")

# yearly mean for each region
mean1 <- apply(out24$Dreg_1, 2, mean)
mean2 <- apply(out24$Dreg_2, 2, mean)
mean3 <- apply(out24$Dreg_3, 2, mean)
mean4 <- apply(out24$Dreg_4, 2, mean)

# yearly lower 95% CI
low1 <- apply(out24$Dreg_1, 2, function(x) quantile(x, probs = 0.025))
low2 <- apply(out24$Dreg_2, 2, function(x) quantile(x, probs = 0.025))
low3 <- apply(out24$Dreg_3, 2, function(x) quantile(x, probs = 0.025))
low4 <- apply(out24$Dreg_4, 2, function(x) quantile(x, probs = 0.025))

# yearly upper 95% CI
high1 <- apply(out24$Dreg_1, 2, function(x) quantile(x, probs = 0.975))
high2 <- apply(out24$Dreg_2, 2, function(x) quantile(x, probs = 0.975))
high3 <- apply(out24$Dreg_3, 2, function(x) quantile(x, probs = 0.975))
high4 <- apply(out24$Dreg_4, 2, function(x) quantile(x, probs = 0.975))

# define year vector for plot
year <- 2000:2024

# store in tibble
pt_reg <- tibble(mean=c(mean1,mean2,mean3,mean4),
                 low=c(low1, low2, low3, low4), high=c(high1, high2, high3, high4), 
                 Year=rep(year, times=4), 
                 reg=c(rep("Indre", times=25), rep("Vest", times=25), rep("Pasvik", times=25), rep("Øst", times=25)))

# save pt_reg
#save(pt_reg, file="data/pt_reg.rds")

# explains why ptarmigan numbers did not decline last year
#(mean1+mean2+mean4)/3

# make regional plot 
pt_reg %>%
  filter(!reg=="Pasvik") %>%
  ggplot( aes(x = Year))  + 
  geom_line(aes(y = mean, color = reg), linewidth=1.5)+
  geom_ribbon(aes(ymin = low, ymax = high, fill=reg), alpha = 0.5) + # so far uncertainty too large for this to make sense
  scale_x_continuous(
    breaks = c(2000:2024),
    limits = c(1999.5, 2024.5)) + 
  scale_colour_manual(values = c("Indre" = "#999999", "Vest" = "#56B4E9", "ost"="#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_y_continuous(
    breaks = c(seq(0,70, by=5)),
    limits = c(0,51)) + 
  ylab(bquote("Average population density" (birds/km^2))) +  
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.75))

ggsave("plot/Pt_Reg_Finnmark_2024_.png",
       width=30, height=18, units="cm")


# make plot including chick production
d_chicks$verbatimLocality[d_chicks$verbatimLocality=="Indre Finnmark"] <- "Indre"
d_chicks$verbatimLocality[d_chicks$verbatimLocality=="Vest Finnmark kyst"] <- "Vest"
d_chicks$verbatimLocality[d_chicks$verbatimLocality=="�st Finnmark"] <- "�st"

names(d_chicks)[1] <- "reg"

dat_reg2 <- left_join(pt_reg, d_chicks)

# make plot including proportion of chicks
dat_reg2 %>%
  filter(!reg=="Pasvik") %>%
  ggplot( aes(x = Year))  + 
  geom_line(aes(y = mean, color = reg), linewidth=1.5)+
  # geom_ribbon(aes(ymin = low, ymax = high, fill=reg), alpha = 0.5) + # so far uncertainty too large for this to make sense
  geom_bar(aes(y=chicks2, fill=reg), position = position_dodge(), stat = 'identity' ) +
  scale_x_continuous(
    breaks = c(2000:2024),
    limits = c(1999.5, 2024.5)) + 
  scale_colour_manual(values = c("Indre" = "#999999", "Vest" = "#56B4E9", "Øst"="#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_y_continuous(
    breaks = c(seq(0,100, by=5)),
    limits = c(0,70)) + 
  ylab(bquote("Chicks per hen & Average population density" (birds/km^2))) +  
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.75))

#ggsave("plot/Pt_Reg_Finnmark_clustermod_chicks.png",
#       width=30, height=18, units="cm")

#------------------------------------------------------------------------------------------------
dat_boxplot <- tibble(value=c(out24$btDD-1, out24$btR, out24$btRD,  out24$btharvt, out24$btDoYft, out24$bttrend), #, out24$btMOTHt, out24$btcarc), 
                      par=c(rep('DD', times=18000), rep('Rod t', times=18000), rep('Rod t-1', times=18000),  rep('Harvest', times=18000), rep('OnsetFall', times=18000), rep('trend', times=18000) )) #,, rep('Moth', times=18000), rep('Carcas', times=18000) ))

ggplot(dat_boxplot, aes(y=value, x=par))+
  geom_violin()+
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.75),
        axis.title = element_blank())

ggsave("plot/Pt_Finnmark_Cov_24_earlymod.png",
       width=30, height=18, units="cm")

#--------------------------------
# look at spatial and temporal components of onset of fall effects
dat_boxplot <- tibble(value = c(out24$btDoYf, out24$btDoYft, out24$btDoYfs), 
                      par=c(rep('Residual', times=18000),rep('temporal', times=18000), rep('spatial', times=18000) ))

ggplot(dat_boxplot, aes(y=value, x=par))+
  geom_violin()+
  ggtitle("Onset of fall") +                    # Add main title
  theme(axis.title.x = element_blank(),         # Remove x-axis title
        axis.title.y = element_blank(),         # Remove y-axis title
        axis.text.x = element_text(size = 14),  # Increase x-axis text size
        axis.text.y = element_text(size = 14),  # Increase y-axis text size
        plot.title = element_text(size = 20))+  # Increase title size
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1)  # Add horizontal dashed line at y = 0 

ggsave("plot/Cov_OnsetOfFall_24.png",
       width=30, height=18, units="cm")

#--------------------------------
# look at spatial and temporal components of harvest effects
dat_boxplot <- tibble(value = c(out24$btharv, out24$btharvt, out24$btharvs), 
                      par=c(rep('Residual', times=18000),rep('temporal', times=18000), rep('spatial', times=18000) ))

ggplot(dat_boxplot, aes(y=value, x=par))+
  geom_violin()+
  ggtitle("Harvest") +                    # Add main title
  theme(axis.title.x = element_blank(),         # Remove x-axis title
        axis.title.y = element_blank(),         # Remove y-axis title
        axis.text.x = element_text(size = 14),  # Increase x-axis text size
        axis.text.y = element_text(size = 14),  # Increase y-axis text size
        plot.title = element_text(size = 20))+  # Increase title size
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1)  # Add horizontal dashed line at y = 0 

#- end of script