## Cod condition - evaluate results 
## Use best model to predict condition for each year / bay combination (for 2018-2020)
## This was used as a diagnostic step to assess modeled results


library(plyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
library(tidyverse)
source("./scripts/stan_utils.R")


## Read in data --------------------------------------------
cod.condition.data <- read.csv("data/condition_data.csv", row.names = 1)

# remove station 117, which does not have cpue data
drop <- cod.condition.data$station==117
cod.condition.data <- cod.condition.data[!drop,]

# and remove Kujulik!
drop <- cod.condition.data$bay=="Kujulik"
cod.condition.data <- cod.condition.data[!drop,]

# check Cooks!
cod.condition.data$bay <- as.character(cod.condition.data$bay)
unique(cod.condition.data$bay)
change <- cod.condition.data$bay=="Cooks"
cod.condition.data$bay[change] <- "Cook Bay"


cod.condition.data$bay_fac <- as.factor(cod.condition.data$bay)
cod.condition.data$year_fac <- as.factor(cod.condition.data$year)
cod.condition.data$site_fac <- as.factor(cod.condition.data$site)
cod.condition.data$bay_site_fac <- as.factor(paste0(cod.condition.data$bay, "_", cod.condition.data$site))
cod.condition.data$date <- as.Date(cod.condition.data$julian,
                                   origin = paste0(cod.condition.data$year, "-01-01"))
cod.condition.data$ln.cpue <- log(cod.condition.data$cod.cpue)

# predict summer temp for each bay/year!

# read in temperature model
seine.temp <- readRDS("./output/temp3_gauss.rds") 
summary(seine.temp)

# create new data frame to predict from best model for Kdry variability

# get mean length, mean ln.cpue, and mean julian day for all sampling events
bay_fac <- unique(cod.condition.data$bay_fac)

new.dat <- data.frame(year_fac=as.factor(rep(2018:2020, each=length(bay_fac))),
                      bay_fac=bay_fac,
                      ln.cpue=mean(cod.condition.data$ln.cpue),
                      length=mean(cod.condition.data$length),
                      julian=mean(cod.condition.data$julian))

fitted.temp <- as.data.frame(fitted(seine.temp, newdata = new.dat))

new.dat$fit.temp.mu <- fitted.temp$Estimate
new.dat$LCI <- fitted.temp$Q2.5
new.dat$UCI <- fitted.temp$Q97.5


# best model for Kdry was cod_cond2
cod_cond2 <- readRDS("./output/cod_cond2.rds")

check_hmc_diagnostics(cod_cond2$fit)
neff_lowest(cod_cond2$fit)
rhat_highest(cod_cond2$fit)
summary(cod_cond2)
bayes_R2(cod_cond2)

predict <- as.data.frame(fitted(cod_cond2, newdata = new.dat))

plot.dat <- cbind(new.dat, predict)

# align bays by longitude
bay.dat <- read.csv("./data/bay_lat_long.csv", row.names = 1)
names(bay.dat)[1] <- "bay_fac"
plot.dat <- left_join(plot.dat, bay.dat)
plot.dat$bay_fac <- reorder(plot.dat$bay_fac, desc(plot.dat$lon))

cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_bw())

ggplot(plot.dat, aes(bay_fac, fit.temp.mu, fill=year_fac)) +
  geom_bar(position = "dodge", stat="identity", color="grey") +
  geom_errorbar(aes(ymin=LCI, ymax=UCI), position="dodge", color="dark grey") +
  scale_fill_manual(values=cb[c(2,4,6)]) +
  scale_color_manual(values=cb[c(2,4,6)]) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=45, vjust = 1, hjust=1))


ggplot(plot.dat, aes(bay_fac, Estimate, fill=year_fac)) +
  geom_bar(position = "dodge", stat="identity", color="grey") +
  geom_errorbar(aes(ymin=Q2.5, ymax=Q97.5), position="dodge", color="dark grey") +
  scale_fill_manual(values=cb[c(2,4,6)]) +
  scale_color_manual(values=cb[c(2,4,6)]) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=45, vjust = 1, hjust=1))

avg.dat <- cod.condition.data %>%
  group_by(year_fac, bay_fac) %>%
  summarise(mean=mean(Kdry),
            sd=sd(Kdry),
            n=n())

avg.dat$LCI <- avg.dat$mean - 1.96*avg.dat$sd/sqrt(avg.dat$n)
avg.dat$UCI <- avg.dat$mean + 1.96*avg.dat$sd/sqrt(avg.dat$n)

avg.dat <- left_join(avg.dat, bay.dat)
avg.dat$bay_fac <- reorder(avg.dat$bay_fac, desc(avg.dat$lon))

ggplot(avg.dat, aes(bay_fac, mean, fill=year_fac)) +
  geom_bar(position = "dodge", stat="identity", color="grey") +
  geom_errorbar(aes(ymin=LCI, ymax=UCI), position="dodge", color="dark grey") +
  scale_fill_manual(values=cb[c(2,4,6)]) +
  scale_color_manual(values=cb[c(2,4,6)]) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=45, vjust = 1, hjust=1))
