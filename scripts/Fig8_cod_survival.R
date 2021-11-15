## Cod abundance between repeat visits to sites within years

library(ggplot2)
library(plyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
library(tidyverse)
source("./scripts/stan_utils.R")


## Read in data --------------------------------------------
cod.data <- read.csv("data/cpue.data.csv", row.names = 1)
cod.data$bay_fac <- as.factor(cod.data$bay)
cod.data$year_fac <- as.factor(cod.data$year)
cod.data$site_fac <- as.factor(cod.data$site)
cod.data$bay_site_fac <- as.factor(paste0(cod.data$bay, "_", cod.data$site))
cod.data$bay_year_fac <- as.factor(paste0(cod.data$bay, "_", cod.data$year))
cod.data$present <- ifelse(cod.data$cod > 0, 1, 0)
cod.data$date <- as.Date(cod.data$julian,
                         origin = paste0(cod.data$year, "-01-01"))

hist(cod.data$cod, breaks=50)

# read in temperature model
seine.temp <- readRDS("./output/temp3_gauss.rds") 

new.dat <- cod.data %>%
  select(bay_fac, year_fac)
new.dat$julian <- mean(cod.data$julian)
nrow(cod.data); nrow(new.dat)

fitted.temp <- as.data.frame(fitted(seine.temp, newdata = new.dat)) %>%
  select(Estimate)
nrow(fitted.temp)

new.dat <- cbind(new.dat, fitted.temp) %>%
  select(-julian)

## this includes multiple bay-year observations that create replicates in left_join!
## limit to unique bay_year combos!
new.dat$bay_year_fac <- as.factor(paste0(new.dat$bay, "_", new.dat$year))
new.dat <- new.dat %>%
  group_by(bay_year_fac) %>%
  summarise(fit.temp.mu = mean(Estimate))

nrow(cod.data); nrow(new.dat)

cod.data <- left_join(cod.data, new.dat)
nrow(cod.data)

# check fitted mu
plot(cod.data$fit.temp.mu)
hist(cod.data$fit.temp.mu, breaks = 20) 

cod.data$fourth.root.cpue <- cod.data$cod^0.25
hist(cod.data$fourth.root.cpue, breaks = 20)


# now restrict to sites sampled > once in a year!
ff <- function(x) length(unique(x))

check <- cod.data %>%
  group_by(site_fac, year_fac) %>%
  summarise(visits=ff(julian))
View(check)

check$site_year_fac <- as.factor(paste0(check$site_fac, "_", as.character(check$year_fac)))

cod.data$site_year_fac <- as.factor(paste0(cod.data$site, "_", cod.data$year))

cod.data <- left_join(cod.data, check)
nrow(cod.data)

surv.data <- cod.data %>%
  filter(visits > 1)
nrow(surv.data)  

## Check distributions
plot(surv.data$fourth.root.cpue)
hist(surv.data$fourth.root.cpue, breaks = 80) 

hist(surv.data$cod, breaks = 80) 

# get initial abundance for each site in each year!
first <- surv.data %>%
  group_by(site_year_fac) %>%
  summarise(julian=min(julian))

# separate out the abundance data for the first visit to each site / year combo
first.dat <- left_join(first, surv.data)
nrow(first.dat)

hist(first.dat$cod, breaks = 30)
hist(first.dat$julian, breaks = 30) 

# discard the first visits that were late in the season
first.dat <- first.dat %>%
  filter(julian < 210)

# get the 4th-root cpue at first visit
first.cpue <- first.dat %>%
  group_by(site_year_fac) %>%
  summarise(first.cpue.4=mean(fourth.root.cpue),
            julian.first=mean(julian))

hist(first.cpue$first.cpue.4, breaks = 30)

# and join it back in
surv.data <- left_join(surv.data, first.cpue)

# get the mean 4th-root cpue at all visits
mean.cpue <- surv.data %>%
  group_by(site_year_fac) %>%
  summarise(mean.cpue.4=mean(fourth.root.cpue))

hist(mean.cpue$mean.cpue.4, breaks = 30)

# and join it back in
surv.data <- left_join(surv.data, mean.cpue)


# and get the change in date and abundance from the first visit
surv.data$delta.date <- surv.data$julian - surv.data$julian.first
hist(surv.data$delta.date)

View(arrange(surv.data, delta.date))

# get change in abundance as proportion of original!
surv.data$delta.abund <- (surv.data$fourth.root.cpue - surv.data$first.cpue.4)/surv.data$first.cpue.4
hist(surv.data$delta.abund)

# use that to calculate survival rate (delta cpue.4/d)
surv.data$surv.rate <- surv.data$delta.abund / surv.data$delta.date

## drop the delta.date = 0 (first visits)
## AND drop delta date <= 5 (too short a duration to resolve signal vs noise!)
## this leaves only visits of 14 days or more separation
## and drop first.cpue.4 = 0 (no catch first visit)
## NB!! this is an important analysis decision, looking to 
## exclude sites where no cod were caught the first time out, with the idea that 
## they're poor estimates of initial abundance
surv.data <- surv.data %>%
  filter(delta.date > 5, first.cpue.4 > 0)
 
min(surv.data$first.cpue.4) 
hist(surv.data$surv.rate, breaks = 100)
## that looks great

# are julian first and delta.date colinear??
# yes! - use only julian.first!
ggplot(surv.data, aes(julian.first, delta.date)) +
  geom_point() +
  theme_bw()

## mgcv fits -----------------------------------------------
gam.surv.0 <- gam(surv.rate ~ s(julian.first, k = 4),
                    data=surv.data)

summary(gam.surv.0) 
## looks like controlling for Julian.first is important
## i.e., early first visits miss some of the initial settlement

plot(gam.surv.0, pages = 1)

gam.surv.1 <- gam(surv.rate ~ s(julian.first, k = 4) + s(fit.temp.mu, k = 4),
                    data=surv.data)

summary(gam.surv.1)
plot(gam.surv.1, pages = 1)


gam.surv.2 <- gam(surv.rate ~ s(julian.first, k = 4) + s(fit.temp.mu, k = 4) + s(mean.cpue.4, k=4),
                    data=surv.data)

summary(gam.surv.2)
plot(gam.surv.2, pages = 1)

gam.surv.3 <- gam(surv.rate ~  s(julian.first, k = 4) + s(fit.temp.mu, k = 4) + s(mean.cpue.4, k=4) + bay_fac,
                    data=surv.data)

summary(gam.surv.3)
plot(gam.surv.3, pages = 1)


MuMIn::AICc(gam.surv.0, gam.surv.1, gam.surv.2, gam.surv.3) # model 3 is best by far

## brms: setup ---------------------------------------------

surv.data$surv.rate_stnd <- as.vector(scale(surv.data$surv.rate))

## Define model formulas
surv0_formula <-  bf(surv.rate_stnd ~ s(julian.first, k = 4))

surv1_formula <-  bf(surv.rate_stnd ~ s(julian.first, k = 4) + s(fit.temp.mu, k = 4))

surv2_formula <-  bf(surv.rate_stnd ~ s(julian.first, k = 4) + s(fit.temp.mu, k = 4) + s(mean.cpue.4, k = 4))

surv3_formula <-  bf(surv.rate_stnd ~ s(julian.first, k = 4) + s(fit.temp.mu, k = 4) + s(mean.cpue.4, k = 4) + bay_fac)

## Show default priors
get_prior(surv0_formula, surv.data)

# using defaults for now!

## fit: brms --------------------------------------
cod_surv0 <- brm(surv0_formula,
                data = surv.data,
                cores = 4, chains = 4, iter = 4000,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99, max_treedepth = 10))
cod_surv0  <- add_criterion(cod_surv0, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_surv0, file = "output/cod_surv0.rds")

cod_surv0 <- readRDS("./output/cod_surv0.rds")
check_hmc_diagnostics(cod_surv0$fit)
neff_lowest(cod_surv0$fit)
rhat_highest(cod_surv0$fit)
summary(cod_surv0)
bayes_R2(cod_surv0)
plot(cod_surv0$criteria$loo, "k")
plot(conditional_smooths(cod_surv0), ask = FALSE)
y <- surv.data$surv.rate_stnd
yrep_cod_surv0  <- fitted(cod_surv0, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_surv0[sample(nrow(yrep_cod_surv0), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_surv0")


cod_surv1 <- brm(surv1_formula,
                   data = surv.data,
                   cores = 4, chains = 4, iter = 4000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.99, max_treedepth = 10))
cod_surv1  <- add_criterion(cod_surv1, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_surv1, file = "output/cod_surv1.rds")

cod_surv1 <- readRDS("./output/cod_surv1.rds")
check_hmc_diagnostics(cod_surv1$fit)
neff_lowest(cod_surv1$fit)
rhat_highest(cod_surv1$fit)
summary(cod_surv1)
bayes_R2(cod_surv1)
plot(cod_surv1$criteria$loo, "k")
plot(conditional_smooths(cod_surv1), ask = FALSE)
y <- surv.data$surv.rate_stnd
yrep_cod_surv1  <- fitted(cod_surv1, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_surv1[sample(nrow(yrep_cod_surv1), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_surv1")


cod_surv2 <- brm(surv2_formula,
                   data = surv.data,
                   cores = 4, chains = 4, iter = 4000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.99, max_treedepth = 10))
cod_surv2  <- add_criterion(cod_surv2, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_surv2, file = "output/cod_surv2.rds")

cod_surv2 <- readRDS("./output/cod_surv2.rds")
check_hmc_diagnostics(cod_surv2$fit)
neff_lowest(cod_surv2$fit)
rhat_highest(cod_surv2$fit)
summary(cod_surv2)
bayes_R2(cod_surv2)
plot(cod_surv2$criteria$loo, "k")
plot(conditional_smooths(cod_surv2), ask = FALSE)
y <- surv.data$surv.rate_stnd
yrep_cod_surv2  <- fitted(cod_surv2, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_surv2[sample(nrow(yrep_cod_surv2), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_surv2")



cod_surv3 <- brm(surv3_formula,
                    data = surv.data,
                    cores = 4, chains = 4, iter = 4000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.9999, max_treedepth = 16))
cod_surv3  <- add_criterion(cod_surv3, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_surv3, file = "output/cod_surv3.rds")

cod_surv3 <- readRDS("./output/cod_surv3.rds")
check_hmc_diagnostics(cod_surv3$fit)
neff_lowest(cod_surv3$fit)
rhat_highest(cod_surv3$fit)
summary(cod_surv3)
bayes_R2(cod_surv3)
plot(cod_surv3$criteria$loo, "k")
plot(conditional_smooths(cod_surv3), ask = FALSE)
y <- surv.data$surv.rate_stnd
yrep_cod_surv3  <- fitted(cod_surv3, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_surv3[sample(nrow(yrep_cod_surv3), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_surv3")


# Model comparison
cod_surv0 <- readRDS("./output/cod_surv0.rds")
cod_surv1 <- readRDS("./output/cod_surv1.rds")
cod_surv2 <- readRDS("./output/cod_surv2.rds")
cod_surv3 <- readRDS("./output/cod_surv3.rds")

model.comp <- loo(cod_surv0, cod_surv1, cod_surv2, cod_surv3) # model 3!

model.comp

# save model comparison for ms.
forms <- data.frame(formula=c(as.character(surv3_formula)[1],
                              as.character(surv2_formula)[1],
                              as.character(surv1_formula)[1],
                              as.character(surv0_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])
write.csv(comp.out, "./output/survival_model_comp.csv")
## temp predictions ##

## 95% CI
ce1s_1 <- conditional_effects(cod_surv3, effect = "fit.temp.mu", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod_surv3, effect = "fit.temp.mu", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod_surv3, effect = "fit.temp.mu", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$fit.temp.mu
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$fit.temp.mu[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$fit.temp.mu[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$fit.temp.mu[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$fit.temp.mu[["lower__"]]
dat_ce[["rug.anom"]] <- c(unique(surv.data$fit.temp.mu), rep(NA, 100-length(unique(surv.data$fit.temp.mu))))

g1 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Summer temperature (ºC)", y = "Abundance change (anomaly)") +
  ylim(-3.7,1.2) +
  geom_rug(aes(x=rug.anom, y=NULL)) +
  theme_bw() 
print(g1)


## mean.fourth.root.cpue predicted effect
## 95% CI
ce1s_1 <- conditional_effects(cod_surv3, effect = "mean.cpue.4", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod_surv3, effect = "mean.cpue.4", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod_surv3, effect = "mean.cpue.4", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$mean.cpue.4
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$mean.cpue.4[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$mean.cpue.4[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$mean.cpue.4[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$mean.cpue.4[["lower__"]]

# dat_ce[["rug.anom"]] <- c(unique(surv.data$mean.cpue.4), rep(NA, 100-length(unique(surv.data$mean.cpue.4))))

rug.anom <- data.frame(mean.cpue.4 = unique(surv.data$mean.cpue.4))

g2 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Fourth root CPUE", y = "Abundance change (anomaly)") +
  geom_rug(data = rug.anom, aes(x=mean.cpue.4, y=NULL)) +
  ylim(-3.7,1.2) +
  theme_bw() 

print(g2)

## and predict bay value!
ce1s_1 <- conditional_effects(cod_surv3, effect = "bay_fac", probs = c(0.025, 0.975))
mod.95 <- ce1s_1$bay_fac %>%
  select(bay_fac, estimate__, lower__, upper__)
names(mod.95)[3:4] <- c("ymin.95", "ymax.95")

# set the bays to plot west - east
order <- read.csv("./data/bay_lat_long.csv")

# fix names
order$Bay[1:2] <- c("Anton Larsen", "Cook")

mod.95$bay_fac <- as.character(mod.95$bay_fac)

change <- grep("Cook", mod.95$bay_fac)
mod.95$bay_fac[change] <- "Cook"

change <- grep("Anton", mod.95$bay_fac)
mod.95$bay_fac[change] <- "Anton Larsen"
mod.95$bay_fac <- as.factor(mod.95$bay_fac)

mod.95$long <- order$lon[match(mod.95$bay_fac, order$Bay)]
mod.95$bay_fac <- reorder(mod.95$bay_fac, desc(mod.95$long))

theme_set(theme_bw())

g3 <- ggplot(mod.95) +
  aes(x = bay_fac, y = estimate__) +
  geom_errorbar(aes(ymin = ymin.95, ymax = ymax.95), width = 0.5) +
  geom_point(size = 3) +
  theme(axis.text.x = element_text(angle=30, vjust=1, hjust=1)) +
  ylab("      Abundance change 
       (anomaly)") +
  xlab("Bay")

print(g3)

png("figs/Fig8_survival_3.png", 10.5, 3, units='in', res=300)
ggpubr::ggarrange(g1, g2, g3,
                  ncol = 3,
                  widths = c(1,1,1.2),
                  labels = "auto")
dev.off()


## add scatter plot for SI -------------------------------

SI.dat <- surv.data

SI.g1 <- ggplot(SI.dat, aes(surv.rate)) +
  geom_histogram(bins=60, fill = "grey", color = "dark grey") +
  xlab("Fourth root CPUE trend (proportion change / d)") +
  ylab("Count") +
  geom_vline(xintercept = 0, lty = 2)
  
SI.g1
  
# fix Anton and Cook
change <- grep("Cook", SI.dat$bay)
SI.dat$bay[change] <- "Cook"

change <- grep("Anton", SI.dat$bay)
SI.dat$bay[change] <- "Anton Larsen"

SI.dat$jitter.surv.rate <- jitter(SI.dat$surv.rate, factor = 2)
SI.dat$jitter.fit.temp.mu <- jitter(SI.dat$fit.temp.mu, factor = 40)

SI.g2 <- ggplot(SI.dat, aes(jitter.fit.temp.mu, jitter.surv.rate)) +
  geom_point(size = 1, alpha = 0.5) +
  labs(x = "Summer temperature (ºC)", y = "Fourth root CPUE trend 
       (proportion change / d)       ") +
  theme_bw() 

SI.g2

SI.dat$jitter.fourth.root.cpue <- jitter(SI.dat$first.cpue.4, factor = 10)

SI.g3 <- ggplot(SI.dat) +
  aes(jitter.fourth.root.cpue, jitter.surv.rate) +
  geom_point(size = 1, alpha = 0.5) +
  labs(x = "Fourth root CPUE", y = "Fourth root CPUE trend 
       (proportion change / d)       ") +
  theme_bw() 

print(SI.g3)

# set the bays to plot west - east
order <- read.csv("./data/bay_lat_long.csv")
order$Bay[1:2] <- c("Anton Larsen", "Cook")

SI.dat$long <- order$lon[match(SI.dat$bay, order$Bay)]
SI.dat$bay <- reorder(SI.dat$bay, desc(SI.dat$long))
SI.dat$bay.number <- as.numeric(SI.dat$bay)
SI.dat$jitter.bay <- jitter(SI.dat$bay.number, factor = 1)

SI.g4 <- ggplot(SI.dat) +
  aes(jitter.bay, jitter.surv.rate) +
  geom_point(size = 1, alpha = 0.5) +
  theme(axis.text.x = element_text(angle=30, vjust=1, hjust=1)) +
  scale_x_continuous(breaks=c(1:11), minor_breaks = NULL, labels = levels(SI.dat$bay)) +
  ylab("Fourth root CPUE trend 
       (proportion change / d)       ") +
  xlab("Bay")

print(SI.g4)

## save

png("figs/SI_Fig_survival_3.png", 9.5, 6.5, units='in', res=300)
ggpubr::ggarrange(SI.g1, SI.g2, SI.g3, SI.g4,
                  ncol = 2, nrow = 2,
                  labels = "auto")
dev.off()
