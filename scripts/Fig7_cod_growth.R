## Cod growth between repeat visits to sites within years

library(ggplot2)
library(plyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
library(tidyverse)
source("./scripts/stan_utils.R")


## Read in data --------------------------------------------
cod.length.data <- read.csv("data/cod length data.csv", row.names = 1)

# remove NAs!
cod.length.data <- na.omit(cod.length.data)

cod.length.data$bay_fac <- as.factor(cod.length.data$bay)
cod.length.data$year_fac <- as.factor(cod.length.data$year)
cod.length.data$site_fac <- as.factor(cod.length.data$site)
cod.length.data$bay_site_fac <- as.factor(paste0(cod.length.data$bay, "_", cod.length.data$site))
cod.length.data$date <- as.Date(cod.length.data$julian,
                                origin = paste0(cod.length.data$year, "-01-01"))

# read in temperature model
seine.temp <- readRDS("./output/temp3_gauss.rds") 
summary(seine.temp)

new.dat <- cod.length.data %>%
  select(bay_fac, year_fac)
new.dat$julian <- mean(cod.length.data$julian)

fitted.temp <- as.data.frame(fitted(seine.temp, newdata = new.dat))

nrow(cod.length.data)
nrow(fitted.temp)

cod.length.data$fit.temp.mu <- fitted.temp$Estimate
cod.length.data$fit.temp.sd <- fitted.temp$Est.Error

# check fitted mu
plot(cod.length.data$fit.temp.mu)
hist(cod.length.data$fit.temp.mu, breaks = 15) 

range(cod.length.data$cod.cpue)

# shouldn't be a 0 cpue! examine
filter(cod.length.data, cod.cpue==0)

cc <- read.csv("data/cpue.data.csv")

# something is wrong with those - remove!
cod.length.data <- cod.length.data %>%
  filter(cod.cpue > 0)

# switching to 4th-root cpue instead of ln(cpue) - better distribution
cod.length.data$cpue.4 <- cod.length.data$cod.cpue^0.25

# now restrict to sites sampled twice in a year!
ff <- function(x) length(unique(x))

check <- cod.length.data %>%
  group_by(site_fac, year_fac) %>%
  summarise(visits=ff(julian))
View(check)

check$site_year_fac <- as.factor(paste0(check$site_fac, "_", as.character(check$year_fac)))

cod.length.data$site_year_fac <- as.factor(paste0(cod.length.data$site, "_", cod.length.data$year))

cod.length.data <- left_join(cod.length.data, check)

growth.data <- cod.length.data %>%
  filter(visits > 1)
  
## Remove age > 0
growth.data <- growth.data %>%
  filter(length %in% 10:149)

## Check distributions
plot(growth.data$length)

# check sample sizes
sum(!is.na(growth.data$length))

check <- growth.data %>%
  group_by(site_fac, year_fac) %>%
  summarise(visits=ff(julian))
View(check)

# save for future reference
write.csv(check, "./output/bay_visits_growth_analysis.csv")

ggplot(filter(cod.length.data, length %in% 10:149), aes(length)) +
  geom_histogram(fill="grey", color="black", bins=40) +
  theme_bw() +
  ggtitle("All length data")

ggplot(growth.data, aes(length)) +
  geom_histogram(fill="grey", color="black", bins=40) +
  theme_bw() +
  ggtitle("Growth data (sites w/ > 1 visit)")


tab <- table(growth.data$length)
plot(tab)
summary(stats::glm(length ~ 1, data = growth.data, family = Gamma))

hist(growth.data$julian, breaks = 30) 


g <- ggplot(growth.data) +
  aes(x = date, y = length, color = site) +
  geom_point() +
  facet_wrap( ~ bay) +
  theme(legend.position = "none")
print(g)

g <- ggplot(growth.data) +
  aes(x = fit.temp.mu, y = length) +
  geom_point()
print(g)

g <- ggplot(growth.data) +
  aes(x = julian, y = length) +
  geom_point()
print(g)

# get initial size for each site in each year!
first <- growth.data %>%
  group_by(site_year_fac) %>%
  summarise(julian=min(julian))

# separate out the length data for the first visit to each site / year combo
first.dat <- left_join(first, growth.data)

hist(first.dat$length, breaks = 30) # pretty normal!
hist(first.dat$julian, breaks = 30) 

# get the mean length at first visit
first.mean <- first.dat %>%
  group_by(site_year_fac) %>%
  summarise(mean.settle.length=mean(length),
            julian.first=mean(julian))

hist(first.mean$mean.settle.length, breaks = 30)

# and join it back in
growth.data <- left_join(growth.data, first.mean)

# and get the mean cpue^0.25 at every visit
density <- growth.data %>%
  group_by(site_year_fac) %>%
  summarise(mean.cpue.4=mean(cpue.4))

# and join it back in
growth.data <- left_join(growth.data, density)

# and get the change in date and length from the first visit
growth.data$delta.date <- growth.data$julian - growth.data$julian.first
growth.data$delta.length <- growth.data$length - growth.data$mean.settle.length

# use that to calculate growth rate (mm/d)
growth.data$growth.rate <- growth.data$delta.length / growth.data$delta.date

View(arrange(growth.data, delta.date))

## drop the delta.date = 0 (first visits)
## and delta.date <= 5 (too short for resolving signal!)
growth.data <- growth.data %>%
  filter(delta.date > 5)
  
hist(growth.data$growth.rate, breaks = 50)

ggplot(growth.data, aes(julian.first, delta.date)) +
  geom_point() +
  theme_bw()

ggplot(growth.data, aes(fit.temp.mu, mean.cpue.4)) +
  geom_point() +
  labs(x = "Summer temp", y = "Mean cpue^0.25") +
  theme_bw()

# Julian first and delta.date are highly colinear! include only Julian date

## mgcv fits -----------------------------------------------
gam.growth.0 <- gam(growth.rate ~ s(julian.first, k = 4),
                    data=growth.data)

summary(gam.growth.0)
plot(gam.growth.0, pages = 1)

gam.growth.1 <- gam(growth.rate ~ s(julian.first, k = 4) + s(fit.temp.mu, k = 4),
                    data=growth.data)

summary(gam.growth.1)
plot(gam.growth.1, pages = 1)

gam.growth.2 <- gam(growth.rate ~ s(julian.first, k = 4) + s(fit.temp.mu, k = 4) + s(mean.cpue.4, k = 4),
                    data=growth.data)

summary(gam.growth.2)
plot(gam.growth.2, pages = 1)

gam.growth.3 <- gam(growth.rate ~ s(julian.first, k = 4) + s(fit.temp.mu, k = 4) + s(mean.cpue.4, k = 4) + bay_fac,
                    data=growth.data)

summary(gam.growth.3)
plot(gam.growth.3, pages = 1)


MuMIn::AICc(gam.growth.0, gam.growth.1, gam.growth.2, gam.growth.3) # model 3 is massively the best!

## brms: setup ---------------------------------------------

# scale growth.rate
growth.data$growth.rate_stnd <- as.vector(scale(growth.data$growth.rate))

## Define model formulas
growth0_formula <-  bf(growth.rate_stnd ~ s(julian.first, k = 4) + bay_fac)

growth1_formula <-  bf(growth.rate_stnd ~ s(julian.first, k = 4) + s(fit.temp.mu, k = 4) + bay_fac)

growth2_formula <-  bf(growth.rate_stnd ~ s(julian.first, k = 4) + s(fit.temp.mu, k = 4) + s(mean.cpue.4, k = 4) + bay_fac)

growth2s_formula <-  bf(growth.rate_stnd ~ s(julian.first, k = 4) + s(fit.temp.mu, k = 4) + s(mean.cpue.4, k = 4) + 
                          (1 | bay_fac/site_fac))

growth3_formula <-  bf(growth.rate_stnd ~ s(fit.temp.mu, k = 4) + s(mean.cpue.4, k = 4) + bay_fac)

growth4_formula <-  bf(growth.rate_stnd ~ s(mean.cpue.4, k = 4) + bay_fac)

growth5_formula <-  bf(growth.rate_stnd ~ s(fit.temp.mu, k = 4) + bay_fac)

## Show default priors
get_prior(growth0_formula, growth.data)

# using defaults for now!

## fit: brms --------------------------------------
cod_growth0 <- brm(growth0_formula,
                   data = growth.data,
                   cores = 4, chains = 4, iter = 4000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 16))
cod_growth0  <- add_criterion(cod_growth0, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_growth0, file = "output/cod_growth0.rds")

cod_growth0 <- readRDS("./output/cod_growth0.rds")
check_hmc_diagnostics(cod_growth0$fit)
neff_lowest(cod_growth0$fit)
rhat_highest(cod_growth0$fit)
summary(cod_growth0)
bayes_R2(cod_growth0)
plot(cod_growth0$criteria$loo, "k")
plot(conditional_smooths(cod_growth0), ask = FALSE)
y <- growth.data$growth.rate_stnd
yrep_cod_growth0  <- fitted(cod_growth0, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_growth0[sample(nrow(yrep_cod_growth0), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_growth0")


cod_growth1 <- brm(growth1_formula,
                   data = growth.data,
                   cores = 4, chains = 4, iter = 4000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 16))
cod_growth1  <- add_criterion(cod_growth1, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_growth1, file = "output/cod_growth1.rds")

cod_growth1 <- readRDS("./output/cod_growth1.rds")
check_hmc_diagnostics(cod_growth1$fit)
neff_lowest(cod_growth1$fit)
rhat_highest(cod_growth1$fit)
summary(cod_growth1)
bayes_R2(cod_growth1)
plot(cod_growth1$criteria$loo, "k")
plot(conditional_smooths(cod_growth1), ask = FALSE)
y <- growth.data$growth.rate_stnd
yrep_cod_growth1  <- fitted(cod_growth1, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_growth1[sample(nrow(yrep_cod_growth1), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_growth1")


cod_growth2 <- brm(growth2_formula,
                   data = growth.data,
                   cores = 4, chains = 4, iter = 4000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.99, max_treedepth = 16))
cod_growth2  <- add_criterion(cod_growth2, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_growth2, file = "output/cod_growth2.rds")

cod_growth2 <- readRDS("./output/cod_growth2.rds")
check_hmc_diagnostics(cod_growth2$fit)
neff_lowest(cod_growth2$fit)
rhat_highest(cod_growth2$fit)
summary(cod_growth2)
bayes_R2(cod_growth2)
plot(cod_growth2$criteria$loo, "k")
plot(conditional_smooths(cod_growth2), ask = FALSE)
y <- growth.data$growth.rate_stnd
yrep_cod_growth2  <- fitted(cod_growth2, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_growth2[sample(nrow(yrep_cod_growth2), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_growth2")


cod_growth2s <- brm(growth2s_formula,
                    data = growth.data,
                    cores = 4, chains = 4, iter = 3000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.999, max_treedepth = 16))
cod_growth2s  <- add_criterion(cod_growth2s, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_growth2s, file = "output/cod_growth2s.rds")

cod_growth2s <- readRDS("./output/cod_growth2s.rds")
check_hmc_diagnostics(cod_growth2s$fit)
neff_lowest(cod_growth2s$fit)
rhat_highest(cod_growth2s$fit)
summary(cod_growth2s)
bayes_R2(cod_growth2s)
plot(cod_growth2s$criteria$loo, "k")
plot(conditional_smooths(cod_growth2s), ask = FALSE)
y <- growth.data$growth.rate_stnd
yrep_cod_growth2s  <- fitted(cod_growth2s, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_growth2s[sample(nrow(yrep_cod_growth2s), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_growth2s")

cod_growth3 <- brm(growth3_formula,
                   data = growth.data,
                   cores = 4, chains = 4, iter = 6000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 16))
cod_growth3  <- add_criterion(cod_growth3, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_growth3, file = "output/cod_growth3.rds")

cod_growth3 <- readRDS("./output/cod_growth3.rds")
check_hmc_diagnostics(cod_growth3$fit)
neff_lowest(cod_growth3$fit)
rhat_highest(cod_growth3$fit)
summary(cod_growth3)
bayes_R2(cod_growth3)
plot(cod_growth3$criteria$loo, "k")
plot(conditional_smooths(cod_growth3), ask = FALSE)
y <- growth.data$growth.rate_stnd
yrep_cod_growth3  <- fitted(cod_growth3, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_growth3[sample(nrow(yrep_cod_growth3), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_growth3")


cod_growth4 <- brm(growth4_formula,
                   data = growth.data,
                   cores = 4, chains = 4, iter = 4000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 12))
cod_growth4  <- add_criterion(cod_growth4, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_growth4, file = "output/cod_growth4.rds")

cod_growth4 <- readRDS("./output/cod_growth4.rds")
check_hmc_diagnostics(cod_growth4$fit)
neff_lowest(cod_growth4$fit)
rhat_highest(cod_growth4$fit)
summary(cod_growth4)
bayes_R2(cod_growth4)
plot(cod_growth4$criteria$loo, "k")
plot(conditional_smooths(cod_growth4), ask = FALSE)
y <- growth.data$growth.rate_stnd
yrep_cod_growth4  <- fitted(cod_growth4, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_growth4[sample(nrow(yrep_cod_growth4), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_growth4")


cod_growth5 <- brm(growth5_formula,
                   data = growth.data,
                   cores = 4, chains = 4, iter = 4000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.99, max_treedepth = 16))
cod_growth5  <- add_criterion(cod_growth5, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_growth5, file = "output/cod_growth5.rds")

cod_growth5 <- readRDS("./output/cod_growth5.rds")
check_hmc_diagnostics(cod_growth5$fit)
neff_lowest(cod_growth5$fit)
rhat_highest(cod_growth5$fit)
summary(cod_growth5)
bayes_R2(cod_growth5)
plot(cod_growth5$criteria$loo, "k")
plot(conditional_smooths(cod_growth5), ask = FALSE)
y <- growth.data$growth.rate_stnd
yrep_cod_growth5  <- fitted(cod_growth5, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_growth5[sample(nrow(yrep_cod_growth5), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_growth5")


# Model comparison
cod_growth0 <- readRDS("./output/cod_growth0.rds")
cod_growth1 <- readRDS("./output/cod_growth1.rds")
cod_growth2 <- readRDS("./output/cod_growth2.rds")
cod_growth2s <- readRDS("./output/cod_growth2s.rds")
cod_growth3 <- readRDS("./output/cod_growth3.rds")
cod_growth4 <- readRDS("./output/cod_growth4.rds")
cod_growth5 <- readRDS("./output/cod_growth5.rds")

model.comp <- loo(cod_growth0, cod_growth1, cod_growth2, cod_growth2s,
                  cod_growth3, cod_growth4, cod_growth5)
model.comp

# save model comparison for ms.
forms <- data.frame(formula=c(as.character(growth2s_formula)[1],
                              as.character(growth2_formula)[1],
                              as.character(growth1_formula)[1],
                              as.character(growth0_formula)[1],
                              as.character(growth3_formula)[1],
                              as.character(growth4_formula)[1],
                              as.character(growth5_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])
write.csv(comp.out, "./output/growth_model_comp.csv")

## temp predictions ##

## 95% CI
ce1s_1 <- conditional_effects(cod_growth2s, effect = "fit.temp.mu", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod_growth2s, effect = "fit.temp.mu", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod_growth2s, effect = "fit.temp.mu", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$fit.temp.mu
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$fit.temp.mu[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$fit.temp.mu[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$fit.temp.mu[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$fit.temp.mu[["lower__"]]
dat_ce[["rug.anom"]] <- c(unique(cod.length.data$fit.temp.mu), rep(NA, 100-length(unique(cod.length.data$fit.temp.mu))))

g1 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Summer temperature (ºC)", y = "Growth rate (anomaly)") +
  theme_bw() 

print(g1)


## mean.cpue.4 predicted effect
## 95% CI
ce1s_1 <- conditional_effects(cod_growth2s, effect = "mean.cpue.4", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod_growth2s, effect = "mean.cpue.4", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod_growth2s, effect = "mean.cpue.4", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$mean.cpue.4
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$mean.cpue.4[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$mean.cpue.4[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$mean.cpue.4[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$mean.cpue.4[["lower__"]]
dat_ce[["rug.anom"]] <- c(unique(cod.length.data$mean.cpue.4), rep(NA, 100-length(unique(cod.length.data$mean.cpue.4))))

g2 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Fourth root CPUE", y = "Growth rate (anomaly)") +
  theme_bw()

print(g2)


png("figs/Fig7_cod_growth_2s_temp_cpue.png", 7, 3, units='in', res=300)
ggpubr::ggarrange(g1, g2,
                  ncol = 2,
                  labels = "auto")
dev.off()

## add scatter plot for SI -------------------------------

SI.dat <- growth.data


SI.dat$jitter.growth.rate <- jitter(SI.dat$growth.rate, factor = 2)
SI.dat$jitter.fit.temp.mu <- jitter(SI.dat$fit.temp.mu, factor = 10)

SI.g1 <- ggplot(SI.dat, aes(jitter.fit.temp.mu, jitter.growth.rate)) +
  geom_point(size = 1, alpha = 0.5) +
  labs(x = "Summer temperature (ºC)", y = "Growth rate (mm / d)") +
  theme_bw() 

SI.g1

SI.dat$jitter.fourth.root.cpue <- jitter(SI.dat$cpue.4, factor = 20)

SI.g2 <- ggplot(SI.dat) +
  aes(jitter.fourth.root.cpue, jitter.growth.rate) +
  geom_point(size = 1, alpha = 0.5) +
  labs(x = "Fourth root CPUE", y = "Growth rate (mm / d)") +
  theme_bw() 

print(SI.g2)


## get predicted growth and CI for ms. --------------------------------
growth2s_formula <-  bf(growth.rate_stnd ~ s(julian.first, k = 4) + s(fit.temp.mu, k = 4) + s(mean.cpue.4, k = 4) + (1 | bay_fac/site_fac))

# predict for a single site in Cook Bay

test <- growth.data %>%
  filter(bay == "Cook Bay")

ggplot(test, aes(length)) +
  geom_histogram(color = "light gray") +
  facet_wrap(~site_fac)

new.dat <- data.frame(bay_fac = "Cook Bay",
                      site_fac = "Eelgrass South",
                      julian.first = mean(growth.data$julian.first),
                      fit.temp.mu = mean(growth.data$fit.temp.mu),
                      mean.cpue.4 = mean(growth.data$cpue.4))

pred.grow <- predict(cod_growth2s, newdata = new.dat)

## compare with cod_growth2
pred.grow <- predict(cod_growth2, newdata = new.dat)
