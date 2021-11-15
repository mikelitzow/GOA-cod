## Cod length
## gamma distribution

library(ggplot2)
library(plyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
library(tidyverse)
source("./scripts/stan_utils.R")


## Read in data --------------------------------------------

## note that this is the same data QA/QC version used in cod_length.R
cod.length.data <- read.csv("data/cod length data.csv", row.names = 1)

# remove NAs!
cod.length.data <- na.omit(cod.length.data)

cod.length.data$log_length <- log(cod.length.data$length)
cod.length.data$bay_fac <- as.factor(cod.length.data$bay)
cod.length.data$year_fac <- as.factor(cod.length.data$year)
cod.length.data$site_fac <- as.factor(cod.length.data$site)
cod.length.data$bay_site_fac <- as.factor(paste0(cod.length.data$bay, "_", cod.length.data$site))
cod.length.data$date <- as.Date(cod.length.data$julian,
                                origin = paste0(cod.length.data$year, "-01-01"))

range(cod.length.data$cod.cpue)

# shouldn't be a 0 cpue! examine
filter(cod.length.data, cod.cpue==0)

cc <- read.csv("data/cpue.data.csv")

test <- cod.length.data %>%
  select(date, year, julian, site, bay, species, length, cod.cpue) %>%
  filter(species == "pcod")

test <- left_join(test, cc)

filter(test, cod.cpue==0)

# confirms that cpue = 0 (or NA) for these
# something is wrong with those - remove!

cod.length.data <- cod.length.data %>%
  filter(cod.cpue > 0)

# fourth-root transform abundance
# note that we aren't controlling abundance for date in this case as the
# date range for analysis is so small (doy 184-200)
cod.length.data$fourth.root.cpue <- cod.length.data$cod.cpue^0.25

# read in temperature - GODAS egg/larval anomalies
temp <- read.csv("./data/godas anomalies.csv", row.names = 1) 
head(temp)

temp <- temp %>%
  select(year, mean.anom)

names(temp)[2] <- "temp.anom"

cod.length.data <- left_join(cod.length.data, temp)

head(cod.length.data)

range(cod.length.data$cod.cpue)

## Remove age > 0
cod.length.data <- cod.length.data %>%
  filter(length %in% 10:149)

# look at data from doy <= 200
cod.length.data <- cod.length.data %>%
  filter(julian <= 200)

# get sample size
nrow(cod.length.data)

plot(cod.length.data$julian, cod.length.data$length) # looks good
plot(cod.length.data$julian, cod.length.data$fourth.root.cpue)
## Check distributions
plot(cod.length.data$length)
hist(cod.length.data$length, breaks = 80) # pretty long tail!
hist(cod.length.data$length[cod.length.data$length < 100], breaks = 80)

# not going to remove the ~100mm fish - will assume they're
# early settling / fast growing age-0!

tab <- table(cod.length.data$length)
plot(tab)
summary(stats::glm(length ~ 1, data = cod.length.data, family = Gamma))

hist(cod.length.data$julian, breaks = 30) 


## brms: setup ---------------------------------------------

## Define model formulas
length0_formula <- bf(length ~ s(julian, k = 4) + s(temp.anom, k = 4))

length1_formula <-  bf(length ~ s(julian, k = 4) + s(temp.anom, k = 4) + s(fourth.root.cpue, k = 4))

length2_formula <-  bf(length ~ s(julian, k = 4) + s(temp.anom, k = 4) + s(fourth.root.cpue, k = 4) + bay_fac)

## Set model distribution
Gamma <- Gamma(link = "log")

## Show default priors
get_prior(length1_formula, cod.length.data)

# using defaults for now!

## fit: brms --------------------------------------
cod_gamma_settle_len0 <- brm(length0_formula,
                             data = cod.length.data,
                             family = Gamma,
                             cores = 4, chains = 4, iter = 4000,
                             save_pars = save_pars(all = TRUE),
                             control = list(adapt_delta = 0.999, max_treedepth = 12))
cod_gamma_settle_len0  <- add_criterion(cod_gamma_settle_len0, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_gamma_settle_len0, file = "output/cod_gamma_settle_len0.rds")

cod_gamma_settle_len0 <- readRDS("./output/cod_gamma_settle_len0.rds")
check_hmc_diagnostics(cod_gamma_settle_len0$fit)
neff_lowest(cod_gamma_settle_len0$fit)
rhat_highest(cod_gamma_settle_len0$fit)
summary(cod_gamma_settle_len0)
bayes_R2(cod_gamma_settle_len0)
plot(cod_gamma_settle_len0$criteria$loo, "k")
plot(conditional_smooths(cod_gamma_settle_len0), ask = FALSE)
y <- cod.length.data $length
yrep_cod_gamma_settle_len0  <- fitted(cod_gamma_settle_len0, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_gamma_settle_len0[sample(nrow(yrep_cod_gamma_settle_len0), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_gamma_settle_len0") # also bad!


cod_gamma_settle_len1 <- brm(length1_formula,
                data = cod.length.data,
                family = Gamma,
                cores = 4, chains = 4, iter = 4000,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.999, max_treedepth = 12))
cod_gamma_settle_len1  <- add_criterion(cod_gamma_settle_len1, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_gamma_settle_len1, file = "output/cod_gamma_settle_len1.rds")

cod_gamma_settle_len1 <- readRDS("./output/cod_gamma_settle_len1.rds")
check_hmc_diagnostics(cod_gamma_settle_len1$fit)
neff_lowest(cod_gamma_settle_len1$fit)
rhat_highest(cod_gamma_settle_len1$fit)
summary(cod_gamma_settle_len1)
bayes_R2(cod_gamma_settle_len1)
plot(cod_gamma_settle_len1$criteria$loo, "k")
plot(conditional_smooths(cod_gamma_settle_len1), ask = FALSE)
y <- cod.length.data $length
yrep_cod_gamma_settle_len1  <- fitted(cod_gamma_settle_len1, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_gamma_settle_len1[sample(nrow(yrep_cod_gamma_settle_len1), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_gamma_settle_len1")


## adding model with bay factor
cod_gamma_settle_len2 <- brm(length2_formula,
                             data = cod.length.data,
                             family = Gamma,
                             cores = 4, chains = 4, iter = 2500,
                             save_pars = save_pars(all = TRUE),
                             control = list(adapt_delta = 0.99999, max_treedepth = 16))
cod_gamma_settle_len2  <- add_criterion(cod_gamma_settle_len2, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_gamma_settle_len2, file = "output/cod_gamma_settle_len2.rds")

cod_gamma_settle_len2 <- readRDS("./output/cod_gamma_settle_len2.rds")
check_hmc_diagnostics(cod_gamma_settle_len2$fit)
neff_lowest(cod_gamma_settle_len2$fit)
rhat_highest(cod_gamma_settle_len2$fit)
summary(cod_gamma_settle_len2)
bayes_R2(cod_gamma_settle_len2)
plot(cod_gamma_settle_len2$criteria$loo, "k")
plot(conditional_smooths(cod_gamma_settle_len2), ask = FALSE)
y <- cod.length.data $length
yrep_cod_gamma_settle_len2  <- fitted(cod_gamma_settle_len2, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_gamma_settle_len2[sample(nrow(yrep_cod_gamma_settle_len2), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_gamma_settle_len2")


# Model comparison
cod_gamma_settle_len0 <- readRDS("./output/cod_gamma_settle_len0.rds")
cod_gamma_settle_len1 <- readRDS("./output/cod_gamma_settle_len1.rds")
cod_gamma_settle_len2 <- readRDS("./output/cod_gamma_settle_len2.rds")


loo(cod_gamma_settle_len0, cod_gamma_settle_len1, cod_gamma_settle_len2)


# save model comparison output
model.comp <- loo(cod_gamma_settle_len0, cod_gamma_settle_len1, cod_gamma_settle_len2)
model.comp

# save model comparison for ms.
forms <- data.frame(formula=c(as.character(length2_formula)[1],
                              as.character(length1_formula)[1],
                              as.character(length0_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])
write.csv(comp.out, "./output/settlement_length_model_comp.csv")

## temp.anom predictions ##

## 95% CI
ce1s_1 <- conditional_effects(cod_gamma_settle_len2, effect = "temp.anom", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod_gamma_settle_len2, effect = "temp.anom", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod_gamma_settle_len2, effect = "temp.anom", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$temp.anom
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$temp.anom[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$temp.anom[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$temp.anom[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$temp.anom[["lower__"]]
dat_ce[["rug.anom"]] <- c(unique(cod.length.data$temp.anom), rep(NA, 100-length(unique(cod.length.data$temp.anom))))

g1 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Egg/larval temperature anomaly", y = "Total length (mm)") +
  theme_bw() +
  geom_rug(aes(x=rug.anom, y=NULL)) +
  ylim(30, 60)

print(g1)

## fourth.root.cpue predicted effect
## 95% CI
ce1s_1 <- conditional_effects(cod_gamma_settle_len2, effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod_gamma_settle_len2, effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod_gamma_settle_len2, effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$fourth.root.cpue
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$fourth.root.cpue[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$fourth.root.cpue[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$fourth.root.cpue[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$fourth.root.cpue[["lower__"]]
dat_ce[["rug.anom"]] <- c(unique(cod.length.data$fourth.root.cpue), rep(NA, 100-length(unique(cod.length.data$fourth.root.cpue))))

g2 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Fourth root CPUE", y = "Total length (mm)") +
  theme_bw() +
  geom_rug(aes(x=rug.anom, y=NULL)) +
  ylim(30, 60)

print(g2)

## julian predicted effect
## 95% CI
ce1s_1 <- conditional_effects(cod_gamma_settle_len2, effect = "julian", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod_gamma_settle_len2, effect = "julian", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod_gamma_settle_len2, effect = "julian", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$julian
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$julian[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$julian[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$julian[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$julian[["lower__"]]
dat_ce[["rug.anom"]] <- c(unique(cod.length.data$julian), rep(NA, 100-length(unique(cod.length.data$julian))))

g3 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Day of year", y = "Total length (mm)") +
  theme_bw() +
  geom_rug(aes(x=rug.anom, y=NULL)) +
  ylim(30, 60)

print(g3)

## and predict bay value!
ce1s_1 <- conditional_effects(cod_gamma_settle_len2, effect = "bay_fac", probs = c(0.025, 0.975))
mod.95 <- ce1s_1$bay_fac %>%
  select(bay_fac, estimate__, lower__, upper__)
names(mod.95)[3:4] <- c("ymin.95", "ymax.95")

# set the bays to plot west - east
order <- read.csv("./data/bay_lat_long.csv")

# remove "Bay" from first two bay names
order$Bay[1:2] <- c("Anton Larsen", "Cook")

mod.95$bay_fac <- as.character(mod.95$bay_fac)
change <- grep("Anton", mod.95$bay_fac)
mod.95$bay_fac[change] <- "Anton Larsen"

change <- grep("Cook", mod.95$bay_fac)
mod.95$bay_fac[change] <- "Cook"

mod.95$long <- order$lon[match(mod.95$bay_fac, order$Bay)]
mod.95$bay_fac <- reorder(mod.95$bay_fac, desc(mod.95$long))

theme_set(theme_bw())

g4 <- ggplot(mod.95) +
  aes(x = bay_fac, y = estimate__) +
  geom_errorbar(aes(ymin = ymin.95, ymax = ymax.95), width = 0.5) +
  geom_point(size = 3) +
  theme(axis.text.x = element_text(angle=30, vjust=1, hjust=1)) +
  ylab("Total length (mm)") +
  xlab("Bay") +
  ylim(30, 60) 

print(g4)


## combine plots

png("figs/Fig5_combined_predicted_length_cod_gamma_settle_len2.png", 8, 6, units='in', res=300)
ggpubr::ggarrange(g1, g2, g3, g4, ncol=2, nrow=2, labels="auto")
dev.off()

## add scatter plot for SI -------------------------------

SI.dat <- cod.length.data %>%
  select(bay, julian, length, temp.anom, fourth.root.cpue)

# fix Anton and Cook
change <- grep("Anton", SI.dat$bay)
SI.dat$bay[change] <- "Anton Larsen"

change <- grep("Cook", SI.dat$bay)
SI.dat$bay[change] <- "Cook"

SI.dat$jitter.length <- jitter(SI.dat$length, factor = 2)
SI.dat$jitter.temp.anom <- jitter(SI.dat$temp.anom, factor = 40)

g1 <- ggplot(SI.dat, aes(jitter.temp.anom, jitter.length)) +
  geom_point(size = 1, alpha = 0.5) +
  labs(x = "Egg/larval temperature anomaly", y = "Total length (mm)") +
  theme_bw() 

g1

SI.dat$jitter.julian <- jitter(SI.dat$julian, factor = 2)

g2 <- ggplot(SI.dat) +
  aes(jitter.julian, jitter.length) +
  geom_point(size = 1, alpha = 0.5) +
  labs(x = "Day of year", y = "Total length (mm)") +
  theme_bw() 

print(g2)

SI.dat$jitter.fourth.root.cpue <- jitter(SI.dat$fourth.root.cpue, factor = 10)
g3 <- ggplot(SI.dat) +
  aes(jitter.fourth.root.cpue, jitter.length) +
  geom_point(size = 1, alpha = 0.5) +
  labs(x = "Fourth root CPUE", y = "Total length (mm)") +
  theme_bw() 

print(g3)

# set the bays to plot west - east
order <- read.csv("./data/bay_lat_long.csv")
order$Bay[1:2] <- c("Anton Larsen", "Cook")

SI.dat$long <- order$lon[match(SI.dat$bay, order$Bay)]
SI.dat$bay <- reorder(SI.dat$bay, desc(SI.dat$long))
SI.dat$bay.number <- as.numeric(SI.dat$bay)
SI.dat$jitter.bay <- jitter(SI.dat$bay.number, factor = 1.5)

g4 <- ggplot(SI.dat) +
  aes(jitter.bay, jitter.length) +
  geom_point(size = 1, alpha = 0.5) +
  theme(axis.text.x = element_text(angle=30, vjust=1, hjust=1)) +
  coord_trans(y = "pseudo_log") +
  scale_x_continuous(breaks=c(1:12), minor_breaks = NULL, labels = levels(SI.dat$bay)) +
  ylab("Total length (mm)") +
  xlab("Bay")

print(g4)

png("./figs/SI_length_settlement_scatter_plots.png", 8, 6, units='in', res=300)
ggpubr::ggarrange(g1, g3, g2, g4, ncol=2, nrow=2, labels=c("a", "b", "c", "d"))
dev.off()
