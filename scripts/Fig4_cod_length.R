## Cod length

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

# save to send a copy to Ben
test <- test %>%
  filter(cod.cpue == 0) %>%
  select(date, year, julian, site, bay, species, length, cod.cpue)

write.csv(test, "./data/lengths_with_zero_cpue.csv")

# something is wrong with those - remove!

cod.length.data <- cod.length.data %>%
  filter(cod.cpue > 0)

# fourth-root transform abundance
cod.length.data$fourth.root.cpue <- cod.length.data$cod.cpue^0.25

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

# check fitted mu
plot(cod.length.data$fit.temp.mu)
hist(cod.length.data$fit.temp.mu, breaks = 15) 

range(cod.length.data$cod.cpue)

# and predict cpue for date!
# using residuals

# first, get a set of dates / cpue values without repeats!
temp.dat <- cod.length.data %>%
  mutate(key = paste(bay_site_fac, julian, year, sep = "_")) %>%
  group_by(key) %>%
  summarize(julian = mean(julian),
            fourth.root.cpue = mean(fourth.root.cpue))

# this is sloppy, but add bay_site_fac and year back in to temp.dat so
# we can include those in the model
other.dat <- cod.length.data %>%
  mutate(key = paste(bay_site_fac, julian, year, sep = "_")) %>%
  select(key, bay_site_fac, year_fac)

temp.dat <- left_join(temp.dat, other.dat)

# use this no-repeat data set to fit the model
c0 <- mgcv::gam(fourth.root.cpue ~ s(julian, k = 4) + bay_site_fac + year_fac, data = temp.dat)
summary(c0)
plot(c0, resid = T, pch = 19)

# so we have a predicted abundance-date relationship
# now we want to calculate residuals from this model for the actual data - 
# BUT! we want to hold the bay_site_fac and year_fac effects constant when 
# calculating the residuals - i.e., is a particular set high or low abundance
# for that particular date ONLY, not that date in that year at that bay_site...

new.dat <- data.frame(julian = temp.dat$julian,
                      fourth.root.cpue = temp.dat$fourth.root.cpue,
                      bay_site_fac = "Agripina_AG-1",
                      year_fac = "2006")

# get predicted values
predicted.cpue <- predict(c0, type = "response", newdata = new.dat)

# and residuals (observed - predicted) - so the sign is correct...
# if resid is positive, then abundance is higher than predicted for that day
cod.length.data$residual.cpue <- cod.length.data$fourth.root.cpue - predicted.cpue

hist(cod.length.data$residual.cpue)

## Remove age > 0
cod.length.data <- cod.length.data %>%
  filter(length %in% 10:149)

plot(cod.length.data$julian, cod.length.data$length) # looks good

## Check distributions
plot(cod.length.data$length)
hist(cod.length.data$length, breaks = 100) 
tab <- table(cod.length.data$length)
plot(tab)
summary(stats::glm(length ~ 1, data = cod.length.data, family = gaussian))
summary(MASS::glm.nb(length ~ 1, data = cod.length.data))

g <- ggplot(cod.length.data) +
  aes(x = date, y = length, color = site) +
  geom_point() +
  facet_wrap( ~ bay) +
  theme(legend.position = "none")
print(g)

g <- ggplot(cod.length.data) +
  aes(x = fit.temp.mu, y = length) +
  geom_point() +
  geom_smooth(method = 'lm')
print(g)
ggsave("./figs/cod_length_temp.png", width = 5, height = 5)

g <- ggplot(cod.length.data) +
  aes(x = fit.temp.mu, y = log_length) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_wrap( ~ bay_fac)
print(g)

g <- ggplot(cod.length.data) +
  aes(x = julian, y = length) +
  geom_point() +
  geom_smooth(method = 'lm')
print(g)
ggsave("./figs/cod_length_julian.png", width = 5, height = 5)

g <- ggplot(cod.length.data) +
  aes(x = julian, y = log_length) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_wrap( ~ bay_fac)
print(g)

g <- ggplot(cod.length.data) +
  aes(x = residual.cpue, y = length) +
  geom_point() +
  geom_smooth(method = 'lm')
print(g)
ggsave("./figs/cod_length_residual.cpue.png", width = 5, height = 5)

## brms: setup ---------------------------------------------

## Define model formula
length_formula <-  bf(length ~ s(julian, k = 4) + s(fit.temp.mu, k = 4) + s(residual.cpue, k = 4) + bay_fac)

## Show default priors
get_prior(length_formula, cod.length.data)


priors_len <- c(set_prior("normal(0, 3)", class = "b"),
                set_prior("normal(0, 3)", class = "Intercept"),
                set_prior("student_t(3, 0, 3)", class = "sds"),
                set_prior("student_t(3, 0, 3)", class = "sigma"))

## fit: brms --------------------------------------
cod_length <- brm(length_formula,
                 data = cod.length.data,
                 cores = 4, chains = 4, iter = 3000,
                 # prior = priors_len,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.9, max_treedepth = 10))
cod_length  <- add_criterion(cod_length, "bayes_R2", moment_match = TRUE)
saveRDS(cod_length, file = "output/cod_length.rds")

cod_length <- readRDS("./output/cod_length.rds")
check_hmc_diagnostics(cod_length$fit)
neff_lowest(cod_length$fit)
rhat_highest(cod_length$fit)
summary(cod_length)
bayes_R2(cod_length)
plot(cod_length$criteria$loo, "k")
plot(conditional_smooths(cod_length), ask = FALSE)
y <- cod.length.data $cod
yrep_cod_length  <- fitted(cod_length, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_length[sample(nrow(yrep_cod_length), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_length")

## SST predictions ##

## 95% CI
ce1s_1 <- conditional_effects(cod_length, effect = "fit.temp.mu", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod_length, effect = "fit.temp.mu", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod_length, effect = "fit.temp.mu", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$fit.temp.mu
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$fit.temp.mu[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$fit.temp.mu[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$fit.temp.mu[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$fit.temp.mu[["lower__"]]

g <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1.5, color = "red3") +
  labs(x = "fit.temp.mu", y = "Cod length") +
  ylim(45,70) +
  theme_bw()
print(g)
ggsave("./figs/sst_predicted_effect_cod_length.png", width = 5, height = 4)


## 95% CI
ce1s_1 <- conditional_effects(cod_length, effect = "residual.cpue", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod_length, effect = "residual.cpue", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod_length, effect = "residual.cpue", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$residual.cpue
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$residual.cpue[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$residual.cpue[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$residual.cpue[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$residual.cpue[["lower__"]]

g <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1.5, color = "red3") +
  labs(x = "residual.cpue", y = "Cod length") +
  ylim(45,70) +
  theme_bw()
print(g)
ggsave("./figs/residual.cpue_predicted_effect_cod_length.png", width = 5, height = 4)

## predict length by day for Julian day ---------------

## 95% CI
ce1s_1 <- conditional_effects(cod_length, effect = "julian", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod_length, effect = "julian", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod_length, effect = "julian", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$julian
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$julian[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$julian[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$julian[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$julian[["lower__"]]

g <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Day of year", y = "Total length (mm)") +
  scale_x_continuous(breaks = seq(190, 240, 10)) +
  theme_bw()
print(g)
ggsave("./figs/Fig4_predicted_length_DOY_cod_lengths.png", width = 5, height = 3.5)

## summarize "linear phase" growth rate for ms.

linear_growth <- dat_ce %>%
  filter((julian >= 205)  & (julian <= 225))

mod <- lm(estimate__ ~ julian, data = linear_growth)
summary(mod) # ~ 0.81 mm/day
