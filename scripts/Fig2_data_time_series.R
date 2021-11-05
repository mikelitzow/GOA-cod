# create a plot of annual estimates of abundance, size, and condition as 
# estimated by brms models

library(plyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
library(tidyverse)
source("./scripts/stan_utils.R")

theme_set(theme_bw())

## annual abundance --------------------------
## this is the best model from the fish-FAR project

# this is the formula, for reference
recr_2_formula <-  bf(cod ~ s(julian, k = 3) + (1 | bay_fac/site_fac) + (year_fac),
                      zi ~ s(julian, k = 3) + (1 | bay_fac/site_fac) + (year_fac))

# load the model object
recr_2_zinb <- readRDS("./output/recr_2_zinb.rds")

## year predictions ##

## 95% CI
ce1s_1 <- conditional_effects(recr_2_zinb, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

# ce1s_1 is far better!
print(ce1s_1)

# make a panel for the fig
plot.abundance <- ce1s_1$year_fac

abundance <- ggplot(plot.abundance, aes(year_fac, estimate__)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymax=upper__, ymin=lower__), width = 0.4) +
  coord_trans(y = "pseudo_log") +
  scale_y_continuous(breaks = c(0,1,5,10,50,100,200,500)) +
  theme(axis.title.x = element_blank()) +
  labs(y = "Fish / set")

abundance

## length------------------------------------
## read in data 
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

hist(cod.length.data$length)

# limit to cod < 150mm
cod.length.data <- cod.length.data %>%
  filter(length < 150)

# set up brms model

# define formula
length_formula <-  bf(length ~ s(julian, k = 4) + (1 | bay_fac/site_fac) + (year_fac))

## Show default priors
get_prior(length_formula, cod.length.data)

# set priors

priors_len <- c(set_prior("normal(0, 3)", class = "b"),
                set_prior("normal(0, 3)", class = "Intercept"),
                set_prior("student_t(3, 0, 2.5)", class = "sd"),
                set_prior("student_t(3, 0, 2.5)", class = "sigma"))

## fit: brms --------------------------------------
cod_length_annual <- brm(length_formula,
                data = cod.length.data,
                prior = priors_len,
                cores = 4, chains = 4, iter = 4000,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.9, max_treedepth = 12))
cod_length_annual  <- add_criterion(cod_length_annual, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_length_annual, file = "output/cod_length_annual.rds")

cod_length_annual <- readRDS("./output/cod_length_annual.rds")
check_hmc_diagnostics(cod_length_annual$fit)
neff_lowest(cod_length_annual$fit)
rhat_highest(cod_length_annual$fit)
summary(cod_length_annual)
bayes_R2(cod_length_annual)
plot(cod_length_annual$criteria$loo, "k")
plot(conditional_smooths(cod_length_annual), ask = FALSE)
y <- cod.length.data$length
yrep_cod_length_annual  <- fitted(cod_length_annual, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_length_annual[sample(nrow(yrep_cod_length_annual), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_length_annual")

# plot
## 95% CI
cod_length_annual <- readRDS("./output/cod_length_annual.rds")
ce1s_1 <- conditional_effects(cod_length_annual, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

print(ce1s_1)

# make a panel for the fig
plot.length <- ce1s_1$year_fac

# plot.dat$year <- as.numeric(as.character(plot.dat$year_fac))

length <- ggplot(plot.length, aes(year_fac, estimate__)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymax=upper__, ymin=lower__), width = 0.4) +
  theme(axis.title.x = element_blank()) +
  labs(y = "Length anomaly (mm)")

length

## fit Kdry model --------------------------------------------
## Read in data
cod.condition.data <- read.csv("data/condition_data.csv", row.names = 1)


# remove Kujulik (not used in main analysis)
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

# limit to Kdry data only
Kdry.data <- cod.condition.data %>%
  select(Kdry, julian, length, bay_fac, site_fac, year_fac) %>%
  na.omit()

# fit brms model to Kdry
# define formula
Kdry_formula <-  bf(Kdry ~ s(julian, k = 4) + s(length, k = 4) + (1 | bay_fac/site_fac) + (year_fac))

## Show default priors
get_prior(Kdry_formula, Kdry.data)

# using defaults for now!

priors_Kdry <- c(set_prior("normal(0, 3)", class = "b"),
                set_prior("normal(0, 3)", class = "Intercept"),
                set_prior("student_t(3, 0, 2.5)", class = "sd"),
                set_prior("student_t(3, 0, 2.5)", class = "sigma"))

## fit: brms --------------------------------------
cod_Kdry_annual <- brm(Kdry_formula,
                         data = Kdry.data,
                         # prior = priors_Kdry,
                         cores = 4, chains = 4, iter = 5000,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.999, max_treedepth = 12))
cod_Kdry_annual  <- add_criterion(cod_Kdry_annual, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_Kdry_annual, file = "output/cod_Kdry_annual.rds")

cod_Kdry_annual <- readRDS("./output/cod_Kdry_annual.rds")
check_hmc_diagnostics(cod_Kdry_annual$fit)
neff_lowest(cod_Kdry_annual$fit)
rhat_highest(cod_Kdry_annual$fit)
summary(cod_Kdry_annual)
bayes_R2(cod_Kdry_annual)
plot(cod_Kdry_annual$criteria$loo, "k")
plot(conditional_smooths(cod_Kdry_annual), ask = FALSE)
y <- cod.length.data$Kdry
yrep_cod_Kdry_annual  <- fitted(cod_Kdry_annual, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_Kdry_annual[sample(nrow(yrep_cod_Kdry_annual), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_Kdry_annual")

# plot
## 95% CI
cod_Kdry_annual <- readRDS("./output/cod_Kdry_annual.rds")
ce1s_1 <- conditional_effects(cod_Kdry_annual, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

print(ce1s_1)

# make a panel for the fig
plot.Kdry <- ce1s_1$year_fac


Kdry <- ggplot(plot.Kdry, aes(year_fac, estimate__)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymax=upper__, ymin=lower__), width = 0.4) +
  theme(axis.title.x = element_blank()) +
  labs(y = "Kdry")

Kdry

## fit model to HSI -----------------------------------
# limit to HSI data only
HSI.data <- cod.condition.data %>%
  select(HSI.dry, julian, length, bay_fac, site_fac, year_fac) %>%
  na.omit()

# fit brms model to HSI
# define formula
HSI_formula <-  bf(HSI.dry ~ s(julian, k = 4) + s(length, k = 4) + (1 | bay_fac/site_fac) + (year_fac))

## Show default priors
get_prior(HSI_formula, HSI.data)

# using defaults for now!

priors_HSI <- c(set_prior("normal(0, 3)", class = "b"),
                 set_prior("normal(0, 3)", class = "Intercept"),
                 set_prior("student_t(3, 0, 2.5)", class = "sd"),
                 set_prior("student_t(3, 0, 2.5)", class = "sigma"))

## fit: brms --------------------------------------
cod_HSI_annual <- brm(HSI_formula,
                       data = HSI.data,
                       prior = priors_HSI,
                       cores = 4, chains = 4, iter = 5000,
                       save_pars = save_pars(all = TRUE),
                       control = list(adapt_delta = 0.999, max_treedepth = 12))
cod_HSI_annual  <- add_criterion(cod_HSI_annual, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_HSI_annual, file = "output/cod_HSI_annual.rds")

cod_HSI_annual <- readRDS("./output/cod_HSI_annual.rds")
check_hmc_diagnostics(cod_HSI_annual$fit)
neff_lowest(cod_HSI_annual$fit)
rhat_highest(cod_HSI_annual$fit)
summary(cod_HSI_annual)
bayes_R2(cod_HSI_annual)
plot(cod_HSI_annual$criteria$loo, "k")
plot(conditional_smooths(cod_HSI_annual), ask = FALSE)
y <- cod.length.data$HSI
yrep_cod_HSI_annual  <- fitted(cod_HSI_annual, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_HSI_annual[sample(nrow(yrep_cod_HSI_annual), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_HSI_annual")

# plot
## 95% CI
cod_HSI_annual <- readRDS("./output/cod_HSI_annual.rds")
ce1s_1 <- conditional_effects(cod_HSI_annual, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

print(ce1s_1)

# make a panel for the fig
plot.HSI <- ce1s_1$year_fac


HSI <- ggplot(plot.HSI, aes(year_fac, estimate__)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymax=upper__, ymin=lower__), width = 0.4) +
  theme(axis.title.x = element_blank()) +
  labs(y = "HSI")

HSI


## combine into one-panel plot----------------

# log-transform abundance results

plot.abundance$log.est <- log(plot.abundance$estimate__, 10)
plot.abundance$log.LCI <- log(plot.abundance$lower__, 10)
plot.abundance$log.UCI <- log(plot.abundance$upper__, 10)

# and scale results to plot
sc.abundance <- plot.abundance %>%
  mutate(estimate = (log.est - mean(plot.abundance$log.est))/sd(plot.abundance$log.est),
         LCI = (log.LCI - mean(plot.abundance$log.est))/sd(plot.abundance$log.est),
         UCI = (log.UCI - mean(plot.abundance$log.est))/sd(plot.abundance$log.est),
         response = "abundance") %>%
  select(year_fac, response, estimate, LCI, UCI)

plot.abundance <- data.frame(year_fac = as.factor(2006:2020),
                             response = "abundance")
plot.abundance <- left_join(plot.abundance, sc.abundance)

# clean up length estimates
# (not sure why there is a constant length in these results, will have to ask Mike Malick!)
# won't matter for this plot as I'm scaling results

sc.length <- plot.length %>%
  mutate(estimate = (estimate__ - mean(plot.length$estimate__))/sd(plot.length$estimate__),
         LCI = (lower__ - mean(plot.length$estimate__))/sd(plot.length$estimate__),
         UCI = (upper__ - mean(plot.length$estimate__))/sd(plot.length$estimate__),
         response = "length") %>%
  select(year_fac, response, estimate, LCI, UCI)


plot.length <- data.frame(year_fac = as.factor(2006:2020),
                          response = "length")
plot.length <- left_join(plot.length, sc.length)

# and Kdry

sc.Kdry <- plot.Kdry %>%
  mutate(estimate = (estimate__ - mean(plot.Kdry$estimate__))/sd(plot.Kdry$estimate__),
         LCI = (lower__ - mean(plot.Kdry$estimate__))/sd(plot.Kdry$estimate__),
         UCI = (upper__ - mean(plot.Kdry$estimate__))/sd(plot.Kdry$estimate__),
         response = "Kdry") %>%
  select(year_fac, response, estimate, LCI, UCI)


plot.Kdry <- data.frame(year_fac = as.factor(2006:2020),
                        response = "Kdry")
plot.Kdry <- left_join(plot.Kdry, sc.Kdry)

# and HSI

sc.HSI <- plot.HSI %>%
  mutate(estimate = (estimate__ - mean(plot.HSI$estimate__))/sd(plot.HSI$estimate__),
         LCI = (lower__ - mean(plot.HSI$estimate__))/sd(plot.HSI$estimate__),
         UCI = (upper__ - mean(plot.HSI$estimate__))/sd(plot.HSI$estimate__),
         response = "HSI") %>%
  select(year_fac, response, estimate, LCI, UCI)

plot.HSI <- data.frame(year_fac = as.factor(2006:2020),
                       response = "HSI")
plot.HSI <- left_join(plot.HSI, sc.HSI)

# combine

plot.all <- rbind(plot.abundance,
                  plot.length,
                  plot.Kdry,
                  plot.HSI)

plot.all$order <- case_when(
  plot.all$response == "abundance" ~ 1,
  plot.all$response == "length" ~ 2,
  plot.all$response == "Kdry" ~ 3,
  plot.all$response == "HSI" ~ 4
)

plot.all$response <- reorder(plot.all$response, plot.all$order)

plot.all$year <- as.numeric(as.character(plot.all$year_fac))

ggplot(plot.all, aes(year, estimate, color = response)) +
  geom_point() +
  geom_line()

# set palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


ggplot(plot.all, aes(year, estimate, fill = response)) +
  geom_bar(stat = "identity", position = "dodge",
           color = "black") +
  scale_x_continuous(breaks = 2006:2020) +
  geom_hline(yintercept = 0) +
  ylab("Difference from mean (units of standard deviation)") +
  theme(legend.title = element_blank(),
        legend.position = c(0.2, 0.8),
        legend.background = element_rect(color = "black"),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = cb[c(2,4,6,8)])

ggsave("./figs/time_series_anomaly_plot.png", width = 6.5, height = 4, units = 'in')


## can also do a four-panel figure in original units -----------------

# re-organize all for outputs

# first length
ce1s_1 <- conditional_effects(recr_2_zinb, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

plot.abundance <- ce1s_1$year_fac

plot.abundance$year <- as.numeric(as.character(plot.abundance$year_fac))

abundance.2 <- ggplot(plot.abundance, aes(year, estimate__)) +
  geom_bar(stat = "identity", position = "dodge", fill = "grey80") +
  geom_errorbar(aes(ymax=upper__, ymin=lower__), width = 0.4) +
  coord_trans(y = "pseudo_log") +
  scale_y_continuous(breaks = c(0,1,5,10,50,100,200,500)) +
  theme(axis.title.x = element_blank()) +
  labs(y = "Fish / set") +
  geom_hline(yintercept = mean(plot.abundance$estimate__), lty = 2)

abundance.2

# length
ce1s_1 <- conditional_effects(cod_length_annual, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

plot.length <- ce1s_1$year_fac

plot.length$year <- as.numeric(as.character(plot.length$year_fac))

plot.length.2 <- data.frame(year = 2006:2020)

plot.length.2 <- left_join(plot.length.2, plot.length)

# add offset to get to actual units!
plot.length.2$estimate__ <- plot.length.2$estimate__ + plot.length.2$length
plot.length.2$lower__ <- plot.length.2$lower__ + plot.length.2$length
plot.length.2$upper__ <- plot.length.2$upper__ + plot.length.2$length


length.2 <- ggplot(plot.length.2, aes(year, estimate__)) +
  geom_bar(stat = "identity", position = "dodge", fill = "grey80") +
  geom_errorbar(aes(ymax=upper__, ymin=lower__), width = 0.4) +
  theme(axis.title.x = element_blank()) +
  labs(y = "Length (mm)") +
  geom_hline(yintercept = mean(plot.length.2$estimate__, na.rm = T), lty = 2) +
  coord_cartesian(ylim = c(40,85))

length.2

# Kdry

ce1s_1 <- conditional_effects(cod_Kdry_annual, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

plot.Kdry.2 <- ce1s_1$year_fac

Kdry.2 <- ggplot(plot.Kdry.2, aes(year_fac, estimate__)) +
  geom_bar(stat = "identity", position = "dodge", fill = "grey80") +
  geom_errorbar(aes(ymax=upper__, ymin=lower__), width = 0.4) +
  theme(axis.title.x = element_blank()) +
  labs(y = "Kdry") +
  coord_cartesian(ylim = c(1.1, 1.28))

Kdry.2

# HSI

ce1s_1 <- conditional_effects(cod_HSI_annual, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

plot.HSI.2 <- ce1s_1$year_fac

HSI.2 <- ggplot(plot.HSI.2, aes(year_fac, estimate__)) +
  geom_bar(stat = "identity", position = "dodge", fill = "grey80") +
  geom_errorbar(aes(ymax=upper__, ymin=lower__), width = 0.4) +
  theme(axis.title.x = element_blank()) +
  labs(y = "HSI") +
  coord_cartesian(ylim = c(3, 7))

HSI.2

# now combine into one plot
library(ggpubr)

png("./figs/four_panel_annual_data.png", 5, 8, units = 'in', res = 300)
ggarrange(abundance.2,
                  length.2,
                  ggarrange(Kdry.2,
                             HSI.2,
                             nrow = 1,
                             ncol = 2,
                             labels = c("c", "d")),
          nrow = 3,
          ncol = 1,
          labels = c("a", "b", ""))

dev.off()
