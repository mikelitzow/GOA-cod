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
## using the best model from the fish-FAR project

## Read in data --------------------------------------------
cod.data <- read.csv("data/cpue.data.csv", row.names = 1)
cod.data$cod <- cod.data$cod.age.0
cod.data$bay_fac <- as.factor(cod.data$bay)
cod.data$year_fac <- as.factor(cod.data$year)
cod.data$site_fac <- as.factor(cod.data$site)

## fit: zero-inflated --------------------------------------

# set the formula
recr_2_formula <-  bf(cod ~ s(julian, k = 3) + (1 | bay_fac/site_fac) + (year_fac),
                      zi ~ s(julian, k = 3) + (1 | bay_fac/site_fac) + (year_fac))


## Set model distributions
zinb <- zero_inflated_negbinomial(link = "log", link_shape = "log", link_zi = "logit")

## Set priors
recr_2_priors_zinb <- c(set_prior("normal(0, 3)", class = "b"),
                        set_prior("normal(0, 3)", class = "Intercept"),
                        set_prior("student_t(3, 0, 3)", class = "sds"),
                        set_prior("gamma(0.01, 0.01)", class = "shape"),
                        set_prior("normal(0, 3)", class = "b", dpar = "zi"),
                        set_prior("logistic(0, 1)", class = "Intercept", dpar = "zi"),
                        set_prior("student_t(3, 0, 3)", class = "sds", dpar = "zi"))


## model 3 for post-settlement paper
recr_2_zinb <- brm(recr_2_formula,
                        data = cod.data,
                        prior = recr_2_priors_zinb,
                        family = zinb,
                        cores = 4, chains = 4, iter = 3000,
                        save_pars = save_pars(all = TRUE),
                        control = list(adapt_delta = 0.999, max_treedepth = 12))
saveRDS(recr_2_zinb, file = "output/recr_2_zinb.rds")

recr_2_zinb  <- add_criterion(recr_2_zinb, "bayes_R2")
saveRDS(recr_2_zinb, file = "output/recr_2_zinb.rds")

recr_2_zinb <- readRDS("./output/recr_2_zinb.rds")
check_hmc_diagnostics(recr_2_zinb$fit)
neff_lowest(recr_2_zinb$fit)
rhat_highest(recr_2_zinb$fit)
summary(recr_2_zinb)
bayes_R2(recr_2_zinb)
plot(conditional_smooths(recr_2_zinb), ask = FALSE)
pdf("./figs/trace_recr_2_zinb.pdf", width = 6, height = 4)
trace_plot(recr_2_zinb$fit)
dev.off()


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
                control = list(adapt_delta = 0.99, max_treedepth = 14))
saveRDS(cod_length_annual, file = "output/cod_length_annual.rds")

cod_length_annual  <- add_criterion(cod_length_annual, "bayes_R2", moment_match = TRUE)

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

# define formula
Kdry_formula <-  bf(Kdry ~ s(julian, k = 3) + (length) +  (year_fac) + (1 | bay_fac/site_fac))

## Show default priors
get_prior(Kdry_formula, Kdry.data)

# set priors

priors_Kdry <- c(set_prior("normal(0, 3)", class = "b"),
                set_prior("normal(0, 3)", class = "Intercept"),
                set_prior("student_t(3, 0, 2.5)", class = "sd"),
                set_prior("student_t(3, 0, 2.5)", class = "sds"),
                set_prior("student_t(3, 0, 2.5)", class = "sigma"))

## fit: brms --------------------------------------
cod_Kdry_annual <- brm(Kdry_formula,
                         data = Kdry.data,
                         prior = priors_Kdry,
                         cores = 4, chains = 4, iter = 3500,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.999999, max_treedepth = 12))
saveRDS(cod_Kdry_annual, file = "output/cod_Kdry_annual.rds")

cod_Kdry_annual  <- add_criterion(cod_Kdry_annual, "bayes_R2")
saveRDS(cod_Kdry_annual, file = "output/cod_Kdry_annual.rds")

cod_Kdry_annual <- readRDS("./output/cod_Kdry_annual.rds")
check_hmc_diagnostics(cod_Kdry_annual$fit) # 2 divergent transitions
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

## getting divergent transitions; trying simpler model
Kdry_formula2 <-  bf(Kdry ~ (julian) + (length) +  (year_fac) + (1 | bay_fac/site_fac))

## Show default priors
get_prior(Kdry_formula2, Kdry.data)

# set priors
priors_Kdry2 <- c(set_prior("normal(0, 3)", class = "b"),
                 set_prior("normal(0, 3)", class = "Intercept"),
                 set_prior("student_t(3, 0, 2.5)", class = "sd"),
                 set_prior("student_t(3, 0, 2.5)", class = "sigma"))

cod_Kdry_annual2 <- brm(Kdry_formula2,
                       data = Kdry.data,
                       prior = priors_Kdry2,
                       cores = 4, chains = 4, iter = 3500,
                       save_pars = save_pars(all = TRUE),
                       control = list(adapt_delta = 0.99999, max_treedepth = 12))
saveRDS(cod_Kdry_annual2, file = "output/cod_Kdry_annual2.rds")

cod_Kdry_annual2  <- add_criterion(cod_Kdry_annual2, "bayes_R2")
saveRDS(cod_Kdry_annual2, file = "output/cod_Kdry_annual2.rds")

cod_Kdry_annual2 <- readRDS("./output/cod_Kdry_annual2.rds")
check_hmc_diagnostics(cod_Kdry_annual2$fit)
neff_lowest(cod_Kdry_annual2$fit)
rhat_highest(cod_Kdry_annual2$fit)
summary(cod_Kdry_annual2)
bayes_R2(cod_Kdry_annual2)
plot(cod_Kdry_annual2$criteria$loo, "k")
plot(conditional_smooths(cod_Kdry_annual2), ask = FALSE)
y <- cod.length.data$Kdry
yrep_cod_Kdry_annual  <- fitted(cod_Kdry_annual2, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_Kdry_annual2[sample(nrow(yrep_cod_Kdry_annual2), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_Kdry_annual2")


# plot
## 95% CI
cod_Kdry_annual2 <- readRDS("./output/cod_Kdry_annual2.rds")

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


## four-panel figure in original units -----------------

# re-organize all for outputs

# first length
ce1s_1 <- conditional_effects(recr_2_zinb, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

plot.abundance <- ce1s_1$year_fac

plot.abundance$year <- as.numeric(as.character(plot.abundance$year_fac))

abundance <- ggplot(plot.abundance, aes(year, estimate__)) +
  geom_bar(stat = "identity", position = "dodge", fill = "grey80") +
  geom_errorbar(aes(ymax=upper__, ymin=lower__), width = 0.4) +
  coord_trans(y = "pseudo_log") +
  scale_y_continuous(breaks = c(0,1,5,10,50,100,200,500)) +
  theme(axis.title.x = element_blank()) +
  labs(y = "# fish / set") +
  geom_hline(yintercept = mean(plot.abundance$estimate__), lty = 2)

abundance

# length
ce1s_1 <- conditional_effects(cod_length_annual, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

plot.length <- ce1s_1$year_fac

plot.length$year <- as.numeric(as.character(plot.length$year_fac))

plot.length.2 <- data.frame(year = 2006:2020)

plot.length <- left_join(plot.length.2, plot.length)

# add offset to get to actual units!
plot.length$estimate__ <- plot.length$estimate__ + plot.length$length
plot.length$lower__ <- plot.length$lower__ + plot.length$length
plot.length$upper__ <- plot.length$upper__ + plot.length$length


length <- ggplot(plot.length, aes(year, estimate__)) +
  geom_bar(stat = "identity", position = "dodge", fill = "grey80") +
  geom_errorbar(aes(ymax=upper__, ymin=lower__), width = 0.4) +
  theme(axis.title.x = element_blank()) +
  labs(y = "Total length (mm)") +
  geom_hline(yintercept = mean(plot.length$estimate__, na.rm = T), lty = 2) +
  coord_cartesian(ylim = c(40,85))

length

# Kdry

ce1s_1 <- conditional_effects(cod_Kdry_annual2, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

plot.Kdry <- ce1s_1$year_fac

Kdry <- ggplot(plot.Kdry, aes(year_fac, estimate__)) +
  geom_bar(stat = "identity", position = "dodge", fill = "grey80") +
  geom_errorbar(aes(ymax=upper__, ymin=lower__), width = 0.4) +
  theme(axis.title.x = element_blank()) +
  labs(y = "Kdry") +
  coord_cartesian(ylim = c(1.14, 1.26))

Kdry

# HSI

ce1s_1 <- conditional_effects(cod_HSI_annual, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

plot.HSI <- ce1s_1$year_fac

HSI <- ggplot(plot.HSI, aes(year_fac, estimate__)) +
  geom_bar(stat = "identity", position = "dodge", fill = "grey80") +
  geom_errorbar(aes(ymax=upper__, ymin=lower__), width = 0.4) +
  theme(axis.title.x = element_blank()) +
  labs(y = "HSI") +
  coord_cartesian(ylim = c(3, 7))

HSI

# now combine into one plot
library(ggpubr)

png("./figs/Fig2_four_panel_annual_data.png", 5, 8, units = 'in', res = 300)
          plot <- ggarrange(abundance,
                  length,
                  ggarrange(Kdry,
                             HSI,
                             nrow = 1,
                             ncol = 2,
                             labels = c("c", "d")),
          nrow = 3,
          ncol = 1,
          labels = c("a", "b", "")) 
          
          annotate_figure(plot, bottom = text_grob("      Year"))

dev.off()
