## Cod analysis
## Negative binomial
## Final set of models used in draft ms.

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
cod.data$cod <- cod.data$cod.age.0
cod.data$bay_fac <- as.factor(cod.data$bay)
cod.data$year_fac <- as.factor(cod.data$year)
cod.data$site_fac <- as.factor(cod.data$site)
cod.data$bay_site_fac <- as.factor(paste0(cod.data$bay, "_", cod.data$site))
cod.data$present <- ifelse(cod.data$cod > 0, 1, 0)
cod.data$date <- as.Date(cod.data$julian,
                         origin = paste0(cod.data$year, "-01-01"))

# add ssb data
ssb.data <- read.csv("./data/cod_pollock_assessment_2020_SAFEs.csv")
ssb.data <- ssb.data %>%
    select(year, codSSB.2020)
names(ssb.data)[2] <- "ssb"
cod.data <- left_join(cod.data, ssb.data)

temp.data <- read.csv("data/godas anomalies.csv")
cod.data$temp.anom <- temp.data$mean.anom[match(as.numeric(as.character(cod.data$year)), temp.data$year)]

## get sample sizes!
nrow(cod.data)
length(unique(cod.data$site))

site.bay <- cod.data %>%
    group_by(bay) %>%
    summarise(n=length(unique(site)))

site.bay

## brms: setup ---------------------------------------------

## Define model formula

# model with population-level bay effect for post-settlement paper

cod3s_sg_formula <-  bf(cod ~ s(julian, k = 3) + s(temp.anom, k = 3) + s(ssb, k = 3) + bay_fac+ (1 | bay_fac/site_fac),
                      zi ~ s(julian, k = 3) + s(temp.anom, k = 3) + s(ssb, k = 3) + bay_fac+ (1 | bay_fac/site_fac))


## Set model distributions
zinb <- zero_inflated_negbinomial(link = "log", link_shape = "log", link_zi = "logit")

## Set priors

mod3sg_priors_zinb <- c(set_prior("normal(0, 3)", class = "b"),
                    set_prior("normal(0, 3)", class = "Intercept"),
                    set_prior("student_t(3, 0, 3)", class = "sds"),
                    set_prior("gamma(0.01, 0.01)", class = "shape"),
                    set_prior("normal(0, 3)", class = "b", dpar = "zi"),
                    set_prior("logistic(0, 1)", class = "Intercept", dpar = "zi"),
                    set_prior("student_t(3, 0, 3)", class = "sds", dpar = "zi"))

## fit: zero-inflated --------------------------------------

## model 3 for post-settlement paper
cod3s_sg_zinb_k3 <- brm(cod3s_sg_formula,
                      data = cod.data,
                      prior = mod3sg_priors_zinb,
                      family = zinb,
                      cores = 4, chains = 4, iter = 3000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 12))
saveRDS(cod3s_sg_zinb_k3, file = "output/cod3s_sg_zinb_k3.rds")

cod3s_sg_zinb_k3  <- add_criterion(cod3s_sg_zinb_k3, "bayes_R2",
                                 moment_match = TRUE, reloo = TRUE,
                                 cores = 4, k_threshold = 0.7)
saveRDS(cod3s_sg_zinb_k3, file = "output/cod3s_sg_zinb_k3.rds")

cod3s_sg_zinb_k3 <- readRDS("./output/cod3s_sg_zinb_k3.rds")
check_hmc_diagnostics(cod3s_sg_zinb_k3$fit)
neff_lowest(cod3s_sg_zinb_k3$fit)
rhat_highest(cod3s_sg_zinb_k3$fit)
summary(cod3s_sg_zinb_k3)
bayes_R2(cod3s_sg_zinb_k3)
plot(conditional_smooths(cod3s_sg_zinb_k3), ask = FALSE)
pdf("./figs/trace_cod3s_sg_zinb_k3.pdf", width = 6, height = 4)
trace_plot(cod3s_sg_zinb_k3$fit)
dev.off()


## Plot predicted effects -------------------------------------

## temp.anom predictions ##

## 95% CI
ce1s_1 <- conditional_effects(cod3s_sg_zinb_k3, effect = "temp.anom", re_formula = NA,
                              probs = c(0.025, 0.975)) 
## 90% CI
ce1s_2 <- conditional_effects(cod3s_sg_zinb_k3, effect = "temp.anom", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod3s_sg_zinb_k3, effect = "temp.anom", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$temp.anom
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$temp.anom[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$temp.anom[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$temp.anom[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$temp.anom[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(cod.data$temp.anom), amount = 0.1),
                          rep(NA, 100-length(unique(cod.data$temp.anom))))


g1 <- ggplot(dat_ce) +
    aes(x = effect1__, y = estimate__) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
    geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
    geom_line(size = 1, color = "red3") +
    labs(x = "Egg/larval temperature anomaly", y = "# fish / set") +
    theme_bw() +
    coord_trans(y = "pseudo_log") +
    geom_rug(aes(x=rug.anom, y=NULL)) +
    scale_y_continuous(breaks=c(0, 1, 5, 10, 20, 50, 100, 200, 400, 600)) 
print(g1)

## Julian predictions ##

## 95% CI
ce1s_1 <- conditional_effects(cod3s_sg_zinb_k3, effect = "julian", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod3s_sg_zinb_k3, effect = "julian", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod3s_sg_zinb_k3, effect = "julian", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$julian

## aside - calculate % decline for post-settlement paper
prop_abund <- dat_ce %>%
    select(julian, estimate__)

prop_abund$percent_decline = 100*(max(prop_abund$estimate__)-prop_abund$estimate__)/max(prop_abund$estimate__)


write.csv(prop_abund, "./output/predicted seasonal abundance cod3s_sg_zinb_k3.csv")
#########################

dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$julian[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$julian[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$julian[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$julian[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(cod.data$julian), amount = 0.1),
                          rep(NA, 100-length(unique(cod.data$julian))))


g2 <- ggplot(dat_ce) +
    aes(x = effect1__, y = estimate__) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
    geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
    geom_line(size = 1.5, color = "red3") +
    labs(x = "Day of year", y = "# fish / set") +
    theme_bw() +
    coord_trans(y = "pseudo_log") +
    geom_rug(aes(x=rug.anom, y=NULL)) +
    scale_y_continuous(breaks=c(0, 1, 5, 10, 20, 50, 100, 200, 500, 1200)) 

print(g2)


## SSB predictions ##

## 95% CI
ce1s_1 <- conditional_effects(cod3s_sg_zinb_k3, effect = "ssb", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod3s_sg_zinb_k3, effect = "ssb", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod3s_sg_zinb_k3, effect = "ssb", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$ssb
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$ssb[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$ssb[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$ssb[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$ssb[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(cod.data$ssb), amount = 0.1),
                          rep(NA, 100-length(unique(cod.data$ssb))))

g3 <- ggplot(dat_ce) +
    aes(x = effect1__, y = estimate__) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
    geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
    geom_line(size = 1.5, color = "red3") +
    labs(x = "SSB (t)", y = "# fish / set") +
    theme_bw() +
    coord_trans(y = "pseudo_log") +
    geom_rug(aes(x=rug.anom, y=NULL)) +
    scale_y_continuous(breaks=c(0, 1, 5, 10, 20, 50, 100, 200, 600))

print(g3)


## and predict bay value!
ce1s_1 <- conditional_effects(cod3s_sg_zinb_k3, effect = "bay_fac", probs = c(0.025, 0.975))
mod.95 <- ce1s_1$bay_fac %>%
    select(bay_fac, estimate__, lower__, upper__)
names(mod.95)[3:4] <- c("ymin.95", "ymax.95")


# correct spelling of Anton Larsen; change Cook Bay to Cook
mod.95$bay_fac <- as.character(mod.95$bay_fac)

change <- grep("Anton", mod.95$bay_fac)
mod.95$bay_fac[change] <- "Anton Larsen"

change <- grep("Cook", mod.95$bay_fac)
mod.95$bay_fac[change] <- "Cook"

mod.95$bay_fac <- as.factor(mod.95$bay_fac)

# set the bays to plot west - east
order <- read.csv("./data/bay_lat_long.csv")
order$Bay[1:2] <- c("Anton Larsen", "Cook")

mod.95$long <- order$lon[match(mod.95$bay_fac, order$Bay)]
mod.95$bay_fac <- reorder(mod.95$bay_fac, desc(mod.95$long))



theme_set(theme_bw())

g4 <- ggplot(mod.95) +
    aes(x = bay_fac, y = estimate__) +
    geom_errorbar(aes(ymin = ymin.95, ymax = ymax.95), width = 0.5) +
    geom_point(size = 3) +
    theme(axis.text.x = element_text(angle=30, vjust=1, hjust=1)) +
    coord_trans(y = "pseudo_log") +
    scale_y_continuous(breaks=c(0, 1, 10, 100, 1000, 2000)) +
    ylab("# fish / set") +
    xlab("Bay")

print(g4)

## compare with k = 4 model---------------------------------------
cod3s_sg_formula_k4 <-  bf(cod ~ s(julian, k = 4) + s(temp.anom, k = 4) + s(ssb, k = 4) + bay_fac+ (1 | bay_fac/site_fac),
                        zi ~ s(julian, k = 4) + s(temp.anom, k = 4) + s(ssb, k = 4) + bay_fac+ (1 | bay_fac/site_fac))


## fit: zero-inflated --------------------------------------

## model 3 for post-settlement paper
cod3s_sg_zinb_k4 <- brm(cod3s_sg_formula_k4,
                        data = cod.data,
                        prior = mod3sg_priors_zinb,
                        family = zinb,
                        cores = 4, chains = 4, iter = 3000,
                        save_pars = save_pars(all = TRUE),
                        control = list(adapt_delta = 0.999, max_treedepth = 12))
saveRDS(cod3s_sg_zinb_k4, file = "output/cod3s_sg_zinb_k4.rds")

cod3s_sg_zinb_k4  <- add_criterion(cod3s_sg_zinb_k4, "bayes_R2",
                                   moment_match = TRUE, reloo = TRUE,
                                   cores = 4, k_threshold = 0.7)
saveRDS(cod3s_sg_zinb_k4, file = "output/cod3s_sg_zinb_k4.rds")

cod3s_sg_zinb_k4 <- readRDS("./output/cod3s_sg_zinb_k4.rds")
check_hmc_diagnostics(cod3s_sg_zinb_k4$fit)
neff_lowest(cod3s_sg_zinb_k4$fit)
rhat_highest(cod3s_sg_zinb_k4$fit)
summary(cod3s_sg_zinb_k4)
bayes_R2(cod3s_sg_zinb_k4)
plot(conditional_smooths(cod3s_sg_zinb_k4), ask = FALSE)
pdf("./figs/trace_cod3s_sg_zinb_k4.pdf", width = 6, height = 4)
trace_plot(cod3s_sg_zinb_k4$fit)
dev.off()



png("./figs/Fig3_combined_predicted_abundance_cod3s_sg_zinb_k3.png", 8, 6, units='in', res=300)
ggpubr::ggarrange(g1, g3, g2, g4, ncol=2, nrow=2, labels=c("a)", "b)", "c)", "d)"))
dev.off()

## add scatter plot for SI -------------------------------

SI.dat <- cod.data %>%
    select(bay, julian, cod.age.0, temp.anom, ssb)

# fix Anton and Cook
change <- grep("Anton", SI.dat$bay)
SI.dat$bay[change] <- "Anton Larsen"

change <- grep("Cook", SI.dat$bay)
SI.dat$bay[change] <- "Cook"

g1 <- ggplot(SI.dat, aes(temp.anom, cod.age.0)) +
    geom_point(size = 1, alpha = 0.5) +
    labs(x = "Egg/larval temperature anomaly", y = "# fish / set") +
    theme_bw() +
    coord_trans(y = "pseudo_log") +
    scale_y_continuous(breaks=c(0, 1, 5, 10, 20, 50, 100, 200, 500, 2000), minor_breaks = NULL)

g1


g2 <- ggplot(SI.dat) +
    aes(julian, cod.age.0) +
    geom_point(size = 1, alpha = 0.5) +
    labs(x = "Day of year", y = "# fish / set") +
    theme_bw() +
    coord_trans(y = "pseudo_log") +
    scale_y_continuous(breaks=c(0, 1, 5, 10, 20, 50, 100, 200, 500, 2000), minor_breaks = NULL)

print(g2)


g3 <- ggplot(SI.dat) +
    aes(ssb/1000, cod.age.0) +
    geom_point(size = 1, alpha = 0.5) +
    labs(x = "SSB (1000 t)", y = "# fish / set") +
    theme_bw() +
    coord_trans(y = "pseudo_log") +
    scale_y_continuous(breaks=c(0, 1, 5, 10, 20, 50, 100, 200, 500, 2000), minor_breaks = NULL) +
    scale_x_continuous(breaks=as.numeric(c(40, 60, 80, 100)))

print(g3)

# set the bays to plot west - east
order <- read.csv("./data/bay_lat_long.csv")
order$Bay[1:2] <- c("Anton Larsen", "Cook")

SI.dat$long <- order$lon[match(SI.dat$bay, order$Bay)]
SI.dat$bay <- reorder(SI.dat$bay, desc(SI.dat$long))

g4 <- ggplot(SI.dat) +
    aes(bay, cod.age.0) +
    geom_point(size = 1, alpha = 0.5) +
    theme(axis.text.x = element_text(angle=30, vjust=1, hjust=1)) +
    coord_trans(y = "pseudo_log") +
    scale_y_continuous(breaks=c(0, 1, 5, 10, 20, 50, 100, 200, 500, 2000), minor_breaks = NULL) +
    ylab("# fish / set") +
    xlab("Bay")

print(g4)

png("./figs/SI_abundance_scatter_plots.png", 8, 6, units='in', res=300)
ggpubr::ggarrange(g1, g3, g2, g4, ncol=2, nrow=2, labels=c("a)", "b)", "c)", "d)"))
dev.off()


## extra comparison model - 2018 / 2020 bay effects ------------------------------

# this analysis evaluates bay differences in abundance in the two years 
# that the entire western Gulf of Alaska study area was sampled *and* 
# cod abundance was high - the goal here was to confirm that our finding of
# modest differences in recruitment across bays was not driven by the ubiquitous
# recruitment failure in 2019

cod_bay1_formula <-  bf(cod ~ s(julian, k = 3) + bay_fac + year_fac,
                        zi ~ s(julian, k = 3) + bay_fac + year_fac)

cod_bay1s_formula <-  bf(cod ~ s(julian, k = 3) + bay_fac + year_fac + (1 | bay_fac/site_fac),
                        zi ~ s(julian, k = 3) + bay_fac + year_fac + (1 | bay_fac/site_fac))


## Set model distributions
zinb <- zero_inflated_negbinomial(link = "log", link_shape = "log", link_zi = "logit")


## Show default priors
# get_prior(cod0_formula, cod.data, family = zinb)
# get_prior(cod0_formula_nb, cod.data, family = nb)

## Set priors
priors_zinb_k3 <- c(set_prior("normal(0, 3)", class = "b"),
                    set_prior("normal(0, 3)", class = "Intercept"),
                    set_prior("student_t(3, 0, 3)", class = "sds"),
                    set_prior("gamma(0.01, 0.01)", class = "shape"),
                    set_prior("normal(0, 3)", class = "b", dpar = "zi"),
                    set_prior("logistic(0, 1)", class = "Intercept", dpar = "zi"),
                    set_prior("student_t(3, 0, 3)", class = "sds", dpar = "zi"))

cod_bay1 <- brm(cod_bay1_formula,
                        data = filter(cod.data, year %in% c(2018, 2020)),
                        prior = priors_zinb_k3,
                        family = zinb,
                        cores = 4, chains = 4, iter = 3000,
                        save_pars = save_pars(all = TRUE),
                        control = list(adapt_delta = 0.999, max_treedepth = 12))
saveRDS(cod_bay1, file = "output/cod_bay1.rds")

cod_bay1  <- add_criterion(cod_bay1, c("loo", "bayes_R2"),
                                   moment_match = TRUE, reloo = TRUE,
                                   cores = 4, k_threshold = 0.7)
saveRDS(cod_bay1, file = "output/cod_bay1.rds")

cod_bay1 <- readRDS("./output/cod_bay1.rds")
check_hmc_diagnostics(cod_bay1$fit)
neff_lowest(cod_bay1$fit)
rhat_highest(cod_bay1$fit)
summary(cod_bay1)
bayes_R2(cod_bay1)
plot(cod_bay1$criteria$loo, "k")
plot(conditional_smooths(cod_bay1), ask = FALSE)
pdf("./figs/trace_cod_bay1.pdf", width = 6, height = 4)
trace_plot(cod_bay_18.20$fit)
dev.off()


cod_bay1s <- brm(cod_bay1s_formula,
                data = filter(cod.data, year %in% c(2018, 2020)),
                prior = priors_zinb_k3,
                family = zinb,
                cores = 4, chains = 4, iter = 3000,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.999, max_treedepth = 12))
saveRDS(cod_bay1s, file = "output/cod_bay1s.rds")

cod_bay1s  <- add_criterion(cod_bay1s, c("loo", "bayes_R2"),
                           moment_match = TRUE, reloo = TRUE,
                           cores = 4, k_threshold = 0.7)
saveRDS(cod_bay1s, file = "output/cod_bay1s.rds")

cod_bay1s <- readRDS("./output/cod_bay1s.rds")
check_hmc_diagnostics(cod_bay1s$fit)
neff_lowest(cod_bay1s$fit)
rhat_highest(cod_bay1s$fit)
summary(cod_bay1s)
bayes_R2(cod_bay1s)
plot(cod_bay1s$criteria$loo, "k")
plot(conditional_smooths(cod_bay1s), ask = FALSE)
pdf("./figs/trace_cod_bay1s.pdf", width = 6, height = 4)
trace_plot(cod_bay_18.20$fit)
dev.off()

## model selection
cod_bay1 <- readRDS("./output/cod_bay1.rds")
cod_bay1s <- readRDS("./output/cod_bay1s.rds")
loo(cod_bay1, cod_bay1s)

## and predict bay value!
ce1s_1 <- conditional_effects(cod_bay1s, effect = "bay_fac", probs = c(0.025, 0.975))
mod.95 <- ce1s_1$bay_fac %>%
    select(bay_fac, estimate__, lower__, upper__)
names(mod.95)[3:4] <- c("ymin.95", "ymax.95")

# set the bays to plot west - east
order <- read.csv("./data/bay_lat_long.csv")
mod.95$long <- order$lon[match(mod.95$bay_fac, order$Bay)]
mod.95$bay_fac <- reorder(mod.95$bay_fac, desc(mod.95$long))

theme_set(theme_bw())

g4 <- ggplot(mod.95) +
    aes(x = long, y = estimate__) +
    geom_errorbar(aes(ymin = ymin.95, ymax = ymax.95), width = 0.5, color="dark grey") +
    geom_point(size = 2) +
    labs(y="Cod abundance", x="°W longitude") +
    coord_trans(y = "pseudo_log") +
    scale_y_continuous(breaks=c(0, 1, 10, 100, 1000, 2000))+
    scale_x_reverse()

print(g4)
