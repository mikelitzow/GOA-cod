## Cod condition

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
cod.condition.data$fourth.root.cpue <- cod.condition.data$cod.cpue^0.25

# predict summer temp for each bay/year!

# read in temperature model
seine.temp <- readRDS("./output/temp3_gauss.rds") 
summary(seine.temp)

new.dat <- cod.condition.data %>%
  select(bay_fac, year_fac) %>%
  mutate(julian=mean(cod.condition.data$julian))

fitted.temp <- as.data.frame(fitted(seine.temp, newdata = new.dat))

cod.condition.data$fit.temp.mu <- fitted.temp$Estimate
cod.condition.data$fit.temp.sd <- fitted.temp$Est.Error

# check fitted mu
plot(cod.condition.data$fit.temp.mu)
hist(cod.condition.data$fit.temp.mu, breaks = 40) 

## Start with Kdry
## Check distributions
plot(cod.condition.data$Kdry)
hist(cod.condition.data$Kdry, breaks = 100) 
tab <- table(cod.condition.data$Kdry)
plot(tab)
summary(stats::glm(Kdry ~ 1, data = cod.condition.data, family = gaussian))

g <- ggplot(cod.condition.data) +
  aes(x = date, y = Kdry, color = site) +
  geom_point() +
  facet_wrap( ~ bay) +
  theme(legend.position = "none")
print(g)

g <- ggplot(cod.condition.data) +
  aes(x = fit.temp.mu, y = Kdry) +
  geom_point()
print(g)

g <- ggplot(cod.condition.data) +
  aes(x = julian, y = Kdry) +
  geom_point()
print(g)

g <- ggplot(cod.condition.data) +
  aes(x = length, y = Kdry) +
  geom_point()
print(g)

g <- ggplot(cod.condition.data) +
  aes(x = fourth.root.cpue, y = Kdry) +
  geom_point()
print(g)


## mgcv fits -----------------------------------------------
gam.Kdry.0 <- gam(Kdry ~ s(length, k = 4) + s(fit.temp.mu, k = 4),
                    data=cod.condition.data)
summary(gam.Kdry.0)
plot(gam.Kdry.0, pages = 1)

gam.Kdry.1 <- gam(Kdry ~ s(length, k = 4) + s(fit.temp.mu, k = 4) + s(fourth.root.cpue, k = 4),
                    data=cod.condition.data)
summary(gam.Kdry.1)
plot(gam.Kdry.1, pages = 1)

gam.Kdry.2 <- gam(Kdry ~ s(length, k = 4) + s(fit.temp.mu, k = 4) + s(fourth.root.cpue, k = 4) + bay_fac,
                  data=cod.condition.data)
summary(gam.Kdry.2)
plot(gam.Kdry.2, pages = 1)

MuMIn::AICc(gam.Kdry.0, gam.Kdry.1, gam.Kdry.2) # model 2


## brms: setup ---------------------------------------------
cod.condition.data$Kdry_stnd <- as.vector(scale(cod.condition.data$Kdry))

## Define model formulas
cond0_formula <- bf(Kdry_stnd ~ s(length, k = 4) + s(fit.temp.mu, k = 4))

cond1_formula <- bf(Kdry_stnd ~ s(length, k = 4) + s(fit.temp.mu, k = 4) +
                    s(fourth.root.cpue, k = 4))

cond2_formula <- bf(Kdry_stnd ~ s(length, k = 4) + s(fit.temp.mu, k = 4) +
                    s(fourth.root.cpue, k = 4) + bay_fac)


## Show default priors
get_prior(cond0_formula, cod.condition.data)

## Set priors
priors_kdry <- c(set_prior("normal(0, 3)", class = "b"),
                 set_prior("normal(0, 3)", class = "Intercept"),
                 set_prior("student_t(3, 0, 2)", class = "sds"),
                 set_prior("student_t(3, 0, 2)", class = "sigma"))



## fit: brms --------------------------------------
cod_cond0 <- brm(cond0_formula,
                 data = cod.condition.data,
                 cores = 4, chains = 4, iter = 3000,
                 save_pars = save_pars(all = TRUE),
                 prior = priors_kdry,
                 seed = 1234,
                 control = list(adapt_delta = 0.999, max_treedepth = 11))
cod_cond0  <- add_criterion(cod_cond0, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_cond0, file = "output/cod_cond0.rds")

cod_cond0 <- readRDS("./output/cod_cond0.rds")
check_hmc_diagnostics(cod_cond0$fit)
neff_lowest(cod_cond0$fit)
rhat_highest(cod_cond0$fit)
summary(cod_cond0)
bayes_R2(cod_cond0)
plot(cod_cond0$criteria$loo, "k")
plot(conditional_smooths(cod_cond0), ask = FALSE)
y <- cod.condition.data$Kdry_stnd
yrep_cod_cond0  <- fitted(cod_cond0, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_cond0[sample(nrow(yrep_cod_cond0), 25), ]) +
  xlim(0, 25) +
  ggtitle("cod_cond0")

## Debug divergent transitions
## Best to refit model without saving all pars
# trace_plot(cod_cond0$fit)
# posterior_cp <- as.array(cod_cond0$fit, pars = "lp__", include = FALSE)
# np_cp <- nuts_params(cod_cond0)
# color_scheme_set("darkgray")
# mcmc_parcoord(posterior_cp, np = np_cp)
# mcmc_pairs(posterior_cp, np = np_cp,
#            off_diag_args = list(size = 0.75))


cod_cond1 <- brm(cond1_formula,
                data = cod.condition.data,
                cores = 4, chains = 4, iter = 3000,
                save_pars = save_pars(all = TRUE),
                prior = priors_kdry,
                seed = 1234,
                control = list(adapt_delta = 0.9999, max_treedepth = 12))
cod_cond1  <- add_criterion(cod_cond1, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_cond1, file = "output/cod_cond1.rds")

cod_cond1 <- readRDS("./output/cod_cond1.rds")
check_hmc_diagnostics(cod_cond1$fit)
neff_lowest(cod_cond1$fit)
rhat_highest(cod_cond1$fit)
summary(cod_cond1)
bayes_R2(cod_cond1)
plot(cod_cond1$criteria$loo, "k")
plot(conditional_smooths(cod_cond1), ask = FALSE)
y <- cod.condition.data$Kdry_stnd
yrep_cod_cond1  <- fitted(cod_cond1, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_cond1[sample(nrow(yrep_cod_cond1), 25), ]) +
  xlim(0, 25) +
  ggtitle("cod_cond1")



cod_cond2 <- brm(cond2_formula,
                data = cod.condition.data,
                cores = 4, chains = 4, iter = 3000,
                save_pars = save_pars(all = TRUE),
                prior = priors_kdry,
                seed = 1234,
                control = list(adapt_delta = 0.9999, max_treedepth = 12))
cod_cond2  <- add_criterion(cod_cond2, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_cond2, file = "output/cod_cond2.rds")

cod_cond2 <- readRDS("./output/cod_cond2.rds")
check_hmc_diagnostics(cod_cond2$fit)
neff_lowest(cod_cond2$fit)
rhat_highest(cod_cond2$fit)
summary(cod_cond2)
bayes_R2(cod_cond2)
plot(cod_cond2$criteria$loo, "k")
plot(conditional_smooths(cod_cond2), ask = FALSE)
y <- cod.condition.data$Kdry_stnd
yrep_cod_cond2  <- fitted(cod_cond2, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_cond2[sample(nrow(yrep_cod_cond2), 25), ]) +
  xlim(0, 25) +
  ggtitle("cod_cond2")


# Model comparison
cod_cond0 <- readRDS("./output/cod_cond0.rds")
cod_cond1 <- readRDS("./output/cod_cond1.rds")
cod_cond2 <- readRDS("./output/cod_cond2.rds")

loo(cod_cond0, cod_cond1, cod_cond2) 

## save model comparison for ms.
model.comp <- loo(cod_cond0, cod_cond1, cod_cond2) 
model.comp

# save model comparison for ms.
forms <- data.frame(formula=c(as.character(cond2_formula)[1],
                              as.character(cond1_formula)[1],
                              as.character(cond0_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])
write.csv(comp.out, "./output/Kdry_model_comp.csv")



## SST predictions ##

## 95% CI
ce1s_1 <- conditional_effects(cod_cond2, effect = "fit.temp.mu", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod_cond2, effect = "fit.temp.mu", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod_cond2, effect = "fit.temp.mu", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$fit.temp.mu
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$fit.temp.mu[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$fit.temp.mu[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$fit.temp.mu[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$fit.temp.mu[["lower__"]]

g1 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1.5, color = "red3") +
  labs(x = "Summer temperature (ºC)", y = "Kdry (anomaly)") +
  theme_bw() +
  coord_cartesian(ylim = c(-1.1, 0.65))
print(g1)

ggsave("./figs/sst_predicted_effect_cod_cond2.png", width = 5, height = 4)


## 95% CI
ce1s_1 <- conditional_effects(cod_cond2, effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod_cond2, effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod_cond2, effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$fourth.root.cpue
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$fourth.root.cpue[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$fourth.root.cpue[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$fourth.root.cpue[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$fourth.root.cpue[["lower__"]]

g2 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1.5, color = "red3") +
  labs(x = "Fourth root CPUE", y = "Kdry (anomaly)") +
  coord_cartesian(ylim = c(-1.1, 0.65)) +
  theme_bw()

print(g2)
ggsave("./figs/fourth_root_cpue_predicted_effect_cod_cond2.png", width = 5, height = 4)

## and predict bay value!
ce1s_1 <- conditional_effects(cod_cond2, effect = "bay_fac", probs = c(0.025, 0.975))
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
mod.95$bay_fac <- as.factor(mod.95$bay_fac)

mod.95$long <- order$lon[match(mod.95$bay, order$Bay)]
mod.95$bay <- reorder(mod.95$bay, desc(mod.95$long))

mod.95$long <- order$lon[match(mod.95$bay_fac, order$Bay)]
mod.95$bay_fac <- reorder(mod.95$bay_fac, desc(mod.95$long))

theme_set(theme_bw())

g3 <- ggplot(mod.95) +
  aes(x = bay, y = estimate__) +
  geom_errorbar(aes(ymin = ymin.95, ymax = ymax.95), width = 0.5) +
  geom_point(size = 3) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=30, vjust=1, hjust=1)) +
  ylab("Kdry (anomaly)") +
  xlab("Bay") +
  coord_cartesian(ylim = c(-1.1, 0.65)) 

print(g3)

png("figs/Kdry_summer_temp.png", 10.5, 3, units='in', res=300)
ggpubr::ggarrange(g1, g2, g3,
          ncol = 3, labels = c("a", "b", "c"),
          widths = c(1,1,1.2))
dev.off()

## add scatter plot for SI -------------------------------

SI.dat <- cod.condition.data

# fix Anton and Cook
change <- grep("Cook", SI.dat$bay)
SI.dat$bay[change] <- "Cook"

SI.dat$jitter.Kdry <- jitter(SI.dat$Kdry, factor = 2)
SI.dat$jitter.fit.temp.mu <- jitter(SI.dat$fit.temp.mu, factor = 40)

g1 <- ggplot(SI.dat, aes(jitter.fit.temp.mu, jitter.Kdry)) +
  geom_point(size = 1, alpha = 0.5) +
  labs(x = "Summer temperature (ºC)", y = "Total length (mm)") +
  theme_bw() 

g1

SI.dat$jitter.fourth.root.cpue <- jitter(SI.dat$fourth.root.cpue, factor = 10)

g2 <- ggplot(SI.dat) +
  aes(jitter.fourth.root.cpue, jitter.Kdry) +
  geom_point(size = 1, alpha = 0.5) +
  labs(x = "Fourth root CPUE", y = "Kdry") +
  theme_bw() 

print(g2)

# set the bays to plot west - east
order <- read.csv("./data/bay_lat_long.csv")
order$Bay[1:2] <- c("Anton Larsen", "Cook")

SI.dat$long <- order$lon[match(SI.dat$bay, order$Bay)]
SI.dat$bay <- reorder(SI.dat$bay, desc(SI.dat$long))
SI.dat$bay.number <- as.numeric(SI.dat$bay)
SI.dat$jitter.bay <- jitter(SI.dat$bay.number, factor = 1)

g3 <- ggplot(SI.dat) +
  aes(jitter.bay, jitter.Kdry) +
  geom_point(size = 1, alpha = 0.5) +
  theme(axis.text.x = element_text(angle=30, vjust=1, hjust=1)) +
  coord_trans(y = "pseudo_log") +
  scale_x_continuous(breaks=c(1:14), minor_breaks = NULL, labels = levels(SI.dat$bay)) +
  ylab("Kdry") +
  xlab("Bay")

print(g3)

png("./figs/SI_Kdry_scatter_plots.png", 10.5, 3, units='in', res=300)
ggpubr::ggarrange(g1, g2, g3, ncol=3, nrow=1, labels=c("a", "b", "c"))
dev.off()

## HSI---------------------


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
cod.condition.data$fourth.root.cpue <- cod.condition.data$cod.cpue^0.25

# predict summer temp for each bay/year!

# read in temperature model
seine.temp <- readRDS("./output/temp3_gauss.rds") 
summary(seine.temp)

new.dat <- cod.condition.data %>%
  select(bay_fac, year_fac) %>%
  mutate(julian=mean(cod.condition.data$julian))

fitted.temp <- as.data.frame(fitted(seine.temp, newdata = new.dat))

cod.condition.data$fit.temp.mu <- fitted.temp$Estimate
cod.condition.data$fit.temp.sd <- fitted.temp$Est.Error

# check fitted mu
plot(cod.condition.data$fit.temp.mu)
hist(cod.condition.data$fit.temp.mu, breaks = 20) 

# compare HSI.dry vs Kdry 
# (raw, not controlled for length!)

# add length factor
cod.condition.data$length_fac <- if_else(cod.condition.data$length <=65, "<=65", ">=66")

hist(cod.condition.data$Kdry/cod.condition.data$length, breaks = 50)

ggplot(cod.condition.data, aes(length, Kdry/length, color=length_fac)) +
  geom_point(alpha=0.5, size = 1.5) +
  theme_bw()

ggplot(cod.condition.data, aes(HSI.dry, Kdry, color=length_fac)) +
  geom_point(alpha=0.5, size = 1.5) +
  geom_smooth(method="gam", se=F) + 
  theme_bw()

ggsave("./figs/HSI_v_Kdry_individually.png", width=6, height=4, units='in')

## Check distributions
sum(!is.na(cod.condition.data$Kdry))
sum(!is.na(cod.condition.data$HSI))

plot(cod.condition.data$HSI.dry)
hist(cod.condition.data$HSI.dry, breaks = 100) 
tab <- table(cod.condition.data$HSI.dry)
plot(tab)
summary(stats::glm(HSI.dry ~ 1, data = cod.condition.data, family = gaussian))

g <- ggplot(cod.condition.data) +
  aes(x = date, y = HSI.dry, color = site) +
  geom_point() +
  facet_wrap( ~ bay) +
  theme(legend.position = "none")
print(g)

g <- ggplot(cod.condition.data) +
  aes(x = fit.temp.mu, y = HSI.dry) +
  geom_point()
print(g)

g <- ggplot(cod.condition.data) +
  aes(x = julian, y = HSI.dry) +
  geom_point()
print(g)

g <- ggplot(cod.condition.data) +
  aes(x = length, y = HSI.dry) +
  geom_point()
print(g)

g <- ggplot(cod.condition.data) +
  aes(x = fourth.root.cpue, y = HSI.dry) +
  geom_point()
print(g)


## mgcv fits -----------------------------------------------
gam.HSI.dry.0 <- gam(HSI.dry ~ s(length, k = 4) + s(fit.temp.mu, k = 4),
                     data=cod.condition.data)
summary(gam.HSI.dry.0)
plot(gam.HSI.dry.0, pages = 1)

gam.HSI.dry.1 <- gam(HSI.dry ~ s(length, k = 4) + s(fit.temp.mu, k = 4) + s(fourth.root.cpue, k = 4),
                     data=cod.condition.data)
summary(gam.HSI.dry.1)
plot(gam.HSI.dry.1, pages = 1)

gam.HSI.dry.2 <- gam(HSI.dry ~ s(length, k = 4) + s(fit.temp.mu, k = 4) + s(fourth.root.cpue, k = 4) + bay_fac,
                     data=cod.condition.data)
summary(gam.HSI.dry.2)
plot(gam.HSI.dry.2, pages = 1)

MuMIn::AICc(gam.HSI.dry.0, gam.HSI.dry.1, gam.HSI.dry.2) # model 2


## brms: setup ---------------------------------------------
# scale HSI data
cod.condition.data$HSI.dry_stnd <- as.vector(scale(cod.condition.data$HSI.dry))

## Define model formulas
HSI0_formula <-  bf(HSI.dry_stnd ~ s(length, k = 4) + s(fit.temp.mu, k = 4))

HSI1_formula <-  bf(HSI.dry_stnd ~ s(length, k = 4) + s(fit.temp.mu, k = 4) + s(fourth.root.cpue, k = 4))

HSI2_formula <-  bf(HSI.dry_stnd ~ s(length, k = 4) + s(fit.temp.mu, k = 4) + s(fourth.root.cpue, k = 4) + bay_fac)

## Show default priors
get_prior(HSI0_formula, cod.condition.data)
# using defaults for now!

## fit: brms --------------------------------------
cod_HSI0 <- brm(HSI0_formula,
                data = cod.condition.data,
                cores = 4, chains = 4, iter = 3000,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.999, max_treedepth = 10))
cod_HSI0  <- add_criterion(cod_HSI0, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_HSI0, file = "output/cod_HSI0.rds")

cod_HSI0 <- readRDS("./output/cod_HSI0.rds")
check_hmc_diagnostics(cod_HSI0$fit)
neff_lowest(cod_HSI0$fit)
rhat_highest(cod_HSI0$fit)
summary(cod_HSI0)
bayes_R2(cod_HSI0)
plot(cod_HSI0$criteria$loo, "k")
plot(conditional_smooths(cod_HSI0), ask = FALSE)
y <- cod.condition.data $cod
yrep_cod_HSI0  <- fitted(cod_HSI0, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_HSI0[sample(nrow(yrep_cod_HSI0), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_HSI0")



cod_HSI1 <- brm(HSI1_formula,
                data = cod.condition.data,
                cores = 4, chains = 4, iter = 3000,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.999, max_treedepth = 10))
cod_HSI1  <- add_criterion(cod_HSI1, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_HSI1, file = "output/cod_HSI1.rds")

cod_HSI1 <- readRDS("./output/cod_HSI1.rds")
check_hmc_diagnostics(cod_HSI1$fit)
neff_lowest(cod_HSI1$fit)
rhat_highest(cod_HSI1$fit)
summary(cod_HSI1)
bayes_R2(cod_HSI1)
plot(cod_HSI1$criteria$loo, "k")
plot(conditional_smooths(cod_HSI1), ask = FALSE)
y <- cod.condition.data $HSI.dry_stnd
yrep_cod_HSI1  <- fitted(cod_HSI1, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_HSI1[sample(nrow(yrep_cod_HSI1), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_HSI1") 



cod_HSI2 <- brm(HSI2_formula,
                data = cod.condition.data,
                cores = 4, chains = 4, iter = 3000,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.999, max_treedepth = 10))
cod_HSI2  <- add_criterion(cod_HSI2, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_HSI2, file = "output/cod_HSI2.rds")

cod_HSI2 <- readRDS("./output/cod_HSI2.rds")
check_hmc_diagnostics(cod_HSI2$fit)
neff_lowest(cod_HSI2$fit)
rhat_highest(cod_HSI2$fit)
summary(cod_HSI2)
bayes_R2(cod_HSI2)
plot(cod_HSI2$criteria$loo, "k")
plot(conditional_smooths(cod_HSI2), ask = FALSE)
y <- cod.condition.data $HSI.dry_stnd
yrep_cod_HSI2  <- fitted(cod_HSI2, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_HSI2[sample(nrow(yrep_cod_HSI2), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_HSI2")


# Model comparison
cod_HSI0 <- readRDS("./output/cod_HSI0.rds")
cod_HSI1 <- readRDS("./output/cod_HSI1.rds")
cod_HSI2 <- readRDS("./output/cod_HSI2.rds")

model.comp <- loo(cod_HSI0, cod_HSI1, cod_HSI2) 
model.comp

# save model comparison for ms.
forms <- data.frame(formula=c(as.character(HSI2_formula)[1],
                              as.character(HSI1_formula)[1],
                              as.character(HSI0_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])
write.csv(comp.out, "./output/HSI_model_comp.csv")

## SST predictions ##

## 95% CI
ce1s_1 <- conditional_effects(cod_HSI2, effect = "fit.temp.mu", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod_HSI2, effect = "fit.temp.mu", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod_HSI2, effect = "fit.temp.mu", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$fit.temp.mu
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$fit.temp.mu[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$fit.temp.mu[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$fit.temp.mu[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$fit.temp.mu[["lower__"]]

g1 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1.5, color = "red3") +
  labs(x = "ºC", y = "HSI.dry_stnd") +
  theme_bw() +
  coord_cartesian(ylim = c(0,9))
print(g1)

ggsave("./figs/sst_predicted_effect_cod_HSI2.png", width = 5, height = 4)


## 95% CI
ce1s_1 <- conditional_effects(cod_HSI2, effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod_HSI2, effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod_HSI2, effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$fourth.root.cpue
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$fourth.root.cpue[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$fourth.root.cpue[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$fourth.root.cpue[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$fourth.root.cpue[["lower__"]]

g2 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1.5, color = "red3") +
  labs(x = "Fourth root CPUE", y = "HSI.dry_stnd") +
  coord_cartesian(ylim = c(0,9)) +
  theme_bw()
print(g2)
ggsave("./figs/fourth_root_cpue_predicted_effect_cod_HSI2.png", width = 5, height = 4)

## and predict bay value!
ce1s_1 <- conditional_effects(cod_HSI2, effect = "bay_fac", probs = c(0.025, 0.975))
mod.95 <- ce1s_1$bay_fac %>%
  select(bay_fac, estimate__, lower__, upper__)
names(mod.95)[3:4] <- c("ymin.95", "ymax.95")

# set the bays to plot west - east
order <- read.csv("./data/bay_lat_long.csv")

# remove "Bay" from first two bay names
order$Bay[1:2] <- c("Anton Larson", "Cook")

mod.95$bay_fac <- as.character(mod.95$bay_fac)
change <- grep("Anton", mod.95$bay_fac)
mod.95$bay_fac[change] <- "Anton Larson"

change <- grep("Cook", mod.95$bay_fac)
mod.95$bay_fac[change] <- "Cook"


mod.95$long <- order$lon[match(mod.95$bay_fac, order$Bay)]
mod.95$bay_fac <- reorder(mod.95$bay_fac, desc(mod.95$long))

theme_set(theme_bw())

g3 <- ggplot(mod.95) +
  aes(x = bay_fac, y = estimate__) +
  geom_errorbar(aes(ymin = ymin.95, ymax = ymax.95), width = 0.5) +
  geom_point(size = 3) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=30, vjust=1, hjust=1)) +
  ylab("HSI.dry_stnd") +
  coord_cartesian(ylim = c(0,9)) 


print(g3)

png("figs/HSIdry_summer_temp.png", 10.5, 3, units='in', res=300)
ggpubr::ggarrange(g1, g2, g3,
                  ncol = 3,
                  widths = c(1,1,1.2))
dev.off()
