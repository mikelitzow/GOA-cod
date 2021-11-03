library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(brms)

theme_set(theme_bw())

# set palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load data sets
bays <- read.csv("./data/bay_lat_long.csv", row.names = 1)

# remove Cooks (repeat) and Kujulik (not resampled)
drop <- bays$Bay %in% c("Cooks", "Kujulik")
bays <- bays[!drop,]
bays$years <- ifelse(bays$Bay %in% c("Anton Larson Bay", "Cook Bay"), "2006-2020", "2018-2020")

# study site plot

ak <- ne_countries(scale = "large", returnclass = "sf", continent="north america")
world <- ne_countries(scale='medium', returnclass = "sf")

# set legend title
legend_title <- "Years sampled"
map.plot <- ggplot(ak) +
  geom_sf(fill="darkgoldenrod3", color=NA) + 
  coord_sf(xlim = c(-162, -151), ylim = c(54.5, 58.5), expand = FALSE) +
  geom_point(data = bays, aes(-lon, lat, fill=years), size=2, shape=21, color="black") +
  theme(axis.title = element_blank(),
        legend.position = c(0.8, 0.25),
        legend.text = element_text(size=8),
        legend.title = element_blank(),
        legend.box.margin = margin(1,1,1,1,unit='pt')) +
  scale_fill_manual(values=cb[c(4,7)]) +
  ggtitle(label = "a) Study site") +
  scale_x_continuous(breaks = c(-160, -156, -152)) +
  scale_y_continuous(breaks = c(55, 56, 57, 58))

map.plot

# now load temperature time series to plot
godas <- read.csv("./climate/godas raw.csv", row.names = 1)

names(godas)[2:3] <- c("Larval", "Egg")

godas <- godas %>%
  pivot_longer(cols = -year)

godas$ymin.95 <- godas$ymax.95 <- NA

summer <- readRDS("./output/temp3_gauss.rds")
## 95% CI
ce1s_1 <- conditional_effects(summer, probs = c(0.025, 0.975))
summer.95 <- ce1s_1$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(summer.95) <- c("year", "Post-settlement", "ymin.95", "ymax.95")
summer.95$year <- as.numeric(as.character(summer.95$year))
summer.95 <- summer.95 %>%
  pivot_longer(cols = c(-year, -ymin.95, -ymax.95))

temps <- rbind(godas, summer.95)

temps <- temps %>%
  filter(year >= 2006)

temps$order <- ifelse(temps$name=="Post-settlement", 1,
                      ifelse(temps$name=="Larval", 2, 3))

temps$name <- reorder(temps$name, temps$order)

temp.plot <- ggplot(temps, aes(year, value, color=name)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=ymin.95, ymax=ymax.95)) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.position = c(0.25, 0.95),
        legend.direction = "horizontal",
        legend.background = element_rect(color="black", size = 0.2)) +
  ylab("°C") +
  ggtitle(label = "b) Temperature") +
  scale_color_manual(values = cb[c(2,4,6)]) +
  ylim(3,13.5)

temp.plot

# alternate version of tempr plot
alt.temp.plot <- ggplot(temps, aes(year, value, color=name)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=ymin.95, ymax=ymax.95)) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=8)) +
  ylab("°C") +
  ggtitle(label = "b) Temperature") +
  scale_color_manual(values = cb[c(2,4,6)]) 
alt.temp.plot


# initial two-panel version
png("figs/study_site_pre_and_post_settlement_temp.png", 7, 3.5, units='in', res=300)
ggarrange(map.plot, alt.temp.plot,
          ncol = 2)
dev.off()


## and FAR
obs_far_fixef <- readRDS("./output/obs_far_fixef.rds")
## 95% CI
ce1s_1 <- conditional_effects(obs_far_fixef, probs = c(0.025, 0.975))
obs.95 <- ce1s_1$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(obs.95)[3:4] <- c("ymin.95", "ymax.95")


obs.95$year <- as.numeric(as.character(obs.95$year_fac))

obs.95 <- obs.95 %>%
  filter(year >= 2006)

far.plot <- ggplot(obs.95) +
  aes(x = year, y = estimate__) +
  geom_point(size=2) +
  geom_line() +
  geom_errorbar(aes(ymin=ymin.95, ymax=ymax.95)) +
  theme(axis.title.x = element_blank()) +
  ylab("FAR") +
  ggtitle(label = "c) Fraction of attributable risk")

print(far.plot)

png("figs/study_site_and_temp.png", 3.5, 7, units='in', res=300)
ggarrange(map.plot, temp.plot, far.plot,
          ncol = 1)
dev.off()

# version with only pre-settlement temps!
pre.temps <- temps %>%
  filter(name != "Juvenile")

pre.temps$order <- ifelse(pre.temps$name=="Egg", 1, 2)
pre.temps$name <- reorder(pre.temps$name, pre.temps$order)  
  
pre.temp.plot <- ggplot(pre.temps, aes(year, value, color=name)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=ymin.95, ymax=ymax.95)) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.position = c(0.25, 0.85),
        legend.direction = "horizontal",
        legend.background = element_rect(color="black", size = 0.2)) +
  ylab("°C") +
  ggtitle(label = "b) Temperature") +
  scale_color_manual(values = cb[c(2,4)]) 

pre.temp.plot

png("figs/study_site_and_pre_settlement_temp.png", 3.5, 7, units='in', res=300)
ggarrange(map.plot, pre.temp.plot, far.plot,
          ncol = 1)
dev.off()

# inset <- ggplot(world) +
#   geom_sf(fill="darkgoldenrod3", color=NA) + 
#   coord_sf(xlim = c(-180, -120), ylim = c(20, 66), expand = FALSE,
#            crs = "+proj=laea +lat_0=35 +lon_0=-100")
# inset
# 
# # ggplot(data = world) +
# #   geom_sf(fill="darkgoldenrod3", color=NA) +
# #   coord_sf(ylim=c(20, 68),
# #            xlim=c(120, -100),
# #     crs = "+proj=laea +lat_0=52 +lon_0=-160 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs ")
# 
# 
# 
# area.box<- data.frame(yy=c(54, 58.5, 58.5, 54, 54), xx=c(-180, -180, -120, -120, -180))
# ggdraw(map.plot) +
#   draw_plot(
#     { 
# 
#         geom_path(aes(xx, yy), area.box, color=cb[6]) +
#         theme(axis.title = element_blank(), axis.text = element_blank())
#     },
#     # The distance along a (0,1) x-axis to draw the left edge of the plot
#     x = 0.58, 
#     # The distance along a (0,1) y-axis to draw the bottom edge of the plot
#     y = 0,
#     # The width and height of the plot expressed as proportion of the entire ggdraw object
#     width = 0.36, 
#     height = 0.46)
