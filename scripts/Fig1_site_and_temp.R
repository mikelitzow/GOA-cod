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
  scale_fill_manual(values=cb[c(3,7)]) +
  ggtitle(label = "a) Study site") +
  scale_x_continuous(breaks = c(-160, -156, -152)) +
  scale_y_continuous(breaks = c(55, 56, 57, 58))

map.plot

# now load temperature time series to plot
godas <- read.csv("./data/godas raw.csv", row.names = 1)

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

# alternate version of temp plot
alt.temp.plot <- ggplot(temps, aes(year, value, color=name)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=ymin.95, ymax=ymax.95)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=8)) +
  ylab("Temperature (°C)") +
  xlab("Year") +
  ggtitle(label = "b) Temperature") +
  scale_color_manual(values = cb[c(2,4,6)]) 

alt.temp.plot


# combine and save
png("figs/Fig1_study_site_pre_and_post_settlement_temp.png", 7, 3.5, units='in', res=300)
ggarrange(map.plot, alt.temp.plot,
          ncol = 2)
dev.off()
