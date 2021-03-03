library('mgcv')
library('scam')
library('ggplot2')
library('readxl')
library('tidyr')
library('cowplot')
library('dplyr')
library('gratia')
theme_set(theme_bw())
YEAR.MAX <- lubridate::decimal_date(as.POSIXlt('2016-07-18')) # coring date

# dating data (using midpoint depths)
dates <- bind_rows(read_xlsx(path = 'data/Kenosee2016-dates.xlsx',
                             sheet = 'age calculation CRS B',
                             range = c('X15:AC26'),
                             col_names = c('mid.depth', 'age', 'uncorrel.se',
                                           'correl.se', 'empty', 'year')) %>%
                     mutate(lake = 'Kenosee'),
                   read_xlsx(path = 'data/WhiteBear2016-dates.xlsx',
                             sheet = 'age calculation CRS B',
                             range = c('X15:AC26'),
                             col_names = c('mid.depth', 'age', 'uncorrel.se',
                                           'correl.se', 'empty', 'year')) %>%
                     mutate(lake = 'White~Bear')) %>%
  select(mid.depth, correl.se, year, lake) %>%
  filter(!is.na(year)) %>%
  mutate(weight = 1 / correl.se,
         weight = weight / mean(weight),
         lake = factor(lake))

# keep k low because SE is as big as Pb measurement in lower depths
m.age <- scam(year ~ s(mid.depth, k = 5, bs = 'mpd', by = lake),
              data = dates, family = gaussian(), weights = weight)

new.ages <- expand_grid(mid.depth = seq(0, 40, length.out = 400),
                        lake = c('Kenosee', 'White~Bear'))
new.ages$fitted <- predict(m.age, newdata = new.ages, type = 'response')

# plot with age models
ggplot() +
  geom_point(aes(year, mid.depth, col = lake, shape = lake), dates) +
  geom_line(aes(fitted, mid.depth, col = lake), new.ages) +
  scale_color_manual('Lake', values = c('forestgreen', 'goldenrod'),
                     labels = c('Kenosee', 'White Bear')) +
  scale_shape_manual('Lake', values = c(19, 4),
                     labels = c('Kenosee', 'White Bear')) +
  scale_y_reverse() +
  labs(x = 'Year CE', y = 'Depth (cm)') +
  theme(legend.position = 'top')

# add dates and temporal weights to data ####
pigments <-
  read_xlsx('data/KWB-pigments.xlsx') %>%
  select(Lake, Depth, Allo, `B-car`, Cantha, Diato, Echine, Lut_Zea, Myxo,
         Okenone, Phaeo_B, Pheo_A) %>%
  rename(lake = Lake, mid.depth = Depth, `beta~car` = `B-car`, Canth = Cantha,
         Echin = Echine, `Lut~Zea` = Lut_Zea, Oken = Okenone, `Pheo~A` = Pheo_A,
         `Pheo~B` = Phaeo_B) %>%
  mutate(lake = if_else(lake == 'WB', 'White~Bear', lake))
sum(is.na(pigments)) # check for NAs
pigments <-
  mutate(pigments,
         year.scam = predict(m.age, pigments, type = 'response'),
         year.scam = if_else(year.scam > YEAR.MAX, YEAR.MAX, year.scam)) %>%
  group_by(lake) %>%
  mutate(interval = if_else(lake == lag(lake),
                            true = lag(year.scam) - year.scam,
                            false = NA_real_),
         mean.interval = mean(interval, na.rm = TRUE),
         min.interval = min(interval, na.rm = TRUE),
         interval = if_else(is.na(interval), min.interval, interval),
         weight = interval / mean.interval) %>%
  select(- mean.interval, - min.interval) %>%
  ungroup()

# estimated weights are reasonable
ggplot() +
  facet_grid(lake ~ ., scales = 'free_y', labeller = label_parsed) +
  geom_point(aes(mid.depth, weight, color = mid.depth == 0), pigments) +
  theme(legend.position = 'top')

# switch to long format
pigments <-
  pivot_longer(pigments,
               names_to = 'pigment',
               values_to = 'conc',
               cols = -c('lake', 'mid.depth', 'weight', 'year.scam','interval')) %>%
  mutate(pigment = factor(pigment),
         lake = factor(lake),
         pigment_lake = interaction(pigment, lake, drop = TRUE))

# check data
unique(pigments$pigment)
apply(pigments, 2, function(x) range(x, na.rm = TRUE))

# save tibble for use in other files
#saveRDS(pigments, 'data/kwb-data.rds')
