library('mgcv')
library('scam')
library('readxl')
library('tidyr')
library('dplyr')
library('gratia')
source('analysis/default-figure-styling.R')
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
         lake = factor(lake),
         lwr.1se = year - correl.se,
         upr.1se = year + correl.se)

# age estimates are based on CRS model => high k; monotone decreasing p spline
m.age <- scam(year ~ s(mid.depth, k = 9, bs = 'mpd', by = lake),
              data = dates, family = gaussian(), weights = weight)

new.ages <- expand_grid(mid.depth = seq(0, 30, length.out = 400),
                        lake = c('Kenosee', 'White~Bear'))
new.ages <- bind_cols(new.ages,
                      predict(m.age, newdata = new.ages, se.fit = TRUE)) %>%
  mutate(lwr = fit - se.fit,
         upr = fit + se.fit)

# plot with age models
COLNAMES <- read_xlsx('data/Kenosee2016-dates.xlsx',
                      sheet = 'age calculation CRS B',
                      range = 'A11:AC13',
                      col_names = FALSE) %>%
  apply(MARGIN = 2,
        FUN = function(x) paste0('col.', paste(x[!is.na(x)], collapse = ' ')))
t(t(COLNAMES))

read.dating.data <- function(lake) {
  read_xlsx(paste0('data/', lake, '2016-dates.xlsx'),
            sheet = 'age calculation CRS B',
            range = 'A15:AC26',
            col_names = FALSE) %>%
    transmute(lake = if_else(lake == 'Kenosee', lake, 'White~Bear'),
              mid.depth = ...3,
              pb.210 = ...6,
              se.pb.210 = ...7,
              cs.137 = ...12 / 60,    # from decays/minute/g to Bq/g
              se.cs.137 = ...13 / 60, # from decays/minute/g to Bq/g
              year = ...29,
              se.year = ...26, # uncorrelated age SE, `...27` is correlated
              pb.210.lwr = pb.210 - se.pb.210, # add 1 SE ranges
              pb.210.upr = pb.210 + se.pb.210,
              cs.137.lwr = cs.137 - se.cs.137,
              cs.137.upr = cs.137 + se.cs.137,
              year.lwr = year - se.year,
              year.upr = year + se.year)
}

# single-core plot ----
kwb.dating <- bind_rows(read.dating.data('Kenosee'),
                        read.dating.data('WhiteBear')) %>%
  filter(mid.depth <= 30)

# top x axis label
isotope.lab <-
  expression(atop(''^{210}~Pb~activity~(Bq~g^{-1}~dry~mass),
                  ''^{137}~Cs~activity~(10^{-2}~Bq~g^{-1}~dry~mass)))

# fixed coefficients for secondary axis
YEAR.A <- 1780
PB.COEF <- 10
CS.COEF <- 100 * PB.COEF

new.ages <- filter(new.ages, lwr > YEAR.A)

kwb.dating.tidy <-
  select(kwb.dating, lake, mid.depth, pb.210, cs.137) %>%
  pivot_longer(-c('mid.depth', 'lake'), values_to = 'est') %>%
  left_join(select(kwb.dating,lake, mid.depth, se.pb.210, se.cs.137) %>%
              pivot_longer(-c('mid.depth', 'lake'), values_to = 'se') %>%
              mutate(name = substr(name, nchar('se.') + 1, nchar(name))),
            by = c('mid.depth', 'name', 'lake')) %>%
  mutate(slope = case_when(name == 'pb.210' ~ PB.COEF,
                           name == 'cs.137' ~ CS.COEF),
         intercept = YEAR.A,
         lwr = (est - se) * slope + intercept,
         upr = (est + se) * slope + intercept,
         est = est * slope + intercept)

p.triple <-
  ggplot(kwb.dating.tidy) +
  facet_grid(lake ~ ., labeller = label_parsed) +
  
  # year line
  geom_ribbon(aes(xmin = lwr, xmax = upr, y = mid.depth), new.ages,
              alpha = 0.2) +
  geom_line(aes(fit, mid.depth), new.ages) +
  geom_point(aes(year, mid.depth), dates, size = 1) +
  geom_errorbarh(aes(xmin = lwr.1se, xmax = upr.1se, y = mid.depth), dates,
                 height = 1) +
  
  # points with SE
  geom_errorbarh(aes(xmin = lwr, xmax = upr, y = mid.depth,
                     color = factor(name, levels = c('pb.210','cs.137'))),
                 height = 1) +
  geom_point(aes(est, mid.depth,
                 shape = factor(name, levels = c('pb.210','cs.137')),
                 color = factor(name, levels = c('pb.210','cs.137')))) +
  
  # other
  labs(x = 'Year C. E.', y = 'Depth (cm)') +
  scale_y_reverse(name = 'Depth (cm)',
                  breaks = seq(0, 35, by = 2.5),
                  labels = c(0, '', 5, '', 10, '', 15, '', 20, '',
                             25, '', 30, '', 35),
                  limits = c(NA, 0)) +
  scale_x_continuous(sec.axis =
                       dup_axis(~ (. - YEAR.A) / PB.COEF,
                                name = isotope.lab,
                                breaks = 0:18 * 12.5 / PB.COEF,
                                labels = c(0, '', 25, '', 50, '', 75, '', 100,
                                           '', 125, '', 150, '', 175, '', 200,
                                           '', 225)),
                     breaks = seq(1775, 2000, by = 25),
                     labels = c('', 1800, '', 1850, '',
                                1900, '', 1950, '', 2000)) +
  scale_shape_manual(NULL, values = c(19, 1)) +
  scale_color_brewer(NULL, type = 'qual', palette = 6) +
  theme(legend.position = c(.3, 1.1),
        legend.text = element_blank(),
        legend.key.height = unit(0.275, 'in'),
        legend.spacing.x = unit(-1, 'in'))
p.triple
#p2pdf('dating-model.pdf', p = p.triple, scale = 2.5)

# add dates and temporal weights to pigment data ####
kwb <-
  read_xlsx('data/KWB-pigments.xlsx') %>%
  select(Lake, Depth, Allo, Diato, Lut_Zea, Cantha, Okenone, Fuco,
         Echine, Phaeo_B, Pheo_A, `B-car`, Chl_Pheo_a, Chl.Pheo, UV_index) %>%
  rename(lake = Lake, mid.depth = Depth, `beta~car` = `B-car`, Canth = Cantha,
         Echin = Echine, `Lut~Zea` = Lut_Zea, Oken = Okenone, `Pheo~A` = Pheo_A,
         `Pheo~B` = Phaeo_B) %>%
  mutate(lake = if_else(lake == 'WB', 'White~Bear', lake))
sum(is.na(kwb)) # check for NAs

kwb <-
  mutate(kwb,
         year.scam = predict(m.age, kwb, type = 'response'),
         year.scam = if_else(year.scam > YEAR.MAX, YEAR.MAX, year.scam)) %>%
  group_by(lake) %>%
  mutate(interval = if_else(lake == lag(lake),
                            true = lag(year.scam) - year.scam,
                            false = NA_real_,
                            missing = YEAR.MAX - year.scam),
         mean.interval = mean(interval, na.rm = TRUE),
         weight = interval / mean.interval) %>%
  select(- mean.interval) %>%
  ungroup()

# weights for top slice are reasonable
ggplot() +
  facet_grid(lake ~ ., scales = 'free_y', labeller = label_parsed) +
  geom_point(aes(mid.depth, interval, color = mid.depth == 0), kwb) +
  theme(legend.position = 'top')

# check data
apply(kwb, 2, function(x) range(x, na.rm = TRUE)) %>% t()

# save tibble for use in other files
#saveRDS(kwb, 'data/kwb-data.rds')

# add dates to isotope data ####
isotopes <-
  read_xlsx('data/KenoseeWB2016-isotopes.xlsx') %>%
  transmute(sample = Sample,
            lake = case_when(grepl('K', sample) ~ 'Kenosee',
                             grepl('WB', sample) ~ 'White~Bear'),
            mid.depth = `Depth (cm)`,
            mg = `Wt.[mg]`,
            d15N = d15NAIR,
            d13C = d13CVPDB,
            #mgN = mgN,
            #mgC = mgC,
            #c.n.ratio = `C/N`,
            percN =  `%N`,
            percC = `%C`)

isotopes <-
  mutate(isotopes,
         year.scam = predict(m.age, isotopes, type = 'response'),
         year.scam = if_else(year.scam > YEAR.MAX, YEAR.MAX, year.scam))

# save tibble for use in other files
#saveRDS(isotopes, 'data/isotopes.rds')
