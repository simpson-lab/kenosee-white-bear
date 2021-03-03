library('mgcv')
library('scam')
library('ggplot2')
library('readxl')
library('tidyr')
library('cowplot')
library('dplyr')
library('gratia')
theme_set(theme_bw())
NCORES <- 4 # number of logical processors/cores

# import data
pigments <-
  readRDS('data/kwb-data.rds') %>%
  select(-Chl_Pheo_a, -Chl.Pheo, -UV_index) %>%
  pivot_longer(names_to = 'pigment',
               values_to = 'conc',
               cols = -c('lake', 'mid.depth', 'weight', 'year.scam',
                         'interval')) %>%
  mutate(pigment = factor(pigment),
         lake = factor(lake),
         pigment_lake = interaction(pigment, lake, drop = TRUE))

# check data
unique(pigments$pigment)
apply(pigments, 2, function(x) range(x, na.rm = TRUE)) %>% t()

# pigment overview (different trends, different smoothness)
ggplot(pigments, aes(mid.depth, conc)) +
  facet_grid(pigment ~ lake, scales = 'free_y', labeller = label_parsed) +
  geom_point(aes(size = weight), alpha = 0.3) +
  geom_smooth(se = FALSE, method = 'gam', formula = y ~ s(x, k = 30)) +
  scale_x_reverse() +
  scale_size('Weight', range = c(0.75, 1.25)) +
  labs(x = 'Depth (cm, reversed)', y = 'Concentration')

# modelling ####
# by = pigment?
m.twlss <- readRDS('models/kwb-twlss.rds')
m.twlss <- gam(list(conc ~ s(year.scam, pigment_lake, k = 30, bs = 'fs'),
                    # constant mean-variance relationship
                    ~ 1,
                    # scale varies over time and by amount of info/slice
                    ~ s(year.scam, pigment_lake, k = 10, bs = 'fs') +
                      s(interval, k = 4, bs = 'cr')), # 3 gives \, 5 gives w
               data = pigments,
               family = twlss(),
               weights = weight,
               method = 'REML',
               control = gam.control(nthreads = NCORES))
# saveRDS(m.twlss, 'models/kwb-twlss.rds')

appraise(m.twlss) # diagnostics ok
summary(m.twlss)
draw(m.twlss)

# residuals
pigments$e <- pigments$conc - fitted(m.twlss)[,1]

# derivatives?
