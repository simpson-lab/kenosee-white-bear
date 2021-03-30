library('readxl') # for data reading
library('dplyr')  # for data wrangling
library('tidyr')  # for data tidying
library('mgcv')   # for dating
source('analysis/default-figure-styling.R')

# import data
isotopes <-
  readRDS('data/isotopes.rds') %>%
  pivot_longer(d15N:percC, names_to = 'name', values_to = 'value') %>%
  mutate(lab = case_when(name == 'd15N' ~ 'delta^{15}~N~\'\211\'',
                         name == 'percN' ~ 'N~\'\045\'',
                         name == 'd13C' ~ 'delta^{13}~C~\'\211\'',
                         name == 'percC' ~ 'C~\'\045\'') %>%
           factor(levels = c('delta^{15}~N~\'\211\'', 'delta^{13}~C~\'\211\'',
                             'N~\'\045\'', 'C~\'\045\'')))

p.isotopes <-
  ggplot(isotopes) +
  facet_grid(lab ~ lake, scales = 'free_y', labeller = label_parsed, switch = 'y') +
  geom_point(aes(mid.depth, value), alpha = 0.5) +
  labs(x = 'Year C.E.', y = NULL) +
  theme(strip.background = element_blank(), strip.placement = 'outside')
#p2pdf('isotopes.pdf', p.isotopes, x.plots = 1.5, y.plots = 2, device = 'pdf')
