library('readxl') # for data reading
library('dplyr')  # for data wrangling
library('tidyr')  # for data tidying
library('mgcv')   # for dating
source('analysis/default-figure-styling.R')

# import data
isotopes <-
  readRDS('data/isotopes.rds') %>%
  filter(year.scam > 1800) %>%
  pivot_longer(d15N:percC, names_to = 'name', values_to = 'value') %>%
  mutate(lab = case_when(name == 'c.n.ratio' ~ 'C:N~ratio',
                         name == 'd15N' ~ 'delta^{15}~N~\'\211\'',
                         name == 'percN' ~ 'N~\'\045\'',
                         name == 'd13C' ~ 'delta^{13}~C~\'\211\'',
                         name == 'percC' ~ 'C~\'\045\'') %>%
           factor(levels = c(
                             'C~\'\045\'', 'N~\'\045\'', 'C:N~ratio',
                             'delta^{13}~C~\'\211\'', 'delta^{15}~N~\'\211\'')))

ghosts <- expand_grid(x = 1800,
                      y = 0,
                      lab = c('N~\'\045\'', 'C~\'\045\'') %>%
                        factor(levels = levels(isotopes$lab)),
                      lake = unique(isotopes$lake))

p.isotopes <-
  ggplot(isotopes) +
  facet_grid(lab ~ lake, scales = 'free_y', labeller = label_parsed, switch = 'y') +
  geom_point(aes(year.scam, value), alpha = 0.5) +
  geom_point(aes(x, y), ghosts, alpha = 0) +
  labs(x = 'Year C.E.', y = NULL) +
  theme(strip.background = element_blank(), strip.placement = 'outside')
#p2pdf('isotopes.pdf', p.isotopes, x.plots = 1.5, y.plots = 2, device = 'pdf')
