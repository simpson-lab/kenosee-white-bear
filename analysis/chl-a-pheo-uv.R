library('mgcv')
library('scam')
library('readxl')
library('tidyr')
library('dplyr')
library('gratia')
source('analysis/default-figure-styling.R')
NCORES <- 4 # number of logical processors/cores

# import data
d <-
  readRDS('data/kwb-data.rds') %>%
  select(lake, year.scam, mid.depth, Chl.Pheo, UV_index, weight) %>%
  mutate(lake = factor(lake)) %>%
  filter(year.scam >= 1800)

d %>%
  pivot_longer(c('Chl.Pheo', 'UV_index')) %>%
  ggplot() +
  facet_grid(name ~ lake, scales = 'free_y') +
  geom_point(aes(year.scam, value), alpha = 0.5) +
  labs(x = 'Estimated year CE', y = NULL) +
  ylim(c(0, NA))

# modeling ####
m.chl.pheo <- gam(Chl.Pheo ~ s(year.scam, by = lake, k = 15, bs = 'ad'),
                  family = Gamma(link = 'log'),
                  data = d,
                  method = 'REML',
                  weights = weight)

m.uv <- gam(UV_index ~ s(year.scam, by = lake, k = 10, bs = 'ad'),
            family = tw(link = 'log'), # some zeros (also, p = 1.028)
            data = d,
            method = 'REML',
            weights = weight)

appraise(m.chl.pheo, method = 'simulate')
appraise(m.uv)

# plots ###
pred <- expand_grid(year.scam = seq_min_max(d$year.scam, n = 400),
                    lake = unique(d$lake)) %>%
  arrange(lake)

pred <-
  bind_rows(bind_cols(pred,
                      predict.gam(m.chl.pheo, pred, se.fit = TRUE),
                      derivatives(m.chl.pheo, newdata = pred) %>%
                        filter(lower != upper)) %>% # REMOVE ###################
            mutate(y.lab = 'Chl:pheo~ratio'),
            bind_cols(pred,
                      predict.gam(m.uv, pred, se.fit = TRUE),
                      derivatives(m.uv, newdata = pred) %>%
                        filter(lower != upper)) %>% # REMOVE ###################
            mutate(y.lab = 'Estimated~UV~index')) %>%
  mutate(muhat = exp(fit),
         lwr = exp(fit - se.fit * 1.96),
         upr = exp(fit + se.fit * 1.96),
         significant = upper < 0 | lower > 0,
         group = 1)

for(i in 2:nrow(pred)) {
  current <- pred[i, ]
  previous <- pred[i - 1, ]
  
  # is the next observation in a different group as the previous one?
  k <- current$significant != previous$significant |
    current$lake != previous$lake | current$y.lab != previous$y.lab
  
  pred$group[i] <- pred$group[i - 1] + k
}

p <-
  d %>%
  pivot_longer(c('Chl.Pheo', 'UV_index'), names_to = 'parameter') %>%
  mutate(y.lab = case_when(parameter == 'Chl.Pheo' ~ 'Chl:pheo~ratio',
                           parameter == 'UV_index' ~ 'Estimated~UV~index')) %>%
  ggplot() +
  facet_grid(y.lab ~ lake, scales = 'free_y', labeller = label_parsed,
             switch = 'y') +
  geom_ribbon(aes(year.scam, ymin = lwr, ymax = upr), pred, alpha = 0.2) +
  geom_line(aes(year.scam, muhat, group = group), filter(pred, significant),
            lwd = 2, col = 'red') +
  geom_line(aes(year.scam, muhat), pred) +
  geom_point(aes(year.scam, value), alpha = 0.5) +
  labs(x = 'Estimated year C.E.', y = NULL) +
  theme(strip.placement = 'outside')
p

p2pdf('chl-a-pheo-uv.pdf', p, scale = 2)
