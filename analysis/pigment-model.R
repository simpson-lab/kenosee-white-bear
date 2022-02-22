library('mgcv')
library('gratia')
library('ggplot2')
library('readxl')
library('tidyr')
library('cowplot')
library('dplyr')
source('analysis/default-figure-styling.R')
source('analysis/simulations-and-derivatives-functions.R')
NCORES <- 4 # number of cores

# import data
pigments <-
  readRDS('data/kwb-data.rds') %>%
  select(-Chl_Pheo_a, -Chl.Pheo, -UV_index) %>%
  pivot_longer(names_to = 'pigment',
               values_to = 'conc',
               cols = -c('lake', 'mid.depth', 'weight', 'year.scam',
                         'interval')) %>%
  mutate(conc = if_else(conc == 0, 0.5, conc), # set zeros to 1/2 mdl
         pigment = factor(pigment, levels = c('Fuco', 'Diato', 'Allo', 'Pheo~B',
                                              'Lut~Zea', 'Echin', 'Canth',
                                              'Oken', 'Pheo~A', 'beta~car')),
         lake = factor(lake),
         pigment_lake = interaction(pigment, lake, drop = TRUE)) %>%
  filter(year.scam > 1800) # start in 1800 to decreasing fitting time

# check data
unique(pigments$pigment)
apply(pigments, 2, function(x) range(x, na.rm = TRUE)) %>% t()

# pigment overview (different trends, different smoothness)
ggplot(pigments, aes(year.scam, conc)) +
  facet_grid(pigment ~ lake, scales = 'free_y', labeller = label_parsed) +
  geom_point(aes(size = weight), alpha = 0.3) +
  geom_smooth(se = FALSE, method = 'gam', formula = y ~ s(x, k = 30)) +
  scale_size('Weight', range = c(0.75, 1.75)) +
  labs(x = 'Estimated year C.E.', y = 'Concentration')

# use log(interval)
plot_grid(ggplot(pigments) + geom_point(aes(interval, conc)),
          ggplot(pigments) + geom_point(aes(log(interval), conc)),
          ncol = 1)

# modeling ####
# using wiggly age estimates: has a hard time fitting if k >= 10
# need pigment_lake because of different smoothness (Canth, Echin)
# adding global smooths decreases fitting time
m.gammals <- readRDS('models/kwb-gammals-wiggly-dating.rds')
if(FALSE) {
  tictoc::tic() # 5-6 minutes
  m.gammals <-
    gam(list(conc ~ s(year.scam, k = 10, bs ='cr') +
               s(year.scam, pigment_lake, k=10, bs='fs', xt = list(bs='cr')),
             ~ s(year.scam, k = 10, bs ='cr') +
               s(year.scam, pigment_lake, k = 5, bs = 'fs',xt = list(bs = 'cr'))+
               s(log(interval), k = 10, bs = 'ad')),
        data = pigments,
        family = gammals(),
        weights = weight,
        method = 'REML',
        control = gam.control(nthreads = NCORES, trace = FALSE))
  tictoc::toc()
  #beepr::beep()
  #saveRDS(m.gammals, 'models/kwb-gammals-wiggly-dating.rds')
  
  # check diagnostics
  appraise(m.gammals)
  draw(m.gammals)
  pigments <- mutate(pigments,
                     e = resid(m.gammals),
                     estar = abs(conc - m.gammals$fitted.values[, 1]))
  
  # residuals for Fuco are odd because of so many zeros
  # but it's not worth increasing k
  cowplot::plot_grid(ggplot(pigments, aes(year.scam, e, group = pigment)) +
                       facet_grid(lake ~ .) +
                       geom_hline(yintercept = 0, color = 'firebrick1') +
                       geom_point(),
                     ggplot(pigments, aes(log(interval), e, group = pigment)) +
                       geom_hline(yintercept = 0, color = 'firebrick1') +
                       facet_grid(lake ~ .) +
                       geom_point(),
                     ggplot(pigments, aes(log(conc), e, group = pigment)) +
                       facet_grid(pigment ~ lake, scale = 'free_x') +
                       geom_point(),
                     ggplot(pigments, aes(e, group = pigment)) +
                       facet_grid(lake ~ .) +
                       geom_density(aes(color = pigment == 'Fuco'),
                                    show.legend = FALSE))
}

# predictions ----
# create new data for regularly-spaced predictions
newd <- with(pigments,
             expand_grid(year.scam = seq(min(year.scam), max(year.scam), by= 1),
                         pigment = unique(pigment),
                         lake = unique(lake))) %>%
  left_join(group_by(pigments, lake) %>% summarise(interval = mean(interval)),
            by = 'lake') %>%
  mutate(pigment_lake = interaction(pigment, lake, drop = TRUE))

pred <- readRDS('data/predictions.rds')

if(FALSE) {
  # find derivatives
  set.seed(1)
  slopes.mu <-
    gammals_mean_deriv(m.gammals, data = newd, var = 'year.scam', nsims = 1e4)%>%
    group_by(year.scam, pigment_lake) %>%
    summarize(mu.deriv = median(derivative),
              lwr.mu.deriv = quantile(derivative, probs = 0.025),
              upr.mu.deriv = quantile(derivative, probs = 0.975),
              .groups = 'drop')
  slopes.s2 <-
    gammals_var_deriv(m.gammals, data = newd, var = 'year.scam', nsims = 1e4) %>%
    group_by(year.scam, pigment_lake) %>%
    summarize(s2.deriv = median(derivative),
              lwr.s2.deriv = quantile(derivative, probs = 0.025),
              upr.s2.deriv = quantile(derivative, probs = 0.975),
              .groups = 'drop')
  
  # predictions
  mu <- gammals_mean(model = m.gammals, data = newd, nsims = 1e4) %>%
    group_by(year.scam, pigment_lake) %>%
    summarize(mu = median(mean),
              lwr.mu = quantile(mean, probs = 0.025),
              upr.mu = quantile(mean, probs = 0.975),
              .groups = 'drop')
  s2 <- gammals_var(model = m.gammals, data = newd, nsims = 1e4) %>%
    group_by(year.scam, pigment_lake) %>%
    summarize(s2 = median(variance),
              lwr.s2 = quantile(variance, probs = 0.025),
              upr.s2 = quantile(variance, probs = 0.975),
              .groups = 'drop') %>%
    mutate(sd = sqrt(s2),
           lwr.sd = sqrt(lwr.s2),
           upr.sd = sqrt(upr.s2))
  
  pred <-
    newd %>%
    left_join(mu, by = c('year.scam', 'pigment_lake')) %>%
    left_join(s2, by = c('year.scam', 'pigment_lake')) %>%
    left_join(slopes.mu, by = c('year.scam', 'pigment_lake')) %>%
    left_join(slopes.s2, by = c('year.scam', 'pigment_lake')) %>%
    mutate(signif.mu = lwr.mu.deriv > 0 | upr.mu.deriv < 0,
           signif.s2 = lwr.s2.deriv > 0 | upr.s2.deriv < 0) %>%
    arrange(pigment_lake) %>%
    mutate(segm.mu.bool = signif.mu != lag(signif.mu) |
             pigment_lake != lag(pigment_lake),
           segm.mu = 1,
           segm.s2.bool = signif.s2 != lag(signif.s2) |
             pigment_lake != lag(pigment_lake),
           segm.s2 = 1)
  
  ## different groups for each significant segment
  for(v in paste0('segm.', c('mu', 's2'))) {
    for(i in 2:nrow(pred)) {
      if(pred[i, paste0(v, '.bool')][[1]]) {
        pred[i, v] <- pred[i - 1, v] + 1
      } else {
        pred[i, v] <- pred[i - 1, v]
      }
    }
  }
  
  #saveRDS(pred, 'data/predictions.rds')
}

# plot predictions ####
# mean
p.mu.k <- ggplot() +
  facet_grid(pigment ~ lake, scales = 'free', labeller = label_parsed) +
  geom_ribbon(aes(year.scam, ymin = lwr.mu, ymax = upr.mu),
              filter(pred, lake == 'Kenosee'), alpha = 0.3) +
  
  # lines
  geom_line(aes(year.scam, mu, group = segm.mu), color = 'red', # mean highlight
            filter(pred, lake == 'Kenosee', signif.mu), lwd = 2) +
  geom_line(aes(year.scam, mu), filter(pred, lake == 'Kenosee')) +
  
  # datapoints
  geom_point(aes(year.scam,conc), filter(pigments, lake=='Kenosee'),alpha = 0.3) +
  
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = 'Year C.E.', y = expression(Mean~concentration~(nmol~g^{-1}~C))) +
  theme(strip.text.y = element_blank())

p.mu.wb <- ggplot() +
  facet_grid(pigment ~ lake, scales = 'free', labeller = label_parsed) +
  geom_ribbon(aes(year.scam, ymin = lwr.mu, ymax = upr.mu),
              filter(pred, lake == 'White~Bear'), alpha = 0.3) +
  # lines
  geom_line(aes(year.scam, mu, group = segm.mu), color = 'red', # mean highlight
            filter(pred, lake == 'White~Bear', signif.mu), lwd = 2) +
  geom_line(aes(year.scam, mu), filter(pred, lake == 'White~Bear')) +
  
  # datapoints
  geom_point(aes(year.scam, conc), filter(pigments, lake=='White~Bear'), alpha=0.3) +
  
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = NULL) +
  theme(strip.text.y = element_blank())

# mean with variance highlight
p.mu.k2 <- ggplot() +
  facet_grid(pigment ~ lake, scales = 'free', labeller = label_parsed) +
  geom_ribbon(aes(year.scam, ymin = lwr.mu, ymax = upr.mu),
              filter(pred, lake == 'Kenosee'), alpha = 0.3) +
  
  # lines
  geom_line(aes(year.scam, mu, group = segm.s2), color = '#5E8BDE', # s2 highlight
            filter(pred, lake == 'Kenosee', signif.s2), lwd = 2.5) +
  geom_line(aes(year.scam, mu, group = segm.mu), color = 'red', # mean highlight
            filter(pred, lake == 'Kenosee', signif.mu), lwd = 2) +
  geom_line(aes(year.scam, mu), filter(pred, lake == 'Kenosee')) +
  
  # datapoints
  geom_point(aes(year.scam, conc), filter(pigments, lake == 'Kenosee'),alpha=0.3)+
  
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = 'Year C.E.', y = expression(Mean~concentration~(nmol~g^{-1}~C)))

p.mu.wb2 <- ggplot() +
  facet_grid(pigment ~ lake, scales = 'free', labeller = label_parsed) +
  geom_ribbon(aes(year.scam, ymin = lwr.mu, ymax = upr.mu),
              filter(pred, lake == 'White~Bear'), alpha = 0.3) +
  # lines
  geom_line(aes(year.scam, mu, group = segm.s2), color = '#5E8BDE', # s2 highlight
            filter(pred, lake == 'White~Bear', signif.s2), lwd = 2.5) +
  geom_line(aes(year.scam, mu, group = segm.mu), color = 'red', # mean highlight
            filter(pred, lake == 'White~Bear', signif.mu), lwd = 2) +
  geom_line(aes(year.scam, mu), filter(pred, lake == 'White~Bear')) +
  
  # datapoints
  geom_point(aes(year.scam, conc), filter(pigments, lake == 'White~Bear'), alpha = 0.3) +
  
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = NULL)

p.mus.nolabs <-
  plot_grid(p.mu.k2 + theme(axis.title.x = element_blank(),
                            strip.text.y = element_text(color = 'transparent')),
            p.mu.wb2,
            nrow = 1) %>%
  plot_grid(get_plot_component(p.mu.k, pattern = 'xlab-b'),
            nrow = 2,
            rel_heights = c(0.95, 0.05))
#p2pdf('mean-predictions-nolabs.pdf', p.mus.nolabs, scale = 2.5, y.plots = 2)

p.mus <- plot_grid(get_plot_component(p.mu.k, pattern = 'ylab-l'),
                   p.mu.k2 + theme(axis.title = element_blank(),
                                   strip.text.y = element_text(color = 'transparent')),
                   NULL,
                   p.mu.wb2,
                   rel_widths = c(0.15, 1, 0, 1),
                   nrow = 1) %>%
  plot_grid(get_plot_component(p.mu.k, pattern = 'xlab-b'),
            nrow = 2,
            rel_heights = c(0.95, 0.05)) + 
  draw_text(LABELS[1:20],
            x = sort(rep(c(.06, 0.53), 10)),
            y = rep(seq(.95, by = -0.089, length.out = 10), 2),
            family = 'serif')

#p2pdf('mean-predictions.pdf', p.mus, scale = 2.5, y.plots = 2)

# variance
# CIs cropped to 2.5 times the max estimated variance
K <- 2
s2.lab <- expression(paste(Concentration~variance~(nmol^2~g^{-2}~C)))
p.s2.k <-
  filter(pred, lake == 'Kenosee') %>%
  group_by(pigment) %>%
  mutate(upr.s2 = if_else(upr.s2 <= max(s2) * K, upr.s2, max(s2) * K)) %>%
  ungroup() %>%
  ggplot() +
  facet_grid(pigment ~ lake, scales = 'free', labeller = label_parsed) +
  geom_ribbon(aes(year.scam, ymin = lwr.s2, ymax = upr.s2), alpha = 0.3) +
  geom_line(aes(year.scam, s2, group = segm.s2), color = '#5E8BDE', lwd = 2,
            filter(pred, lake == 'Kenosee', signif.s2)) +
  geom_line(aes(year.scam, s2)) +
  labs(x = NULL, y = s2.lab) +
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())

p.s2.wb <-
  filter(pred, lake == 'White~Bear') %>%
  group_by(pigment) %>%
  mutate(upr.s2 = if_else(upr.s2 <= max(s2) * K, upr.s2, max(s2) * K)) %>%
  ungroup() %>%
  ggplot() +
  facet_grid(pigment ~ lake, scales = 'free', labeller = label_parsed) +
  geom_ribbon(aes(year.scam, ymin = lwr.s2, ymax = upr.s2), alpha = 0.3) +
  geom_line(aes(year.scam, s2, group = segm.s2), color = '#5E8BDE', lwd = 2,
            filter(pred, lake == 'White~Bear', signif.s2)) +
  geom_line(aes(year.scam, s2)) +
  labs(x = NULL, y = NULL)

# full figure
x.axis <- scale_x_continuous(name = NULL,
                             breaks = c(1800, 1900, 2000))
y.axis <- scale_y_continuous(n.breaks = 4, limits = c(0, NA))
p.full <- 
  plot_grid(plot_grid(p.mu.k + x.axis + y.axis,
                      p.mu.wb + x.axis + y.axis,
                      p.s2.k + x.axis + y.axis + coord_cartesian(ylim=c(0, NA)),
                      p.s2.wb + x.axis + y.axis + coord_cartesian(ylim=c(0,NA)),
                      rel_widths = c(.89, .77, .95, .96),
                      nrow = 1),
            get_plot_component(p.mu.k, pattern = 'xlab-b'),
            nrow = 2,
            rel_heights = c(0.95, 0.05))
#p2pdf('mean-variance-predictions.pdf', p.full, scale = 2.5, width = 3)

# plot coefficient of variation
mutate(pred, cv = sqrt(s2) / mu) %>%
  ggplot() +
  facet_grid(pigment ~ ., scales = 'free_y', labeller = label_parsed) +
  geom_line(aes(year.scam, cv, color = lake), lwd = 1) +
  scale_color_brewer('Lake', type = 'qual', palette = 1,
                     labels = c('Kenosee', 'White Bear')) +
  labs(x = 'Year C.E.', y = 'Coefficient of variation') +
  theme(legend.position = 'top')
 
