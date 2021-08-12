library('ggplot2')
library('readxl')
library('dplyr')
library('mgcv')
library('gratia')

# strong collinearity between depths
depths <-
  left_join(read_xlsx('data/KenoseeWB2016-max depths_ calc.xlsx',
                      range = 'A1:J54') %>%
              transmute(year = Year, k = `Kenosee (m)`) %>%
              filter(!is.na(k)),
            read_xlsx(
              'data/WhiteBear2016-water_level_(Cullimore_and_Griffin).xlsx',
              range = 'A2:C4884', col_names = c('year', 'x', 'depth')) %>%
              select(-x) %>%
              mutate(year = floor(year)) %>%
              group_by(year) %>%
              summarise(w = mean(depth)),
            by = 'year')
m.depths <- gam(k ~ s(w), depths, family = Gamma('log'), method = 'REML')
summary(m.depths)
draw(m.depths, residuals = TRUE)

# isotopes and depth
i.w <- readRDS('data/isotopes.rds') %>%
  filter(lake == 'White~Bear') %>%
  mutate(year = round(year.scam),
         frac.c = percC / 100,
         frac.n = percN / 100) %>%
  left_join(select(depths, -k), by = 'year')

m.c.w <- gam(frac.c ~ s(w), i.w, family = betar(link = 'logit'), method = 'REML')
summary(m.c.w)
draw(m.c.w, residuals = TRUE)

m.n.w <- gam(frac.n ~ s(w), i.w, family = betar(link = 'logit'), method = 'REML')
summary(m.n.w)
draw(m.n.w, residuals = TRUE)

m.c.n.w <- gam(c.n.ratio ~ s(w), i.w, family = Gamma('log'), method = 'REML')
summary(m.c.n.w)
draw(m.c.n.w, residuals = TRUE)

# isotopes and pigments, UV
iso.pigm <- left_join(readRDS('data/isotopes.rds'),
                      readRDS('data/kwb-data.rds'),
                      by = c('mid.depth', 'lake', 'year.scam')) %>%
  mutate(lake = factor(lake))

m.n.o <- gam(d15N ~ s(Oken, by = lake), data = iso.pigm)
summary(m.n.o)
draw(m.n.o, residuals = TRUE)

m.n.c <- gam(d15N ~ s(Canth, by = lake), data = iso.pigm)
summary(m.n.c)
draw(m.n.c, residuals = TRUE)

# UV and depth
u.d <-
  readRDS('data/kwb-data.rds') %>%
  mutate(year = round(year.scam)) %>%
  left_join(depths %>%
              tidyr::pivot_longer(k:w) %>%
              transmute(year,
                        lake = if_else(name == 'k', 'Kenosee', 'White~Bear'),
                        depth = value), by = c('year', 'lake')) %>%
  mutate(lake = factor(lake))

m.n.u <- gam(UV_index ~ s(depth, by = lake), data = u.d)
summary(m.n.u)
draw(m.n.u, residuals = TRUE)

