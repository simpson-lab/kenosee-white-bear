library('mgcv')
library('ggplot2')
library('readxl')
library('tidyr')
library('cowplot')
library('dplyr')
library('gratia')
source('analysis/default-figure-styling.R')
NCORES <- 4 # number of cores

kenosee <-
  read_xlsx('data/KenoseeWB2016-max depths_ calc.xlsx',
            range = 'A1:J54') %>%
  transmute(year = Year,
            depth = `Kenosee (m)`,
            max = `Est WB max depth`,
            lake = 'Kenosee') %>%
  filter(!is.na(depth))
wb <-
  read_xlsx('data/WhiteBear2016-water_level_(Cullimore_and_Griffin).xlsx',
            range = 'A2:C4884', col_names = c('year', 'x', 'depth')) %>%
  select(-x) %>%
  mutate(year = floor(year)) %>%
  group_by(year) %>%
  summarise(depth = mean(depth)) %>%
  mutate(lake = 'White Bear')

# relative to 1964 level
kenosee$depth <- kenosee$depth - filter(kenosee, year == 1964)$depth
wb$depth <- wb$depth - filter(wb, year == 1964)$depth

p <- ggplot(water) +
  geom_hline(yintercept = 0, color = 'grey80') +
  geom_line(aes(year, depth, color = lake), wb) +
  geom_line(aes(year, depth, color = lake), kenosee) +
  labs(x = 'Year C.E.', y = 'Relative change in lake level (m)') +
  scale_color_brewer(NULL, type = 'qual', palette = 6) +
  theme(legend.position = 'top')
p

p2pdf('lake-levels.pdf', p = p, width = 3, height = 4.5, scale = 1.5)
