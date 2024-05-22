library(tidyverse)


df <- read_csv("/Users/julianashwin/Documents/GitHub/EcoNNet.jl/figures/pw_linear/sim_data/export_data.csv")

df %>%
  mutate(y_rounded = plyr::round_any(y, 0.5)) %>%
  group_by(y_rounded) %>%
  summarise(Epi_lead_mean = mean(Epi_lead)) %>%
  filter(abs(y_rounded) < 4) %>%
  print(n = 100)


df <- read_csv("/Users/julianashwin/Documents/GitHub/EcoNNet.jl/figures/pw_linear/sim_data/many_sims.csv")

df %>%
  filter(version != 0) %>%
  filter(progress <= 100) %>%
  ggplot() + theme_bw() + 
  geom_line(aes(x = progress, y = log(Rsq_Epi), group = version, color = version))



df %>%
  filter(version != 0) %>%
  group_by(version) %>%
  filter(progress == max(progress)) %>%
  ggplot() + theme_bw() + 
  geom_point(aes(x = y_lag, y = Epi_lead, group = version, color = factor(version)))

