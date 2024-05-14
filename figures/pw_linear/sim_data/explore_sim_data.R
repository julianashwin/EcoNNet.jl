library(tidyverse)


df <- read_csv("/Users/julianashwin/Documents/GitHub/EcoNNet.jl/figures/pw_linear/sim_data/export_data.csv")

df %>%
  mutate(y_rounded = plyr::round_any(y, 0.5)) %>%
  group_by(y_rounded) %>%
  summarise(Epi_lead_mean = mean(Epi_lead)) %>%
  filter(abs(y_rounded) < 4) %>%
  print(n = 100)