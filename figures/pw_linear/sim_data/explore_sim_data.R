library(tidyverse)


df <- read_csv("/Users/julianashwin/Documents/GitHub/EcoNNet.jl/figures/pw_linear/sim_data/export_data.csv")

df %>%
  mutate(y_rounded = plyr::round_any(y, 0.5)) %>%
  group_by(y_rounded) %>%
  summarise(Epi_lead_mean = mean(Epi_lead)) %>%
  filter(abs(y_rounded) < 4) %>%
  print(n = 100)


df <- read_csv("/Users/julianashwin/Documents/GitHub/EcoNNet.jl/figures/pw_linear/sim_data/many_sims_eta0p75.csv")

df_converged <- df %>%
  group_by(version) %>%
  filter(progress == max(progress))

df %>%
  filter(progress <= 100) %>%
  ggplot() + theme_bw() + 
  geom_line(aes(x = progress, y = log(Rsq_Epi), group = version, color = version))

df_converged %>%
  ggplot() + theme_bw() + 
  geom_point(aes(x = y_lag, y = pi, color = factor(version)))



mean_df <- df_converged %>%
  group_by(y_lag) %>%
  summarise(y = mean(y), pi = mean(pi), Epi_lead = mean(Epi_lead), Rsq_Epi = mean(Rsq_Epi)) %>%
  mutate(version = "Average") 

mean_df %>%
  mutate(direction = sign((y - y_lag))) %>%
  ggplot() + theme_bw() + facet_wrap(~version) + 
  geom_point(aes(x = y_lag, y = pi, color = direction))



df_converged %>%
  group_by(version) %>%
  arrange(y_lag) %>%
  mutate(direction = sign((y - y_lag))) %>%
  ggplot() + theme_bw() + facet_wrap(~version) + 
  geom_point(aes(x = y_lag, y = pi, color = direction))

mean_df %>%
  rbind(mutate(select(df_converged, -progress), version = as.character(version))) %>%
  mutate(direction = sign((y - y_lag))) %>%
  mutate(path_number = case_when(direction == -1 & y_lag > 0.5 ~ 1,
                                 direction == 1 & y_lag < -0.5 ~ 4,
                                 direction == 1 ~ 2,
                                 direction == -1 ~ 3)) %>%
  ggplot() + theme_bw() + facet_wrap(~version) + 
  #geom_point(aes(x = y_lag, y = Epi_lead, color = path_number)) 
  geom_point(aes(x = y_lag, y = path_number, color = pi))


mean_df %>%
  rbind(mutate(select(df_converged, -progress), version = as.character(version))) %>%
  mutate(direction = sign((y - y_lag))) %>%
  mutate(path_number = case_when(direction == -1 & y_lag > 0.5 ~ 1,
                                 direction == 1 & y_lag < -0.5 ~ 4,
                                 direction == 1 ~ 2,
                                 direction == -1 ~ 3)) %>%
  select(-direction) %>%
  write_csv("nnet_learning_paths_eta0p75.csv")




"
ZLB experiments
"

df <- read_csv("/Users/julianashwin/Documents/GitHub/EcoNNet.jl/figures/pw_linear/sim_data/zlb_sims_lower.csv")


df %>%
  #filter(progress == max(progress)) %>% 
  mutate(direction = sign((y - y_lag))) %>%
  ggplot() + theme_bw() + 
  facet_wrap(~pistar_low) + 
  geom_line(aes(x = y_lag, y = pi, color = factor(direction), group = interaction(progress, version)))


df_converged <- df %>%
  group_by(version) %>%
  filter(progress == max(progress))

df %>%
  filter(progress <= 100) %>%
  ggplot() + theme_bw() + 
  geom_line(aes(x = progress, y = log(Rsq_Epi), group = version, color = version))

df_converged %>%
  ggplot() + theme_bw() + 
  geom_point(aes(x = y_lag, y = pi, color = factor(version)))



mean_df <- df_converged %>%
  group_by(y_lag) %>%
  summarise(y = mean(y), pi = mean(pi), Epi_lead = mean(Epi_lead), Rsq_Epi = mean(Rsq_Epi)) %>%
  mutate(version = "Average") 

mean_df %>%
  mutate(direction = sign((y - y_lag))) %>%
  ggplot() + theme_bw() + facet_wrap(~version) + 
  geom_point(aes(x = y_lag, y = pi, color = direction))



df_converged %>%
  group_by(version) %>%
  arrange(y_lag) %>%
  mutate(direction = sign((y - y_lag))) %>%
  ggplot() + theme_bw() + facet_wrap(~version) + 
  geom_point(aes(x = y_lag, y = pi, color = direction))

mean_df %>%
  rbind(mutate(select(df_converged, -progress), version = as.character(version))) %>%
  mutate(direction = sign((y - y_lag))) %>%
  mutate(path_number = case_when(direction == -1 & y_lag > 0.5 ~ 1,
                                 direction == 1 & y_lag < -0.5 ~ 4,
                                 direction == 1 ~ 2,
                                 direction == -1 ~ 3)) %>%
  ggplot() + theme_bw() + facet_wrap(~version) + 
  #geom_point(aes(x = y_lag, y = Epi_lead, color = path_number)) 
  geom_point(aes(x = y_lag, y = path_number, color = pi))


mean_df %>%
  rbind(mutate(select(df_converged, -progress), version = as.character(version))) %>%
  mutate(direction = sign((y - y_lag))) %>%
  mutate(path_number = case_when(direction == -1 & y_lag > 0.5 ~ 1,
                                 direction == 1 & y_lag < -0.5 ~ 4,
                                 direction == 1 ~ 2,
                                 direction == -1 ~ 3)) %>%
  select(-direction) %>%
  write_csv("nnet_learning_paths_eta0p75.csv")







