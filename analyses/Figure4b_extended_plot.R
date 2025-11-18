library(tidyverse)

d <- readRDS("data/sim_results_20250812/p_100_0.2_2_4")

ggplot(d, aes(x = o_index, y = extinction_gen, group = o_index, color = o_index)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 150)) +
  scale_x_continuous(breaks = c(1.0, 1.2, 1.4, 1.6)) +
  scale_color_gradient(low = "#fcbba1", high = "#cb181d") +
  theme_classic() +
  labs(x = "index of dispersion",
       y = "time to extinction") +
  theme(legend.position = "none")

dir_list <- list.files(path = "data/sim_results_20250821/",
full.names = TRUE)
names(dir_list) <- list.files(path = "data/sim_results_20250821/")
files_df <- map_dfr(dir_list, read_rds, .id = "source")

mmr2 <- files_df |> 
  filter(str_detect(source, "2_4")) |> 
  mutate(source = str_remove(source, "_2_4"),
         source = str_remove(source, "p_"))

ggplot(mmr2, aes(x = o_index, y = extinction_gen, group = o_index, color = o_index)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 90)) +
  scale_x_continuous(breaks = c(1.0, 1.2, 1.4, 1.6)) +
  scale_color_gradient(low = "#fcbba1", high = "#cb181d") +
  theme_classic() +
  labs(x = "index of dispersion",
       y = "time to extinction",
      title = "mean mutation rate = 2, fecundity = 4") +
  theme(legend.position = "none") +
  facet_wrap(~source)

mmr3 <- files_df |> 
  filter(str_detect(source, "3.17_10")) |> 
  mutate(source = str_remove(source, "_3.17_10"),
         source = str_remove(source, "p_"))

ggplot(mmr3, aes(x = o_index, y = extinction_gen, group = o_index, color = o_index)) +
  geom_jitter(alpha = 0.5) +
  scale_y_continuous(limits = c(0, 90)) +
  scale_x_continuous(breaks = c(1.0, 1.2, 1.4, 1.6)) +
  scale_color_gradient(low = "#fcbba1", high = "#cb181d") +
  theme_classic() +
  labs(x = "index of dispersion",
       y = "time to extinction",
      title = "mean mutation rate = 3.17, fecundity = 10") +
  theme(legend.position = "none") +
  facet_wrap(~source)
