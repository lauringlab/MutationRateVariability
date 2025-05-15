library(tidyverse)
library(scales)

fec_plot_poisson <- tibble(
  mean_mut_rate = c(1:10),
  fecundity = log10(1/exp(-mean_mut_rate))
)

ggplot(fec_plot_poisson, aes(x = 10^fecundity,
                     y = mean_mut_rate)) +
  geom_line() +
  theme_minimal() +
  scale_y_continuous(breaks = c(2,4,6,8,10)) +
  scale_x_log10(labels = label_log()) +
  labs(x = "fecundity", y = "mutation rate") +
  annotate(geom = "text", x = 20, y = 8, label = "Viral extinction", fontface = "italic") +
  annotate(geom = "text", x = 1500, y = 4, label = "Viral survival", fontface = "italic") +
  theme_classic()
