library(tidyverse)
library(readxl)
visher <- read_xlsx("data/ppat.1005856.s007.xlsx",
                    sheet = 1) |> 
  select(`Fitness Mean`, `Fitness SD`, `Clone ID`, Dataset) |> 
  mutate(`Fitness Mean` = as.numeric(`Fitness Mean`),
         `Fitness SD` = as.numeric(`Fitness SD`)) |> 
  filter(Dataset == "Genomewide" & !is.na(`Fitness Mean`)) |> 
  mutate(class = case_when(
    `Fitness Mean` >= 0.85 & `Fitness Mean` <= 1.15 ~ "Neutral",
    `Fitness Mean` < 0.85 ~ "Deleterious",
    `Fitness Mean` > 1.15 ~ "Beneficial"
  ))

ggplot(visher, aes(x = `Fitness Mean`, fill = class)) +
  geom_histogram(color = "black", binwidth = 0.02, alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("darkred", "steelblue")) +
  labs(title = "Distribution of fitness values",
       subtitle = "Genomewide dataset") +
  annotate(geom = "text",
           x = 0.5, y = 6, label = "Mean Fitness  = 0.798",
          size = 2.8) +
  annotate(geom = "text",
           x = 0.5, y = 5, label = "Mean Deleterious Fitness = 0.608",
          size = 2.8, color = "darkred")

mean(visher$`Fitness Mean`)
visher |> 
  filter(class == "Deleterious") |> 
  pull(`Fitness Mean`) |> 
  mean()


