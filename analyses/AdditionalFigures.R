library(tidyverse)
source("analyses/SequencingAnalysis.R")

# Overdispersion in polymerase mutants vs everyone else
  
mut_polymerase <- mut_freq_fixed |> 
  mutate(polymerase = ifelse(REGION %in% c("PB2", "PB1", "PA"), TRUE, FALSE)) |> 
  group_by(sample) |> 
  summarize(n_mutations = n(), condition = unique(condition), polymerase = ifelse(max(polymerase >0), "At least 1 polymerase mutation", "Other mutations"))

#zero mutants
zero_mutants <- all_samples |> 
  anti_join(mut_polymerase)

#Make sure zero mutation samples get counted
mut_all <- mut_polymerase |> 
  bind_rows(zero_mutants) %>%
  replace_na(list(n_mutations = 0, polymerase = "Other mutations")) %>% 
  bind_rows(zero_mutants) %>%
  replace_na(list(n_mutations = 0, polymerase = "At least 1 polymerase mutation")) %>%
  #Exclude empty and multi-infected samples
  filter(!(sample %in% empty_samples$sample) &
         !(sample %in% multi_infected$sample)) |> 
    mutate(class = ifelse(n_mutations == 0, TRUE, FALSE))

mut_summary <- mut_all |> 
  filter(!class) |> 
  select(condition, polymerase, n_mutations) |> 
  mutate(condition = str_sub(condition, start = 1, end = 3)) |> 
  group_by(condition, polymerase) |> 
  summarize(mean_mut = mean(n_mutations))


ggplot(mut_all |> filter(str_detect(condition, "Ctl")), aes(x = n_mutations,
       fill = class)) +
  geom_bar(position = "identity", alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("steelblue", "grey")) +
  labs(x = "number of mutations per sample", title = "Control") +
  facet_wrap(~polymerase) +
  theme(legend.position = "none")


ggplot(mut_all |> filter(str_detect(condition, "5um")), aes(x = n_mutations,
  fill = class)) +
geom_bar(position = "identity", alpha = 0.5) +
theme_classic() +
scale_fill_manual(values = c("steelblue", "grey")) +
labs(x = "number of mutations per sample", title = "5uM ribavirin") +
facet_wrap(~polymerase) +
theme(legend.position = "none")

ggplot(mut_all |> filter(str_detect(condition, "10um")), aes(x = n_mutations,
  fill = class)) +
geom_bar(position = "identity", alpha = 0.5) +
theme_classic() +
scale_x_continuous(breaks = c(0:4), limits = c(-0.5,4)) +
scale_fill_manual(values = c("steelblue", "grey")) +
labs(x = "number of mutations per sample", title = "10uM ribavirin") +
facet_wrap(~polymerase) +
theme(legend.position = "none")

#Probability of viable given X = m

source("analyses/ConditionalProbability.R")
p_viable <- tibble(
  n_mutations = c(0:10),
  p = p_alive_given_m(n_mutations, 0.3)
)

ggplot(p_viable, aes(x = n_mutations, y = p)) +
  geom_col() +
  theme_classic() +
  scale_x_continuous(breaks = c(0:10))

p_m_given_alive_poisson <- function(m, m_mean, p_lethal){
  dpois(m, m_mean) * p_alive_given_m(2, p_lethal)
}

p_m_viable <- tibble(
  n_mutations = rep(c(0:10), 3),
  mu = c(rep(1, 11), rep(2, 11), rep(3, 11)),
  p = p_m_given_alive_poisson(n_mutations, mu, 0.3)
)

ggplot(p_m_viable, aes(x = n_mutations, y = p)) +
  geom_col() +
  theme_classic() +
  scale_x_continuous(breaks = c(0:10)) +
  facet_wrap(~mu)
