library(tidyverse)

#Find empty samples
coverage <- read_csv("data/coverage.csv") %>% 
  rename(sample = ID)

mean_coverage <- coverage %>% 
  group_by(sample) %>% 
  summarize(mean_cov = mean(cov))

all_samples <- tibble(
  condition = c(rep("Ctl_rep1", 82),
                rep("Ctl_rep2", 82),
                rep("5um_rep1", 62),
                rep("5um_rep2", 62),
                rep("10um_rep1", 48),
                rep("10um_rep2", 48)),
  sample = c(1:82, 1:82,
             1:62, 1:62,
             1:48, 1:48)
) %>% 
  mutate(sample = paste(condition, sample, sep = "_"))

empty_samples <- anti_join(all_samples, mean_coverage, by = "sample")

#Read mutation frequency file
mut_freqs <- read_delim("data/all_variants_filtered", delim = " ") %>% 
  select(REGION, POS, ALT_FREQ, mutation, sample) %>% 
  unique() %>% #because NS1 and NS2 coding regions overlap
  mutate(sample = str_remove(sample, ".variants.tsv"),
         condition = str_extract(sample, "[0-9]*[A-z]*_rep[0-9]"))

#Look for reference sequence errors
mut_freq_1 <- mut_freqs %>% 
  filter(ALT_FREQ == 1) %>%
  group_by(mutation) %>% 
  summarize(n_samples = n())

#Only need to mask PB2_a1700G, the rest are in 3 or fewer samples
mut_freqs <- mut_freqs %>% 
  filter(mutation != "PB2_a1700G")

#check for multiple infections
#(should have multiple mutations > ~15% AND no mutations near fixed)
min_freq = 0.15
max_freq = 0.95

multi_infected <- mut_freqs %>% 
  mutate(multi = ifelse(
    ALT_FREQ > min_freq & ALT_FREQ < max_freq, TRUE, FALSE)) %>% 
  group_by(sample) %>% 
  summarize(multi_s = sum(multi)) %>% 
  filter(multi_s > 0)

#Select all mutations above 95%
mut_freq_fixed <- mut_freqs %>% 
  filter(ALT_FREQ >= max_freq)

#Count fixed mutations per sample
mut_per_sample <- mut_freq_fixed %>% 
  group_by(sample) %>% 
  summarize(n_mutations = n(), condition = unique(condition)) %>% 
  #Make sure zero mutation samples get counted
  full_join(all_samples, by = c("sample", "condition")) %>%
  replace_na(list(n_mutations = 0)) %>% 
  #Exclude empty and multi-infected samples
  filter(!(sample %in% empty_samples$sample) &
           !(sample %in% multi_infected$sample))

overdispersion_full <- mut_per_sample %>% 
  mutate(condition = str_split_i(condition, "_", 1),
         condition = factor(condition, levels = c("Ctl", "5um", "10um"))) %>% 
  ungroup() %>% 
  group_by(condition) %>% 
  summarize(mean_mut = mean(n_mutations),
            var_mut = var(n_mutations),
            oi = var_mut/mean_mut,
            n_obs = n())
