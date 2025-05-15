library(tidyverse)
library(gridExtra)

## Figure 3A: IAV growth curve ##
gc <- read_csv("data/audrey_growthcurve_rep.csv") %>%
  mutate(Replicate = as.character(Replicate))

gc_summary <- gc %>% 
  group_by(`Timepoint (hpi)`) %>% 
  summarize(mean_titer = mean(Titer),
            sd_titer = sd(Titer))

ggplot(gc_summary, aes(x = `Timepoint (hpi)`, y = mean_titer)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_titer - sd_titer,
                    ymax = mean_titer + sd_titer),
                width = 0.2) +
  scale_x_continuous(breaks = c(3:12)) +
  scale_y_log10(limits = c(1e3, 3e5),
                breaks = c(1e3, 1e4, 1e5),
                labels = label_log()) +
  labs(x = "timepoint (hpi)",
       y = "titer (TCID50)") +
  theme_classic() +
  annotation_logticks(sides = "l") +
  geom_vline(xintercept = 8, linewidth = 10, color = "steelblue",
             alpha = 0.2)


## Figure 3B: Mutation accumulation experiment ##

#Load initial processing of sequencing data
source("analyses/SequencingAnalysis.R")

#Separate conditions into different tables
#control
ctl <- mut_per_sample %>% 
  filter(str_detect(condition, "Ctl")) %>% 
  pull(n_mutations)

#5uM ribavirin
r5 <- mut_per_sample %>% 
  filter(str_detect(condition, "5um")) %>% 
  pull(n_mutations)

#10uM ribavirin
r10 <- mut_per_sample %>% 
  filter(str_detect(condition, "10um")) %>% 
  pull(n_mutations)

#Function to calculate bootstrapped confidence intervals
bootstrap <- function(n, mut_counts, lethal_frac){
  #pre-allocate
  bs_mean <- vector(mode = "numeric", length = n)
  bs_var <- vector(mode = "numeric", length = n)
  for(i in 1:n){
    s <- sample(mut_counts, n, replace = TRUE)
    bs_mean[i] <- mean(s)
    bs_var[i] <- var(s)
  }
  return(tibble(m = bs_mean, v = bs_var, oi = v/m))
}

#calculate bootstrapped 95% confidence intervals
lethal_fraction = 0.3

ctl_bs <- bootstrap(1000, ctl, lethal_fraction)
r5_bs <- bootstrap(1000, r5, lethal_fraction)
r10_bs <- bootstrap(1000, r10, lethal_fraction)

ctl_95 <- quantile(ctl_bs$oi,c(0.05, 0.95))
r5_95 <- quantile(r5_bs$oi, c(0.05, 0.95))
r10_95 <- quantile(r10_bs$oi, c(0.05, 0.95))
ci <- tibble(
  condition = factor(c("Ctl", "5um", "10um"), levels = c("Ctl", "5um", "10um")),
  q5 = c(ctl_95[1], r5_95[1], r10_95[1]),
  q95 = c(ctl_95[2], r5_95[2], r10_95[2])
)

#plot data and confidence intervals
overdispersion_ci <- overdispersion_full %>% 
  right_join(ci, by = "condition")

ggplot(overdispersion_ci, aes(x = condition, y = oi)) +
  geom_point(size = 3) +
  geom_text(aes(label = round(oi,2)), hjust = -0.5, size = 6) +
  geom_linerange(aes(ymin = q5, ymax = q95)) +
  geom_hline(yintercept = 1, alpha = 0.5, linetype = 5) +
  theme_classic() +
  scale_y_continuous(limits = c(0.8, 1.5)) +
  labs(x = "ribavirin",
       y = "index of dispersion",
      title = "Sampled (viable)")

#Infer values for full population

#function for unconditional bootstrap
bootstrap_uc <- function(n, mut_counts, lf){
  #pre-allocate
  bs_mean <- vector(mode = "numeric", length = n)
  bs_var <- vector(mode = "numeric", length = n)
  for(i in 1:n){
    s <- sample(mut_counts, n, replace = TRUE)
    bs_mean[i] <- exp_x_uc(mean(s), var(s), lf)
    bs_var[i] <- var_x_uc(mean(s), var(s), lf)
  }
  return(tibble(m = bs_mean, v = bs_var, oi = v/m))
}

#calculate inferred values
ctl_bs_uc <- bootstrap_uc(1000, ctl, lethal_fraction)
r5_bs_uc <- bootstrap_uc(1000, r5, lethal_fraction)
r10_bs_uc <- bootstrap_uc(1000, r10, lethal_fraction)

ctl_95_uc <- quantile(ctl_bs_uc$oi,c(0.05, 0.95), na.rm = TRUE)
r5_95_uc <- quantile(r5_bs_uc$oi, c(0.05, 0.95), na.rm = TRUE)
r10_95_uc <- quantile(r10_bs_uc$oi, c(0.05, 0.95), na.rm = TRUE)
ci_uc <- tibble(
  condition = factor(c("Ctl", "5um", "10um"), levels = c("Ctl", "5um", "10um")),
  q5 = c(ctl_95_uc[1], r5_95_uc[1], r10_95_uc[1]),
  q95 = c(ctl_95_uc[2], r5_95_uc[2], r10_95_uc[2])
)

#plot inferred data and confidence intervals
overdispersion_ci_uc <- overdispersion_full %>% 
  right_join(ci_uc, by = "condition") %>% 
  mutate(oi_uc = var_x_ucV(mean_mut, var_mut, lethal_fraction)/
           exp_x_ucV(mean_mut, var_mut, lethal_fraction),
         q5 = ifelse(is.na(q5), 0.8, q5))

ggplot(overdispersion_ci_uc, aes(x = condition, y = oi_uc)) +
  geom_point(size = 3) +
  geom_text(aes(label = round(oi_uc,2)), hjust = -0.5, size = 6) +
  geom_linerange(aes(ymin = q5, ymax = q95)) +
  geom_hline(yintercept = 1, alpha = 0.5, linetype = 5) +
  theme_classic() +
  scale_y_continuous(limits = c(0.8, 1.5)) +
  labs(x = "ribavirin",
       y = "index of dispersion",
       title = "Inferred")


## Figure 3C: Comparing distributions ##

#Fit gamma-Poisson for control samples
shape = shape_gpois(mean(ctl), var(ctl)/mean(ctl))
rate = shape/mean(ctl)

gpois_compare <- as_tibble(table(ctl, dnn = "n_mutations")) %>% 
  mutate(condition = "IAV", n_mutations = as.numeric(n_mutations)) %>%
  bind_rows(tibble(condition = rep("gamma-Poisson", 5),
                   n_mutations = 0:4,
                   n = round(dgpois(0:4,shape = shape, rate = rate) * length(ctl)))) %>% 
  mutate(condition = factor(condition, levels = c("IAV", "gamma-Poisson")))

gpois_plot <- ggplot(gpois_compare, aes(x = n_mutations, y = n, fill = condition)) +
  geom_col(position = "identity", alpha = 0.5) +
  theme_classic() +
  theme(legend.position = "top") +
  labs(x = "number of mutations", y = "count", fill = "") +
  scale_fill_brewer(palette = "Set1")

#Fit poisson for control samples
pois_compare <- as_tibble(table(ctl, dnn = "n_mutations")) %>% 
  mutate(condition = "IAV", n_mutations = as.numeric(n_mutations)) %>%
  bind_rows(tibble(condition = rep("Poisson", 5),
                   n_mutations = 0:4,
                   n = round(dpois(0:4, mean(ctl)) * length(ctl)))) %>% 
  mutate(condition = factor(condition, levels = c("IAV", "Poisson")))

pois_plot <- ggplot(pois_compare, aes(x = n_mutations, y = n, fill = condition)) +
  geom_col(position = "identity", alpha = 0.5) +
  theme_classic() +
  theme(legend.position = "top") +
  labs(x = "number of mutations", y = "count", fill = "") +
  scale_fill_brewer(palette = "Set1")

#Calculate difference in expected vs observed
gpois_dif <- gpois_compare %>% 
  pivot_wider(names_from = condition, values_from = n) %>% 
  mutate(dif = IAV-Gammapoisson)

pois_dif <- pois_compare %>% 
  pivot_wider(names_from = condition, values_from = n) %>% 
  mutate(dif = IAV-`Poisson`)

#plot differences
gpois_dif_plot <- ggplot(gpois_dif, aes(x = n_mutations, y = dif)) +
  geom_col() +
  theme_classic() +
  geom_hline(yintercept = 0) +
  scale_y_continuous(limits = c(-5, 5), breaks = c(-5, 0, 5)) +
  labs(y = "Difference", x = NULL) +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

pois_dif_plot <- ggplot(pois_dif, aes(x = n_mutations, y = dif)) +
  geom_col() +
  theme_classic() +
  geom_hline(yintercept = 0) +
  scale_y_continuous(limits = c(-5, 5), breaks = c(-5, 0, 5)) +
  labs(y = "difference", x = NULL) +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

#combine plots
grid.arrange(pois_plot, gpois_plot,
  pois_dif_plot, gpois_dif_plot,
  ncol = 2, heights = c(4.5, 1))

## Calculate P(alive) ##
palive_ctl <- 1/calc_sum(mean(ctl), var(ctl), lethal_fraction)
