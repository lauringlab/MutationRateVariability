library(tidyverse)
library(gridExtra)
library(scales)


## Figure 3B: Mutation accumulation experiment ##

#Load initial processing of sequencing data
source("analyses/SequencingAnalysis.R")
source("analyses/ConditionalProbability.R")

#Separate conditions into different tables
#control
ctl <- mut_per_sample %>% 
  filter(str_detect(condition, "Ctl")) %>% 
  pull(n_mutations)

#calculate bootstrapped 95% confidence intervals
lethal_fraction = 0.4

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

ctl_95_uc <- quantile(ctl_bs_uc$oi,c(0.05, 0.95), na.rm = TRUE)

ci_uc <- tibble(
  condition = factor(c("Ctl"), levels = c("Ctl")),
  q5 = c(ctl_95_uc[1]),
  q95 = c(ctl_95_uc[2])
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


