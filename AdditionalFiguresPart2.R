library(tidyverse)
library(extraDistr)
library(data.table)
library(bench)
source("analyses/ConditionalProbability.R")
source("analyses/Figure1.R")

## Figure 4A: Extinction threshold by index of dispersion ##

#function to get extinction fecundity given shape, rate
fecundity <- function(shape, rate){
  return((1+(1/rate))^shape)
}

#mean mutation rates and indices of dispersion to plot
mean_mut_rate <- 1:10
overdispersion_index <- c(1.16, 1.24)

#calculate thresholds
fec_plot <- expand_grid(mean_mut_rate, overdispersion_index) %>%
  mutate(shape = shape_gpois(mean_mut_rate, overdispersion_index),
         rate = shape/mean_mut_rate,
         fecundity = fecundity(shape, rate))

#plot thresholds
ggplot(fec_plot, aes(x = fecundity,
                     y = mean_mut_rate,
                     group = overdispersion_index)) +
  geom_line() +
  geom_line(data = fec_plot_poisson, aes(x = fecundity,
                                         y = mean_mut_rate,
                                         group = 1),
            color = "black") +
  scale_y_continuous(breaks = c(2,4,6,8,10)) +
  scale_x_log10(labels = label_log()) +
  labs(x = "fecundity", y = "mean mutation rate",
       color = "index of dispersion") +
  annotate(geom = "text", x = 20, y = 9, label = "Viral extinction", fontface = "italic") +
  annotate(geom = "text", x = 1500, y = 4, label = "Viral survival", fontface = "italic") +
  theme_classic()

#initial population with 100 individuals, no mutations, fitness = 1
population0 <- tibble(
  mu = rep(0, 100),
  fitness = rep(1, 100)
)

#Function to generate progeny of a population
#no heredity of mutation rate, but there is accumulation of mutations
progeny_gamma <- function(pop, fecundity, shape, rate, s){
  #vector of progeny
  new_pop <- vector(mode = "list", length = nrow(pop)*fecundity)
  
  #for every individual in the current generation, make a number of offspring
  #equal to the fecundity times the fitness of that generation.
  #Each offpsring has mutations according to rgpois().
  #Then calculate fitness for each offspring based on mutation number.
  for(i in 1:nrow(pop)){
    progeny_mu = rgpois(floor(fecundity*pop$fitness[i]),
                        shape = shape, rate = rate) + pop$mu[i]
    progeny_fitness = (1-s)^progeny_mu
    non_zero_fitness <- progeny_fitness > 0
    new_pop[[i]] <- list(mu = progeny_mu[non_zero_fitness],
                         fitness = progeny_fitness[non_zero_fitness])
  }
  return(rbindlist(new_pop))
}

#Simulate populations according to gamma-Poisson, agent-based
#assume viruses only produce offspring once and die after producing offspring
mu_sim_gamma <- function(sim, initial_pop, n_generations, fecundity, shape, rate, s){
  #Monitor progress
  if(sim %% 10 == 0){
    print(paste("Shape:", shape, "Sim:", sim))
  }
  
  #initialize population
  generations = 1:n_generations
  pop_results <- list()
  pop_results[[1]] <- initial_pop
  
  #For each generation, calculate the progeny and add to results
  #If no progeny, break
  for(i in generations){
    pop_results[[i+1]] = progeny_gamma(pop_results[[i]], fecundity,shape, rate,s)
    #print(nrow(pop_results[[i+1]]))
    if(nrow(pop_results[[i+1]]*fecundity*mean(pop_results[[i+1]]$fitness)) < 1)
      break
  }
  
  #Summarize results into a list of tibbles
  pop_results_summary <- list()
  for(i in 1:length(pop_results)){
    pop_results_summary[[i]] = 
      tibble(pop_size = nrow(pop_results[[i]]),
             mean_mu = mean(pop_results[[i]]$mu),
             mean_fitness = mean(pop_results[[i]]$fitness))
  }
  
  #Collapse list into a single tibble 
  pop_results_summary <- rbindlist(pop_results_summary, fill = TRUE) %>%
    rownames_to_column() %>%
    mutate(generation = as.numeric(rowname)-1) %>%
    select(-rowname)
  
  return(pop_results_summary)
}

#get shape parameters for mean of 3, oi of 1.0 to 1.6
shapes <- shape_gpois(3, c(1, 1.16,1.24))

#set shape for mean=3, oi=0 to arbitrarily large value bc can't actually
#calculate shape when oi=0
shapes[1] <- 300

#function to run simulations for a given shape/rate and record time to extinction
run_sim_extinction <- function(n_sims, shape, rate){
  sim_death <- c()
  for(i in 1:n_sims){
    sim <- mu_sim_gamma(i, population0, 150, 10, shape, rate, 0.4)
    sim_death[i] <- sim %>% 
      summarize(max_gen = max(generation)) %>%
      pull(max_gen)
  }
  return(sim_death)
}

#Run 100 simulations per parameter set
extinction_gens <- list()
for(i in 1:length(shapes)){
  print(paste("Running sim with shape = ", shapes[i]))
  extinction_gens[[i]] <- run_sim_extinction(100, shapes[i], shapes[i]/3)
  print(mean(extinction_gens[[i]]))
  
}

#Combine simulation results into a tibble
extinction_gens_comb <- tibble(
  o_index = sort(rep(c(1.0,1.16,1.24), 100)),
  extinction_gen = unlist(extinction_gens)
)

#Save results
#saveRDS(extinction_gens_comb, "data/extinction_time_sim_extra.rds")
extinction_summary <- extinction_gens_comb %>% 
  group_by(o_index) %>% 
  summarize(mean_extinction = mean(extinction_gen))

ggplot(extinction_gens_comb, aes(x = o_index, y = extinction_gen, group = o_index)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 50)) +
  scale_x_continuous(breaks = c(1.0, 1.16, 1.24)) +
  annotate("text", x = c(1, 1.16, 1.24), y = 35,
           label = extinction_summary$mean_extinction) +
  theme_classic() +
  labs(x = "index of dispersion",
       y = "time to extinction") +
  theme(legend.position = "none")

