library(tidyverse)
library(extraDistr)
library(data.table)
library(bench)
source("analyses/ConditionalProbability.R")
source("analyses/Figure1.R")
## Figure 4B: time to extinction simulations ##

#Parameters
  #Starting pop size: 10, 100, 1000
    starting_pop_size <- c(10, 100, 1000)
  #Fitness cost: 0.2, 0.4, 1
    fitness_cost <- c(0.2, 0.4, 1)
  #Mean mutation rate: 2, 3, 6
    mean_mutation_rate <- c(2, 3, 6)
  #Fecundity: 4, 10, 100, 1000
    fecundity <- c(4, 10, 100)
  #Index of dispersion: 1.0 - 1.6
    #specified in function

parameter_combos <- expand.grid(starting_pop_size, fitness_cost, mean_mutation_rate) |> 
  mutate(Var4 = case_when(
    Var3 == 2 ~ 4,
    Var3 == 3 ~ 10,
    Var3 == 6 ~ 100
  ))

run_parameter_combo <- function(spz, fc, mmr, f){
  population0 <- tibble(
    mu = rep(0, spz),
    fitness = rep(1, spz)
  )

  #get shape parameters for mean of 3, oi of 1.0 to 1.6
  shapes <- shape_gpois(mmr, seq(1.0, 1.6, by = 0.2))

  #set shape for mean=3, oi=0 to arbitrarily large value bc can't actually
  #calculate shape when oi=0
  shapes[1] <- mmr*1000

  extinction_gens <- list()
  for(i in 1:length(shapes)){
    print(paste("Running sim with shape = ", shapes[i]))
    extinction_gens[[i]] <- run_sim_extinction(100, population0, shapes[i], shapes[i]/3, f, fc, 200)
    print(mean(extinction_gens[[i]]))

  }

  extinction_gens_comb <- tibble(
    o_index = sort(rep(seq(1.0, 1.6, 0.2), 100)),
    extinction_gen = unlist(extinction_gens)
  )

  return(extinction_gens_comb)
}

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
mu_sim_gamma <- function(sim, initial_pop, max_generations, fecundity, shape, rate, s){
  #Monitor progress
  if(sim %% 10 == 0){
    print(paste("Shape:", shape, "Sim:", sim))
  }

  #initialize population
  generations = 1:max_generations
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

#function to run simulations for a given shape/rate and record time to extinction
run_sim_extinction <- function(n_sims, population0, shape, rate, fecundity, s, max_generations){
  sim_death <- c()
  for(i in 1:n_sims){
    sim <- mu_sim_gamma(i, population0, max_generations, fecundity, shape, rate, s)
    sim_death[i] <- sim %>% 
      summarize(max_gen = max(generation)) %>%
      pull(max_gen)
  }
  return(sim_death)
}

for(i in 1:nrow(parameter_combos)){
  p <- parameter_combos[i,]
  result <- run_parameter_combo(p$Var1, p$Var2, p$Var3, p$Var4)
  filename <- paste("data/extinction_time_simulations/p", p$Var1, p$Var2, p$Var3, p$Var4, sep = "_")
  saveRDS(result, filename)
}