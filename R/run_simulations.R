# Simulate mngmt strategies -----------------------------------------------


# packages & data ---------------------------------------------------------

# packages
library(tidyverse)
library(mgcv)

source("R/simulation_functions.R")
impRes <- readRDS("data/H_saturated_markovian_nostrat_RES_MODS")
impData <- readRDS("data/H_saturated_markovian_nostrat_RES_DATA")
stratList <- readRDS("data/all_sim_strats_updated.rds")

# build strats ------------------------------------------------------------

sim_list <- list()
for(i in seq_along(stratList)){

  strat_i <- stratList[[i]]  
    
  sim_list[[i]] <- build_strategy( 
                        rotation     = strat_i$crop ,
                        soil_group   = strat_i$soil_group,
                        n_glyphosate = strat_i$a_gly,
                        spray_days   = strat_i$h_bc_app,
                        cultivation  = strat_i$cult_cat,
                        d_season     = strat_i$d_season , 
                        c_diff       = strat_i$c_diff,
                        d_diff       = strat_i$d_diff,
                        mean_mort    = strat_i$mean_mort,
                        name         = strat_i$strategy[1])
  
  initDens <- strat_i$initial_bg_density[1]
  
   switch(initDens , 
          "very high" = Nt <- matrix(c(0.1,0.1,0.2,0.3,0.3),ncol=5),
          "high"      = Nt <- matrix(c(0.2,0.2,0.2,0.2,0.2),ncol=5),
          "medium"    = Nt <- matrix(c(0.2,0.2,0.6,0,0),ncol=5),
          "low"       = Nt <- matrix(c(0.2,0.8,0,0,0),ncol=5))
  
  sim_list[[i]]$Nt <- Nt

}



# simulate strategies ----------------------------------------------------


thin   <- seq(10 , 300 , by = 3)
states <- c(1,2,3,4,5)

out <- simulate_strategy(parList = sim_list , impRes = impRes , impData = impData , thin = thin ) %>% 
        bind_rows() %>% 
        group_by(imp) %>% 
        mutate(step = row_number()) 


# summarise & plot --------------------------------------------------------



out <- out %>% select(imp , iter , mds = mds.ms , a , l , m ,h , v , strategy , year = step)
saveRDS(out , "data/simulation_results.rds")

simSum <- out %>% 
  dplyr::group_by(strategy , iter) %>% 
  dplyr::summarise(mdsm = mean(mds) , 
                   ci  = 1.96 * (sd (mds) / sqrt(n()))) %>% 
  dplyr::ungroup() %>% 
  mutate(init = str_extract(strategy , "HD|LD|MD|VD"))



simSum %>% filter(init == "HD") %>% 
    ggplot(aes(iter , mdsm,group=strategy)) +
      geom_line()+
      geom_point()+
      geom_errorbar(aes(min = mdsm - ci , max = mdsm + ci) , width = 0.2)+
      scale_y_continuous(limits  = c(0,5))+
      facet_wrap(~strategy)
