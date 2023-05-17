# Simulate mngmt strategies -----------------------------------------------


# packages & data ---------------------------------------------------------

# packages
library(tidyverse)
library(mgcv)

source("R/simulation_functions.R")
impRes  <- readRDS("data/alexa_regularised_cropping_nostrat_RES.rds_MODS")
impData <- readRDS("data/alexa_regularised_cropping_nostrat_RES.rds_DATA")
stratList_updated  <- readRDS("data/all_sim_strats_updated.rds")
stratList_HD       <- readRDS("data/all_sim_strats_HD.rds")
stratList_VD       <- readRDS("data/all_sim_strats_VHD.rds")


strats <- list(stratList_updated , stratList_HD, stratList_VD)
saveList  <- list("data/up_simulation_results.rds", "data/up_simulation_results_HD" , "data/up_simulation_results_VD")

thin   <- seq(10 , 120, by = 3)
states <- c(1,2,3,4,5)

# build strats ------------------------------------------------------------

for(j in seq_along(strats)){
  
stratList <- strats[[j]]  
saveFile  <- saveList[[j]] 
  
  sim_list <- list()
  # add cropping and rotation vars
  
  for(i in seq_along(stratList)){
  
    strat_i <- stratList[[i]]
  
    sim_list[[i]] <- build_strategy(
                          cropping     = strat_i$crop,
                          rotation     = strat_i$rotation ,
                          soil_group   = strat_i$soil_group,
                          n_glyphosate = strat_i$a_gly,
                          spray_days   = strat_i$spray_days,
                          cultivation  = strat_i$cult_cat,
                          d_season     = strat_i$d_season ,
                          c_diff       = strat_i$c_diff,
                          d_diff       = strat_i$d_diff,
                          mean_mort    = strat_i$mean_mort,
                          name         = strat_i$strategy[1])
  
    initDens <- strat_i$initial_bg_density[1]
  
     switch(initDens ,
            "VD" = Nt <- matrix(c(0.1,0.1,0.2,0.3,0.3),ncol=5),
            "HD" = Nt <- matrix(c(0.2,0.2,0.2,0.2,0.2),ncol=5),
            "MD" = Nt <- matrix(c(0.2,0.2,0.6,0,0),ncol=5),
            "LD" = Nt <- matrix(c(0.2,0.8,0,0,0),ncol=5))
  
    sim_list[[i]]$Nt <- Nt
  
  }
  
# Simulate
out <- simulate_strategy(parList = sim_list , impRes = impRes , impData = impData , thin = thin ) %>% 
        bind_rows() %>% 
        group_by(imp) %>% 
        mutate(step = row_number()) 

saveRDS(out , saveFile)

}

library(beepr)
beepr::beep()

# summarise & plot --------------------------------------------------------

out <- readRDS("data/simulation_results.rds")
out <- out %>% select(imp , iter , mds = mds.ms , a , l , m ,h , v , strategy , year = step) 

simSum <- out %>% 
  dplyr::group_by(strategy , iter) %>% 
  dplyr::summarise(mdsm = mean(mds) , 
                   ci  = 1.96 * (sd (mds) / sqrt(n()))) %>% 
  dplyr::ungroup() %>% 
  mutate(init = str_extract(strategy , "HD|LD|MD|VD"),
         init_res = str_extract(strategy , "HR|LR"))


simSum %>% filter(str_detect(strategy , "FF")) %>% 
    ggplot(aes(iter , mdsm,group=strategy)) +
      geom_line(aes(colour=init_res))+
      geom_point()+
      geom_errorbar(aes(min = mdsm - ci , max = mdsm + ci) , width = 0.2)+
      scale_y_continuous(limits  = c(0,5))+
      facet_grid(init_res~init)


