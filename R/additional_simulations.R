# Simulate mngmt strategies -----------------------------------------------

# Contained are some additional simulations requested by a reviewer.
# They include 
# 1) Additional runs at 100% high or very high densities across a filed
# 2) Sequential runs over 18 years with high resistance levels, high & very high densities, for all BAU, MIT, and BAU MIT combinations
# See supplementary material for diagramatic explanations of the simulation strategy. 


# packages & data ---------------------------------------------------------

# packages
library(tidyverse)
library(mgcv)
library(stringdist)

source("R/simulation_functions.R")
impRes  <- readRDS("data/alexa_regularised_cropping_nostrat_RES.rds_MODS")
impData <- readRDS("data/alexa_regularised_cropping_nostrat_RES.rds_DATA")
stratList_updated  <- readRDS("data/all_sim_strats_updated.rds")
stratList_HD       <- readRDS("data/all_sim_strats_HD.rds")
stratList_VD       <- readRDS("data/all_sim_strats_VHD.rds")


strats    <- list(stratList_updated , stratList_HD, stratList_VD)
saveList  <- list("data/up_simulation_results.rds", "data/up_simulation_results_HD" , "data/up_simulation_results_VD")

thin   <- seq(10 , 120, by = 3)
states <- c(1,2,3,4,5)

# get new strats + initial conditions ---------------------------------------------------------
ind              <- str_detect(names(stratList_updated) , "FF_LD_HR|FF_HD_HR")
new_inits_strats <- stratList_updated[ind] 

# build strats ------------------------------------------------------------

stratList <- new_inits_strats
save_file <- "data/up_additional_HD_VD_100_runs.rds"

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
    print(initDens)
    
    switch(as.character(initDens) ,
           "LD" = Nt <- matrix(c(0,0,0,1,0),ncol=5),
           "HD" = Nt <- matrix(c(0,0,0,0,1),ncol=5))

    sim_list[[i]]$Nt <- Nt
    
  }
  
  # Simulate
out <- simulate_strategy(parList = sim_list , impRes = impRes , impData = impData , thin = thin ) %>% 
    bind_rows() %>% 
    group_by(imp) %>% 
    mutate(step = row_number()) 
  
saveRDS(out , save_file)
  

  


# summarise & plot --------------------------------------------------------

out <- readRDS("data/up_additional_HD_VD_100_runs.rds")
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


# build  sequential rotations ---------------------------------------------------------------------


# Get above rotations from strategy names
MIT_names <- names(stratList_updated)[str_detect(names(stratList_updated) , "FF")]
BAU_names <- names(stratList_updated)[str_detect(names(stratList_updated) , "BAU")]

# MIT 18 year rotation
MIT_18       <- rep(MIT_names , each = 3) 
MIT_18_split <- split(MIT_18,ceiling(seq_along(MIT_18)/3))

# BAU 18 year rotation
BAU_18       <- rep(BAU_names , each = 3) 
BAU_18_split <- split(BAU_18,ceiling(seq_along(BAU_18)/3))

#BAU_BAU_MIT
BAU_2_MIT   <- rep(BAU_names,each = 2)
BAU_2_split <- split(BAU_2_MIT,ceiling(seq_along(BAU_2_MIT)/2))
BAU_BAU_MIT <- lapply(seq_along(BAU_2_split) , function(x) c(BAU_2_split[[x]] , MIT_names[[x]]))

#BAU_MIT_MIT
MIT_2         <- rep(MIT_names,each = 2)
BAU_MIT_split <- split(MIT_2,ceiling(seq_along(MIT_2)/2))
BAU_MIT_MIT   <- lapply(seq_along(BAU_2_split) , function(x) c(BAU_names[[x]] , BAU_MIT_split[[x]]))

# All new rotations
new_seq_list <- list(MIT_18_split , BAU_18_split , BAU_BAU_MIT , BAU_MIT_MIT)

# ---------------------------------------------------------------------------------------------


# Loop over current data and assemble new 18 year rotations
new_stratList <- list()
for(j in seq_along(new_seq_list)) {
  new_rotations <- new_seq_list[[j]] # New rotations
    rList <- list()
    for(i in seq_along(new_rotations)) {
      rList[[i]] <- stratList_updated [new_rotations[[i]]] |> 
                    bind_rows() |> ungroup() |> 
                    mutate(rotation_year = 1:n() , 
                           initial_resistance = "HR") # Set initial rotations
  
    }
    new_stratList[[j]] <- rList
}

# ---------------------------------------------------------------------------------------------

# Index of rotational coefficients that we have in the data
valid_rtns <- readRDS("data/valid_rotations.rds") |> names()

sim_list  <- list()
out_list  <- list()
dens_list <- list()


# Loop over new initial densities , strategy groups, & individual 18 year rotations
for(k in c("HD_100" , "VD_100")){
 # Switch between high and very high densities
  switch(k ,
         "HD_100" = Nt <- matrix(c(0,0,0,1,0),ncol=5),
         "VD_100" = Nt <- matrix(c(0,0,0,0,1),ncol=5))
  
for(j in seq_along(new_stratList)){
  
  # Select rotation set
  stratList <- new_stratList[[j]]
  
  for(i in seq_along(stratList)){
    
    # Build rotation
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
      name         = paste((strat_i$strategy[c(1,7,13)]) , collapse = "~"))
    
    sim_list[[i]]$Nt <- Nt
    
    # Need to match end / beginning rotational coefficients at year 8 / 12
    # List of rotations we have
    
    switch_1  <- sim_list[[i]]$rotation[c(6,7)] 
    switch_2  <- sim_list[[i]]$rotation[c(12,13)]
    
    c1 <- str_extract(switch_1[1] , "(?<=_).+")
    c2 <- str_extract(switch_1[2] , ".+(?=_)")
    c3 <- str_extract(switch_2[1] , "(?<=_).+")
    c4 <- str_extract(switch_2[2] , ".+(?=_)")

    # get appropriate swap if it exists
    if(c1!=c2)  {appR  <- valid_rtns[str_detect(valid_rtns, paste(c1,c2 , sep = "_"))]} else{appR  <- switch_1[2]}
    if(c3!=c4)  {appR2 <- valid_rtns[str_detect(valid_rtns, paste(c3,c4 , sep = "_"))]} else{appR2 <- switch_2[2]}
    if(length(appR)==0)  {print(paste("In rotation" , i , "of set" , j,  "first rotation to swap not in valid set - " , switch_1[1] ,"/", switch_1[2] , " - swapping for linseed_bar"))}
    if(length(appR2)==0) {print(paste("In rotation" , i , "of set" , j,  "second rotation to swap not in valid set - " , switch_2[1] ,"/", switch_2[2]," - swapping for linseed_bar"))}
    
    # Switch if appropriate rotations exist
     if(length(appR)  >0)   sim_list[[i]]$rotation[c(7)]   <- appR
     if(length(appR2) >0)   sim_list[[i]]$rotation[c(13)]  <- appR2
    
    # Only time this fails is where the rotation includes wheat_linseed
    if(switch_1[1]=="wheat_linseed")  sim_list[[i]]$rotation[c(7)] <- "linseed_bar"
    if(switch_2[1]=="wheat_linseed")  sim_list[[i]]$rotation[c(13)] <- "linseed_bar"
    

  }
  
  # Simulate
  print(paste0("Simulating strategies 1:", length(sim_list)  , " for set " , j , " for density " , k))
  
  out_list[[j]] <- simulate_strategy(parList = sim_list , impRes = impRes , impData = impData , thin = thin ) %>%
                    bind_rows() %>%
                    group_by(imp) %>%
                    mutate(step = row_number())
  }

  dens_list[[k]] <- out_list
  
}

# Save
saveRDS(dens_list , "data/sequential_rotation_results.rds")



# plot ----------------------------------------------------------------------------------------

simSum <- dens_list$VD_100[[4]] %>% 
  mutate(iter=as.numeric(iter)) |> 
  dplyr::group_by(strategy , iter) %>% 
  dplyr::summarise(mdsm = mean(mds.ms) , 
                   ci  = 1.96 * (sd (mds.ms) / sqrt(n()))) %>% 
  dplyr::ungroup() %>% 
  mutate(init = str_extract(strategy , "HD|LD|MD|VD"),
         init_res = str_extract(strategy , "HR|LR"))


simSum %>% 
  ggplot(aes(iter , mdsm,group=strategy)) +
  geom_line(aes(colour=init_res))+
  geom_point()+
  geom_errorbar(aes(min = mdsm - ci , max = mdsm + ci) , width = 0.2)+
  scale_y_continuous(limits  = c(0,5))+
  facet_grid(init_res~init)

