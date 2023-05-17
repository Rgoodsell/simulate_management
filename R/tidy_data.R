
# packages ----------------------------------------------------------------

library(tidyverse)
library(janitor)
library(lubridate)

source("R/simulation_functions.R")

# data --------------------------------------------------------------------

# General tidying
init_data <- readRDS("data/rotation_data_for_Rob_2022-09-28.rds") %>% clean_names() %>%
             mutate(strategy = strategy_id_for_rob , 
                    soil_group = soil_group_superseded) %>% 
            mutate(crop = case_when(rotation_year == 0 & crop == "winter OSR" ~ "winter wheat" , TRUE ~ crop)) %>% 
  filter(rotation_year>0)
  


# tidy --------------------------------------------------------------------

tidy_data <-   init_data %>% 
              mutate(d_season = case_when( str_detect(crop , "winter") ~ "winter" , 
                                           str_detect(crop , "spring") ~ "spring" , 
                TRUE ~ NA_character_)) %>%
                mutate(crop = recode(crop , 
                                   "18 month fallow" = "fallow" , 
                                   "spring barley" = "bar",
                                   "winter barley" = "bar",
                                    "spring beans" = "beans" , 
                                   "spring linseed" = "linseed" , 
                                   "spring oats" = "oats" , 
                                   "sugar beet" = "beet" , 
                                   "winter OSR" = "osr" , 
                                   "winter wheat" = "wheat",
                                   "maize" = "osr") , 
                       cult_cat = recode(cult_cat , "none" = "surface")) %>% 
  mutate(d_season = factor(case_when(d_date < 183 ~ "spring" , d_date > 183 ~ "winter")),
         d_date = as.numeric(d_date) , 
         c_date = as.numeric(c_date) ,
         h_date = as.numeric(h_date),
         a_gly = as.numeric(a_gly),
        spray_days = as.numeric(spray_days)) %>% 
  mutate(c_date = case_when(cult_cat == "none" & d_season == "spring" ~ 251,
                            cult_cat == "none" & d_season == "winter" ~ 240,
                            TRUE ~ c_date)) %>% 
  mutate(d_diff = case_when(d_season == "spring" ~ d_date - 82, 
                            d_season == "winter" ~ d_date - 273),
         c_diff = case_when(d_season == "spring" ~ c_date - 251, 
                            d_season == "winter" ~ c_date - 240),
         mean_mort = case_when(initial_resistance  == "HR" ~ 13 , 
                               initial_resistance == "LR"   ~ 63 , 
                               is.na(initial_resistance) ~ 90)) %>% 
  ungroup() %>% 
  mutate(strategy_split = interaction(strategy , initial_bg_density)) %>% droplevels() %>% 
  group_by(strategy) %>%  mutate(fcs = paste(crop , collapse = "_")) %>% distinct() %>% ungroup()



# replace missing coefficients
tidy_data <- tidy_data %>%
          group_by(strategy) %>%
          mutate(rotation = paste0(lag(crop , default = "wheat")  , "_" , crop)) %>%
          mutate(rotation = case_when(rotation == "fallow_osr"    ~ "fallow_wheat",
                                      rotation == "osr_bar"       ~ "osr_wheat"   ,
                                      rotation == "bar_osr"       ~ "bar_beans"   ,
                                      rotation == "linseed_beans" ~ "linseed_bar" ,
                                      rotation == "linseed_peas"  ~ "linseed_bar" ,
                                      rotation == "beet_oats"     ~ "beet_wheat"  ,
                                      rotation == "beat_oats"     ~ "beet_wheat"  ,
                                      rotation == "peas_bar"      ~ "peas_wheat"  ,
                                      rotation == "bar_peas"      ~ "bar_beans"   ,
                                      TRUE ~ rotation))


# Split 
sList    <- split(tidy_data , tidy_data$strategy_split)
sList_HD <- tidy_data %>% mutate(initial_bg_density = "HD") %>% split(. , .$strategy_split)
sList_VHD <- tidy_data %>% mutate(initial_bg_density = "VD") %>% split(. , .$strategy_split)


saveRDS(sList , "data/all_sim_strats_updated.rds")
saveRDS(sList_HD , "data/all_sim_strats_HD.rds")
saveRDS(sList_VHD , "data/all_sim_strats_VHD.rds")


