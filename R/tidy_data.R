
# packages ----------------------------------------------------------------

library(tidyverse)
library(janitor)
library(lubridate)

source("R/simulation_functions.R")

# data --------------------------------------------------------------------

my_data <- read_csv("data/management_data.csv") %>% clean_names() 

old_soil <- read_csv("data/mgmtdata_strats.csv") %>% select(soil_group , strategy = scenario) %>% distinct()
new_soil <- my_data %>% select(strategy ,s4:s8) %>% distinct()

soil_lookup <- full_join(old_soil , new_soil) %>%
                    select(-strategy) %>% distinct() %>% arrange(soil_group) %>% drop_na() %>% 
                    mutate(iid = interaction(s4,s5,s7, s8)) %>% 
                    select(-s4:-s8)
                    

# tidy --------------------------------------------------------------------


tidy_data <-   my_data %>% 
                mutate(iid = interaction(s4,s5,s7, s8)) %>% left_join(soil_lookup) %>% 
                select(-s4:-s8 , iid) %>% 
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
         mean_mort = case_when(initial_resistance  == "high" ~ 20 , 
                               initial_resistance == "low"   ~ 90 , 
                               is.na(initial_resistance) ~ 100)) %>% 
  ungroup() %>% 
  mutate(strategy_split = interaction(strategy , initial_bg_density)) %>% droplevels() 


# Split 
sList <- split(tidy_data , tidy_data$strategy_split)

saveRDS(sList , "data/all_sim_strats_updated.rds")


