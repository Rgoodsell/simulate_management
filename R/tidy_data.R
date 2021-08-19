
# packages ----------------------------------------------------------------

library(tidyverse)
library(janitor)
library(lubridate)

source("R/simulation_functions.R")

# data --------------------------------------------------------------------

my_data <- read_csv("data/mgmtdata_strats.csv") %>% clean_names() 

# tidy --------------------------------------------------------------------


tidy_data <-   my_data %>%
              filter(strategy == "FF") %>% 
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
                                   "maize" = "cover") , 
                       soil_regroup = recode(soil_group , 
                                           "Denchworth" = "S7" , 
                                           "Beccles" = "S5" , 
                                           "Dunnington Heath" = "S5",
                                           "Evesham" = "S4" , 
                                           "Foggathorpe" = "S7" , 
                                           "Holderness" = "S7" , 
                                           "Rivington" = "S5" , 
                                           "Wickham" = "S7" , 
                                           "Wisbech" = "S8")) %>% 
  mutate(d_season = factor(case_when(d_date < 183 ~ "spring" , d_date > 183 ~ "winter")),
         d_date = as.numeric(d_date) , 
         c_date = as.numeric(c_date) ,
         h_date = as.numeric(h_date),
         a_gly = as.numeric(a_gly),
        h_bc_app = as.numeric(h_bc_app)) %>% 
  mutate(c_date = case_when(cult_cat == "none" & d_season == "spring" ~ 251,
                            cult_cat == "none" & d_season == "winter" ~ 240,
                            TRUE ~ c_date)) %>% 
  mutate(d_diff = case_when(d_season == "spring" ~ d_date - 82, 
                            d_season == "winter" ~ d_date - 273),
         c_diff = case_when(d_season == "spring" ~ c_date - 251, 
                            d_season == "winter" ~ c_date - 240)) %>% 
  ungroup()



# Split 
sList <- split(tidy_data , tidy_data$strategy_initcond)

saveRDS(sList , "data/all_sim_strats.rds")


