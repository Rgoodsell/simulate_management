library(tidyverse)
library(ggplot2)
library(ggpubr)
library(mgcv)
library(patchwork)

# data --------------------------------------------------------------------


modsAll <- readRDS("data/alexa_regularised_cropping_nostrat_RES.rds_MODS")
dataAll <- readRDS("data/alexa_regularised_cropping_nostrat_RES.rds_DATA")
thin <- seq(10 , 120 , by = 3)
mods <- modsAll$modimp[thin]

# simulate ----------------------------------------------------------------

simCoefs <- lapply(mods , function(x) gam.mh(x , ns = 100 , thin = 1))
saveRDS(simCoefs , "data/simulated_slim_coeffs.rds")

# plot --------------------------------------------------------------------

simCoefs <- readRDS("data/simulated_slim_coeffs.rds")

coefs    <- lapply(1:length(simCoefs) , function(x) data.frame(simCoefs[[x]]$bs , check.names = FALSE)) %>% bind_rows(.id = "imp")

# set levels for regularised cropping coefficients to rename later
crop_t1_levels  <- levels(dataAll$all_data[[1]]$crop_t1)
crop_t2_levels  <- levels(dataAll$all_data[[1]]$crop_t2)
rotation_levels <- levels(dataAll$all_data[[1]]$rotation)

names(coefs)[str_detect(names(coefs) , "crop_t1")] <- paste0(crop_t1_levels , "_crop_t1")
names(coefs)[str_detect(names(coefs) , "crop_t2")] <- paste0(crop_t2_levels , "_crop_t2")
names(coefs)[str_detect(names(coefs) , "rotation")] <- paste0(rotation_levels , "_rotation")


coefs <- coefs %>% pivot_longer(-imp , names_to = "coef")

coef_sum <- coefs %>%  bind_rows( .id = "iter") %>% 
  filter(!str_detect(coef , "accept")) %>% 
  group_by(coef) %>% 
  summarise( value_m = median(value) , 
             upr90 = quantile(value , 0.95),
             lwr90 = quantile(value , 0.05),
             upr50 = quantile(value , 0.75),
             lwr50 = quantile(value , 0.25)) %>% 
  mutate(type = case_when(str_detect(coef , "Intercept") ~ "Intercept",
                          str_detect(coef , "dens") ~ "Density" , 
                          str_detect(coef , "cult") ~ "Cultivation" , 
                          str_detect(coef , "diff|season") ~ "Timing" , 
                          str_detect(coef , "hBC|gly|gw|mort") ~ "Herbicide" , 
                          str_detect(coef , "soil|pc") ~ "Soil" , 
                          str_detect(coef , "transition_year") ~ "Year", 
                          str_detect(coef , "rotation") ~ "Rotation",
                          str_detect(coef , "crop_t1") ~ "Crop_t1",
                          str_detect(coef , "crop_t2") ~ "Crop_t2")) %>% 
  mutate(coef = str_remove_all(coef , ".x|factor|\\(|\\)|cult_cat.x|bs\\.|\\.")) 


coef_sum$type <- factor(coef_sum$type, 
                        levels = factor(c("Intercept" , "Density" , "Year" , 
                                          "Cultivation" , "Herbicide","Timing",
                                          "Soil"  , "Rotation" , "Crop_t1" , "Crop_t2")))

# pretty coef names
coef_tidy <- coef_sum %>% mutate(
  coef =   recode(coef  , 
                  "Intercept" = "Intercept" , 
                  "dens_t12" = "2",
                  "dens_t13" = "3",
                  "dens_t14" = "4",
                  "dens_t15" = "5",
                  "transition_year2015_2016" = "2015 \n -> \n 2016",
                  "cult_cat_t2inversion" = "Inv.",
                  "cult_cat_t2none"      = "Min-till",
                  "cult_cat_t2subsoil"   = "Subs.",
                  "cult_cat_t2surface"   = "Surf.",
                  "soil_group_t1S3"      = "Litho.",
                  "soil_group_t1S5"      = "Brown",
                  "soil_group_t1S8"      = "Gw. \n gley",
                  "soil_group_t1S7"      = "Sw. \n gley", 
                  "gw_spec_t2"  = "Herb.",
                  "a_gly_t2"  = "A. gly",
                  "d_season_t2winter"   = "D.season",
                  "c_diff_t2" = "C.date",
                  "d_diff_t2" = "D.date",
                  "d_season_t2winter:d_diff_t2" = "D.s \n * \n D.date",
                  "d_season_t2winter:c_diff_t2" = "D.s \n * \n C.date")) 


# plot --------------------------------------------------------------------

pm <- coef_tidy %>%   
  filter(!str_detect(coef , "crop|FF|rotation")) %>% 
  mutate(coef = str_remove(coef , "^bs\\.")) %>% 
  ggplot(aes(coef , value_m))+
  geom_errorbar(aes(ymin = upr90 , ymax = lwr90), width = 0)+
  geom_errorbar(aes(ymin = upr50 , ymax = lwr50), width = 0 , lwd=1.5)+
  theme_classic()+
  theme(strip.background = element_blank())+
  labs(y = "Value" , x = "Coefficient")+
  geom_hline(yintercept = 0  , lty = 3)+
  geom_point(size = 2 , pch=21 ,  fill = "black")+
  geom_point(size = 1 , pch=21 , fill = "white")+
  facet_wrap(~type , scale = "free",ncol=1)+
  coord_flip()

pr <- coef_tidy %>%   
  filter(str_detect(coef , "crop|rotation")) %>% 
  mutate(coef = str_remove_all(coef , "_rotation|_t1|_t2|_crop")) %>% 
  ggplot(aes(coef , value_m))+
  geom_errorbar(aes(ymin = upr90 , ymax = lwr90), width = 0)+
  geom_errorbar(aes(ymin = upr50 , ymax = lwr50), width = 0 , lwd=1.5)+
  theme_classic()+
  labs(y = "Value" , x = "Coefficient")+
  geom_hline(yintercept = 0  , lty = 3)+
  geom_point(size = 2 , pch=21 ,  fill = "black")+
  geom_point(size = 1 , pch=21 , fill = "white")+
  facet_wrap(~type , scale = "free_x" , nrow = 2)+
  theme(strip.background = element_blank())+
  facet_wrap(~type , scales = "free")+
  coord_flip()

#plot 
pm 

dev.new()
pr



