library(tidyverse)
library(ggplot2)
library(ggpubr)
library(mgcv)

source("imputation_helpers.R")



# data --------------------------------------------------------------------



modsAll <- readRDS("data/H_saturated_core_comp_nostrat_RES_MODS")
thin <- seq(10 , 300 , by = 3)
mods <- modsAll$modimp[thin]

# summarise ---------------------------------------------------------------

simCoefs <- lapply(mods , function(x) gam.mh(x , ns = 1000 , thin = 2))
saveRDS(simCoefs , "data/simulated_coeffs.rds")
simCoefs <- readRDS("data/simulated_coeffs.rds")

coefs    <- lapply(1:length(simCoefs) , function(x) data.frame(simCoefs[[x]]$bs , check.names = FALSE)) %>% 
  bind_rows(.id = "imp") %>% pivot_longer(-imp , names_to = "coef")



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
                          str_detect(coef , "rotation|crop") ~ "Rotation")) %>% 
  mutate(coef = str_remove_all(coef , ".x|factor|\\(|\\)|cult_cat.x|bs\\.|\\.")) 


coef_sum$type <- factor(coef_sum$type, 
                        levels = factor(c("Intercept" , "Density" , "Year" , 
                                          "Cultivation" , "Herbicide","Timing",
                                          "Soil"  , "Rotation")))

unique(coef_sum$coef)

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
  filter(!str_detect(coef , "rotation|FF")) %>% 
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
  facet_wrap(~type , scale = "free" , nrow = 2 , ncol=4)+
  coord_flip()


pm

pr <- coef_tidy %>%   
  filter(str_detect(coef , "rotation")) %>% 
  mutate(coef = str_remove_all(coef , "^bs\\.|rotation")) %>% 
  mutate(coef = str_replace_all(coef , "_" , "\n -> \n")) %>% 
  ggplot(aes(reorder(coef, -value_m) , value_m))+
  geom_errorbar(aes(ymin = upr90 , ymax = lwr90), width = 0)+
  geom_errorbar(aes(ymin = upr50 , ymax = lwr50), width = 0 , lwd=1.5)+
  theme_classic()+
  labs(y = "Value" , x = "Coefficient")+
  geom_hline(yintercept = 0  , lty = 3)+
  geom_point(size = 2 , pch=21 ,  fill = "black")+
  geom_point(size = 1 , pch=21 , fill = "white")+
  facet_wrap(~type , scale = "free_x" , nrow = 2)+
  theme(strip.background = element_blank())


library(patchwork)
tiff(filename="figures/fig2a.tif",height=4600,width=6200,units="px",res=800,compression="lzw", type="cairo")  
pm / pr
dev.off()


tiff(filename="figures/figS2.tif",height=4600,width=6200,units="px",res=800,compression="lzw", type="cairo")  

pf <- coef_tidy %>%   
  filter(str_detect(coef , "FF")) %>% 
  mutate(coef = str_remove_all(coef , "^bs\\.|sFF")) %>% 
  ggplot(aes(reorder(coef , -value_m) , value_m))+
  geom_errorbar(aes(ymin = upr90 , ymax = lwr90), width = 0)+
  geom_errorbar(aes(ymin = upr50 , ymax = lwr50), width = 0 , lwd=1.5)+
  theme_classic()+
  theme(strip.background = element_blank(),
        axis.text.x = element_blank())+
  labs(y = "Value" , x = "Field ID")+
  geom_hline(yintercept = 0  , lty = 3)+
  geom_point(size = 2 , pch=21 ,  fill = "black")+
  geom_point(size = 1 , pch=21 , fill = "white")

dev.off()
