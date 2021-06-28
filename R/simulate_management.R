# Simulations
library(tidyverse)
library(mgcv)

source("R/simulation_functions.R")
# models & data -----------------------------------------------------------


# imputated models
impRes     <- readRDS("data/H_saturated_markovian_nostrat_RES.rds")
impmodList <- impRes$modimp[c(41:50)]  # take final 10 imputations

# dry run
nMod <- c(41,50) 
Nt <- list(matrix(c(0.5,0.25,0.25,0,0),ncol=5) )
states <- c(1,2,3,4,5)


# simulation parameters ---------------------------------------------------


# start with a list of the parameters for a single simulation
parList_A <- list(
  rotation     = c("wheat" , "bar","wheat","bar","wheat"),
  d_season = c("winter","spring","winter","spring","winter"),
  soil         = rep("S4" , 5),
  nHerb        = c(10,0,10,0,10),
  nGWS         = c(0,0,0,0,0),
  aGly         = c(0,0,0,0,0),
  cultivation  = c("conventional", "conventional","conventional","conventional","conventional"),
  c_diff       = c(0,0,0,0,0),
  d_diff       = c(5,5,5,5,5)
)

# start with a list of the parameters for a single simulation
parList_B <- list(
  rotation     = c("wheat" , "wheat","wheat","wheat","wheat"),
  d_season = c("winter","winter","winter","winter","winter"),
  soil         = rep("S4" , 5),
  nHerb        = c(10,0,10,0,10),
  nGWS         = c(0,0,0,0,0),
  aGly         = c(0,0,0,0,0),
  cultivation  = c("conventional", "conventional","conventional","conventional","conventional"),
  c_diff       = c(0,0,0,0,0),
  d_diff       = c(-5,-5,-5,-5,-5)
)


parList <- list(parList_A , parList_B)
# Run simulations ---------------------------------------------------------


simList <- list()
for(i in seq_along(parList)){
  resList     <- unpack_impMod(impRes , nMod)
  tf_matrices <- lapply(resList ,construct_Tf, parList[[i]])
  simList[[i]]     <- lapply(tf_matrices ,project_sim, Nt) %>% bind_rows(.id = "imp") %>% mutate(imp = as.numeric(imp))
}


# Summarise & plot --------------------------------------------------------

# Summarise data
projRes     <- simList  %>% bind_rows(.id = "strategy") 
projRes_sum <- projRes %>% 
  group_by(iter,strategy) %>%
  summarise(mss = median(ms) , 
            upr90 = quantile(ms , 0.90), 
            lwr90 = quantile(ms , 0.1), 
            lwr50 = quantile(ms , 0.25), 
            upr50 = quantile(ms , 0.75)) %>% 
  rename(ms = mss)


projRes %>% 
  filter(field == "DN40_3PX~glebe.24") %>%                                                               # Choose field
  mutate(fimp = interaction(field , imp,strategy)) %>% 
  ggplot(aes(iter, ms))+
  geom_line(alpha=0.8,aes(group=fimp,colour=factor(strategy)) , show.legend = FALSE,lty=3)+              #
  geom_line(data = projRes_sum , aes(iter,ms,group=strategy))+
  geom_point(data = projRes_sum , aes(iter,ms),size=2)+
  geom_point(data = projRes_sum , aes(iter,ms , fill = strategy ), pch = 21 ,size=3,colour="black")+ 
  theme_classic()+
  labs(x = "Step" , y = "Mean density state")





