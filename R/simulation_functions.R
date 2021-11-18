# functions  ------------------------------------------------------
unpack_parList <- function(parList){
  
  # Parlist = list of parameters for simulation
  
  # Number of transitions
  rlen <- length(parList[[2]])-1
  
  datList <- list()
  for(i in seq_len(rlen)){
    
    datList[[i]] <- data.frame(
      # Constant variables
      dens_t1         = factor(1:5),
      
      crop_t1        = parList$rotation[i],
      crop_t2        = parList$rotation[i+1],
      
      
      cult_cat_t2     = parList$cultivation[i],
      soil_group_t1   = parList$soil_group[i],
      
      spray_days_gw_t2 = parList$spray_days[i+1],
      # n_prod_gw_t2      = parList$n_products[i+1],
      a_gly_t2        = parList$n_glyphosate[i+1],
      
      d_season_t2     = parList$d_season[i+1],
      d_diff_t2       = parList$d_diff[i+1],
      c_diff_t2       = parList$c_diff[i+1],
      
      mean_mort_t1    = parList$mean_mort[i],
      
      transition_year = "2015_2016",
      
      weights         = rep(1,5)
    ) 
    
    
  }
  
  return(datList)
  
}              # function to unpack parList into chunks of new data
unpack_impMod  <- function(impRes , impData , thin){
  # impRes = imputation models
  # impData = imputation data
  # thin = sequence to thin by
  
  mods <- impRes$modimp[thin]
  data <- impData$all_data[thin]
  
  outList <- list()
  for(i in seq_along(mods)){
    
    modObj <-mods[[i]]     # Models
    modDat <- data[[i]]   # Data fit to model
    
    # Scale predictors
    toScale   <- modDat %>% select(matches("gly_t2|diff_t2|spray_days_gw_t2|mort_t1|^S[0-9]_t1")) %>% colnames()
    scaleAttr <- attributes(scale(modDat[toScale] , scale = TRUE , center = TRUE)) 
    
    outList[[i]] <- list(modObj  = modObj , scaleAttr = scaleAttr)
    
  }
  
  return(outList)
}

# function to extract models & data scaling attributes 
construct_Tf   <- function(modList, parList){
  # fitList = # list containing
  # - $modObj    = model fit to data
  # - $scaleAttr = scaling attributes from data fitted to modObj 
  
  # parList = List of management parameters for a single transition.
  
  # Get models , scaling attributes, & parameter lists
  modObj          <- modList$modObj
  scaleAttr       <- modList$scaleAttr
  transition_list <- unpack_parList(parList) # Field-level data for each step of the rotation
  
  
  projection_list <- list()
  for(i in seq_along(transition_list)){
    
    td      <- transition_list[[i]]    # Get step i
    fields  <- unique(modObj$model$FF)[1] %>% droplevels() # Get field 
    newData <- expand_grid(td,fields) %>% mutate(FF = fields) # Expand data
    newData[scaleAttr$dimnames[[2]]] <- scale(newData[scaleAttr$dimnames[[2]]], scaleAttr$`scaled:center`, scaleAttr$`scaled:scale`) # scale to same attributes as fitted data
    
    projection_list[[i]] <- predict(object = modObj , newdata = newData , type = "response" , exclude ="s(FF)") %>% t() # Predict.
  }
  
  return(projection_list)
  
}      # build field level matrices
project_matrix <- function(matrix,Nt){
  # Matrix = matrix to be projected 
  # Nt = matrix of initial ds distributions by row
  Nt1 <- apply(Nt,1, function(x) matrix %*% x  )
  return(t(Nt1))
}            # project a single matrix a single time-step
mean_ds        <- function(ds,states){
  # ds = counts of how many observations in each density state OR the density state distribution as a proportion
  # states = vector of values for each state 
  
  dd <- rep(0,5) 
  dd[1:length(ds)] <- ds # account for missing 0s in ds
  fx <- dd/sum(dd)# get proportions
  ms <- fx%*%states # get mean density state
  vs <- fx%*%((states-c(ms)))^2 # get variance associated with distribution
  return(data.frame(ms = ms , vs = vs))
  
}            # calculate mean density state 
p_inc_ds       <- function( ds , target_state = 3 , states = 1:5){
  sum(ds[states[states  > target_state]])
  
} 
p_dec_ds       <- function( ds , target_state = 3 , states = 1:5){
  sum(ds[states[states  < target_state]])
  
} 
project_sim    <- function(projection_list , Nt){
  # projection_list = list of matrices for each transition of the projection
  # Nt = 1 row matrix of starting state probabilities (as a single element list)
  
  
  p_len  <- length(projection_list)   # length of projection
  

  
  res_list <- list()
  
  res_list[[1]] <- data.frame(pinc = p_inc_ds(Nt,  target_state = 3 , states),
                              pdec = p_dec_ds(Nt , target_state = 3 , states),
                              mds  = mean_ds(Nt , states), 
                              a = Nt[1] , l = Nt[2] , m = Nt[3] , h = Nt[4] , v = Nt[5])
  
  for(i in seq_len(p_len)+1){
    mSet          <- projection_list[[i-1]]  # matrix 
    Nt            <- project_matrix(mSet , Nt) # project field-level matrices
    res_list[[i]] <- data.frame(pinc = p_inc_ds(Nt,  target_state = 3 , states), # summarise to mean density state
                                pdec = p_dec_ds(Nt , target_state = 3 , states),
                                mds  = mean_ds(Nt , states) , 
                                a = Nt[1] , l = Nt[2] , m = Nt[3] , h = Nt[4] , v = Nt[5]) # summarise to mean density state
  }
  
  out   <- bind_rows(res_list , .id = "iter")  
  
} # project entire simulation from starting state Nt 
simulate_strategy <- function(parList  , impRes , impData ,  thin){
  simList <- list()
  for(i in seq_along(parList)){
    print(i)
    Nt <- parList[[i]]$Nt
    resList     <- unpack_impMod(impRes = impRes, impData = impData ,  thin = thin)
    tf_matrices <- lapply(resList ,construct_Tf, parList[[i]])
    simList[[i]]     <- lapply(tf_matrices , project_sim , Nt) %>%
      bind_rows(.id = "imp") %>% mutate(imp = as.numeric(imp) , strategy = parList[[i]]$strategy)
  }
  
  return(simList)
}
build_strategy    <- function(rotation = c("wheat_wheat" ,"wheat_wheat") , d_season = "winter",soil_group = "S5",cultivation = "conventional",spray_days = 0 , n_glyphosate = 0 , mean_mort = 100, c_diff = 0 , d_diff = 0, nstep = 3 , name){
  list(
    rotation     = rotation,
    d_season     = d_season,
    soil_group         = soil_group,
    cultivation  = cultivation,
    spray_days        = spray_days,
    n_glyphosate         = n_glyphosate ,
    mean_mort = mean_mort,
    c_diff       = c_diff,
    d_diff       = d_diff ,
    strategy         = name)
  
  
}



# functions using resimulted coefficients ---------------------------------

construct_Tf_rs      <- function(modList , Beta , parList){
  
  # Model = model object from imputation
  # beta = resimulated coefficients
  # parList =  dataframe of parameters for single transition
  
  
  # Set up model matrix for simulation
  
  # Get models , scaling attributes, & parameter lists
  modObj          <- modList$modObj
  scaleAttr       <- modList$scaleAttr
  transition_list <- unpack_parList(parList) # Field-level data for each step of the rotation
  
  Beta      <- Beta$bs                           # trim resim diags
  dlevels   <- modObj$xlevels                     # Levels for categorical predictors
  
  projection_list <- list()
  for(i in seq_along(transition_list)){
    fseq      <- which(sapply(transition_list[[i]],is.factor))  # positions in data
    transition_list[[i]][,fseq] <-   lapply(seq_along(fseq)  , function(x) {transition_list[[i]][,fseq[x]] <- factor(transition_list[[i]][,fseq[x]] , levels = dlevels[[x]])}) # relevel
    
    newData <- transition_list[[i]]
    newData[scaleAttr$dimnames[[2]]] <- scale(newData[scaleAttr$dimnames[[2]]], scaleAttr$`scaled:center`, scaleAttr$`scaled:scale`) # scale to same attributes as fitted data
    
    
    # Get model matrix & coefficients
    form <- formula(str_remove_all(modObj$formula , 'dens_t2|s\\(FF, bs = \"re\"\\) \\+|factor\\(|\\)')) # trim uneeded strings
    X    <- model.matrix(form , newData) # get Xvars matrix
    Beta <- Beta[,1:ncol(X)]       # trim field intercepts
    
    # Calculate transition probabilities
    cp <-  modObj$family$getTheta(TRUE)
    eta   <- sapply(1:5 , function(x) Beta %*% X[x,])  # Linear predictor for resim coefficients - row = coef vector , col = source state / new data
    
    
    thetaL <- array(dim = c(nrow(eta) , 5,5))  # Array of transition probabilities , row = coef vector , col = destination state , slice = destination state
    for(j in 1:5){
      theta <- matrix(0 , ncol = 5 , nrow = nrow(Beta))
      theta[,1] <- 1- plogis(eta[,i] - cp[1]) 
      for(k in 2:4){
        theta[,k] <- plogis(eta[,i] - cp[k-1]) - plogis(eta[,i] - cp[k])
      }
      theta[,5] <- plogis(eta[,i] - cp[4]) 
      thetaL[,,j] <- theta
    }
    
    # thetaL[i,,] is the transition matrix for a single vector of resimulated coefficients (i). 
    
    # Transition matrices for each resimulated coef vector
    projection_list[[i]] <- lapply(1:nrow(thetaL) , function(x) thetaL[x,,])
    
  }
  return(projection_list)
}
project_sim_rs       <- function(projection_list , Nt){
  # projection_list = list of matrices for each transition of the projection
  # Nt = 1 row matrix of starting state probabilities (as a single element list)
  
  
  p_len  <- length(projection_list)   # length of projection
  
  res_list <- list()
  Nt <- rep(Nt , length(projection_list[[1]]))
  for(i in seq_len(p_len)){
    mSet          <- projection_list[[i]]  # matrix 
    Nt            <- lapply(1:length(mSet) , function(x) project_matrix(mSet[[x]], Nt[[x]]))  # project field-level matrices
    
    
    res_list[[i]] <-  Nt %>% map(~data.frame(pinc = p_inc_ds(.,target_state = 3 , states),
                                             pdec = p_dec_ds(. , target_state = 3 , states),
                                             mds  = mean_ds(. , states))) %>% 
      bind_rows(.id="parsim")
    
    
    
  }
  
  out   <- bind_rows(res_list , .id = "iter")  
  
} # project entire simulation from starting state Nt 
simulate_strategy_rs <- function(parList_comp  , impRes , betaList ,   thin , Nt){
  simList <- list()
  
  for(i in seq_along(parList_comp)){
    print(noquote(paste("Simulating strategy" , i , " of " , length(parList_comp),"...")))
    resList     <- odod(impRes , impData ,  thin)
    tf_matrices <- mapply(construct_Tf_rs , modList = resList , Beta = betaList , parList = list(parList_comp[[i]]) , SIMPLIFY = FALSE)
    simList[[i]]     <- lapply(tf_matrices ,project_sim_rs, Nt) %>%
      bind_rows(.id = "imp") %>%
      mutate(imp = as.numeric(imp) , strategy = parList_comp[[i]]$strategy)
  }
  
  return(simList)
}


# summary & plotting ------------------------------------------------------

sum_fun   <- function(simList){
  projRes_sum <- simList %>% 
    filter(iter == max(iter)) %>% 
    pivot_longer(pinc:mds.ms , names_to = "metric") %>% 
    group_by(iter , strategy,metric) %>%
    summarise(mss = mean( value),
              ci = 1.96 * (sd(value) / sqrt(n())),
              upr90 = quantile( value , 0.90), 
              lwr90 = quantile( value , 0.1), 
              lwr50 = quantile( value , 0.25), 
              upr50 = quantile( value , 0.75)) %>% 
    rename( ms = mss) 
  return(projRes_sum)
}
nice_plot <- function(p){
  dev.off()
  gp <- ggplotGrob(p)
  facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]
  x.var <- sapply(ggplot_build(p)$layout$panel_scales_x,
                  function(l) length(l$range$range))
  gp$widths[facet.columns] <- gp$widths[facet.columns] * x.var
  grid::grid.draw(gp)
} # squishes ggplot facets


