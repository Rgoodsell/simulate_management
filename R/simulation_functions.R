# functions  ------------------------------------------------------

unpack_parList <- function(parList){
  
  # Parlist = list of parameters for simulation
  
  # Number of transitions
  rlen <- length(parList[[1]])-1
  
  datList <- list()
  for(i in seq_len(rlen)){
    
    datList[[i]] <- data.frame(
                        # Constant variables
                        dens_t1         = factor(1:5),
                        transition_year = "2015_2016",
                        weights         = rep(1,5),
                        
                        # Year 1 variables
                        crop_t1         = parList$rotation[i],
                        d_season_t1     = parList$d_season[i],
                        hBC_app_t1      = parList$nHerb[i],
                        gw_spec_t1      = parList$nGWS[i],
                        cult_cat_t1     = parList$cultivation[i],
                        c_diff_t1       = parList$c_diff[i],
                        d_diff_t1       = parList$d_diff[i], 
                        a_gly_t1        = parList$aGly[i],
                        soil_group_t1   = parList$soil[i],
                        
                        # Year 2 variables
                        crop_t2         = parList$rotation[i+1],
                        d_season_t2     = parList$d_season[i+1],
                        hBC_app_t2      = parList$nHerb[i+1],
                        gw_spec_t2      = parList$nGWS[i+1],
                        cult_cat_t2     = parList$cultivation[i+1],
                        c_diff_t2       = parList$c_diff[i+1],
                        a_gly_t2        = parList$aGly[i+1],
                        d_diff_t2       = parList$d_diff[i+1]) 
      
    
  }
  
  return(datList)
  
}              # function to unpack parList into chunks of new data
unpack_impMod  <- function(impRes, nMod){
  # impRes = imputation results
  # nMod   = range of models to take from results 
  
  outList <- list()
  
  for(i in seq(nMod[1] , nMod[2] , by = 1)){
    modObj <- impRes$modimp[[i]]     # Models
    modDat <- impRes$all_data[[i]]   # Data fit to model
  
    # Scale predictors
    toScale   <- modDat %>% select(matches("gly|diff|hBC|gw")) %>% colnames()
    scaleAttr <- attributes(scale(modDat[toScale] , scale = TRUE , center = TRUE)) 
    
    outList[[length(impRes$modimp)-(i-1)]] <- list(modObj  = modObj , scaleAttr = scaleAttr)
    
  
    }
  
  return(outList)
  
  } # function to extract models & data scaling attributes 
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
    fields  <- unique(modObj$model$FF) # Get all fields
    newData <- expand_grid(td, fields) %>% mutate(FF = fields)  # Expand data
    newData[scaleAttr$dimnames[[2]]] <- scale(newData[scaleAttr$dimnames[[2]]], scaleAttr$`scaled:center`, scaleAttr$`scaled:scale`) # scale to same attributes as fitted data
    
    projection_list[[i]] <- predict(object = modObj , newdata = newData , type = "response") %>% # Predict.
                            split(.,newData$FF) %>%                                              # Split by field.
                            map(~matrix(.,nrow=5,byrow = TRUE))                                  # Build transition matrices.
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
project_sim    <- function(projection_list , Nt){
  # projection_list = list of field level matrices for each transition of the projection
  # Nt = 1 row matrix of starting state probabilities (as a single element list)
  
  # get lists of inital density states
  initC <- mean_ds(Nt[[1]] , states = states) %>%
            slice(rep(1:n(), times = length(projection_list[[1]])))  %>% 
            mutate(field = names(projection_list[[1]]))
  
  p_len  <- length(projection_list)   # length of projection

  res_list <- list()
  for(i in seq_len(p_len)){
    mSet          <- projection_list[[i]]  # matrix set for projection 
    Nt            <- mapply(project_matrix , matrix = mSet , Nt = Nt,SIMPLIFY = FALSE) # project field-level matrices
    res_list[[i]] <- Nt %>%  map(~mean_ds(.,states)) %>% do.call(rbind,.) %>% rownames_to_column(var = "field") # summarise to mean density state
  }
  
  out   <- bind_rows(initC ,res_list , .id = "iter")  
  
} # project entire simulation from starting state Nt 

