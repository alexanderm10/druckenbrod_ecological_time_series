# ar_order.R ################################################################
# This subfunction determines the
# autoregressive parameters for the best model order as calculated using
# AIC criteria. The function returns the residuals, and AR model
# estimate of the best order found with AIC criteria.
# Returns out_res and out_est
ar_order = function(series){
  
  source("arburg.R")
  
  # Calculate Autoregressive parameters, residuals, and variance for orders 1 through 10.
  ar_param <- array(0,c(10,10))
  for(ar.order in 1:10){
    ar_param[ar.order,1:ar.order] <- -arburg(series,ar.order)[[1]]
  }
  
  residuals <- array(0,c(length(YEARS),10))
  ar_estimate <- array(NA,c(length(YEARS),10))
  resid_var <- c()
  #mod.burg <- ar.burg(series,order.max=10,aic=TRUE,var.method=2)
  for(i in 1:10){
    for(g in 1:length(YEARS)){ # g = observation year
      ar <- 0 # ar model estimate for order i, year g
      for(k in 1:length(ar_param[1,])){ # kth parameter of order i
        if((g-k)>0){# & !is.na(ar_param[ar.order,k])){ # ensure obs yr > model order
          ar <- (ar_param[i,k]*(series[g-k]))+ar
        }
      }
      if((g-i)>0){ # calculate residuals
        residuals[g,i] <- (series[g]-ar)
        ar_estimate[g,i] <- ar
      }
    }
    # Calculate the total variance of the residuals by model order
    resid_var[i] <- var(residuals[,i],na.rm=TRUE)   
  }
    
  # Calculate variance of the residuals of a particular AR order model.
  # Reference Box & Jenkins & Reinsel 1994 pp. 200-201.
  # Using Akaike Information Criteria (A-key-key)
  # Equation now uses natural log and simply 'n' in the denominator.
  # t+1 (or # of params+1) is a penalty factor for estimating the mean.
  aic <- c()
  for(t in 1:length(resid_var)){
    aic[t] <- log(resid_var[t])+(2*(t+1))/length(YEARS)
  }
  # Find the first minimum AIC order (ie first saw-tooth-shaped dip).
  # If AIC values monotonically decrease, set best_order=9.
  best_order <- 0
  for(s in 2:length(aic)){
    if((aic[s]>=aic[s-1]) & (best_order==0)){
      best_order <- s-1
    }
  }
  
  
  # print('AIC Values:')
  # print(toString(aic))
  
  #best_order <- mod.burg$order
  
  if(best_order==0){
    best_order <- 9
  }
  
  # print('AIC Values:')
  # print(toString(aic))
  
  # NOTE: Assigns to Global Environment
  assign("ORDER",best_order,envir=.GlobalEnv)
  assign("PARAM",ar_param[best_order,1:best_order],envir=.GlobalEnv)
  
  # print(paste(c('Best order = ',toString(ORDER)),collapse=""))
  print(paste(c('AR Model Parameters: ',toString(PARAM)),collapse=""))
  
  out_res <- residuals[,best_order]
  out_est <- series-out_res
  return(list(out_res,out_est))
}
