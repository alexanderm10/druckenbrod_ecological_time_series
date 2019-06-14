# backcast.R  #############################################################################
# This subfunction estimates the first elements of a series for which the   
# residuals could not be calculated owing to the use of ar modeling.       
backcast <- function(seriesb){
  
  # Invert time series for backcasting.
  flipped <- rev(seriesb)
  
  # Add in backcasted AR estimates as new values at end of inverted series.
  for(g in (length(YEARS)+1):(length(YEARS)+ORDER)){ # g = backcasted years
    ar <- 0 # ar model estimate for order i, year g
    for(k in 1:ORDER){ # kth parameter of order ORDER
      ar <- PARAM[k]*(flipped[g-k])+ar
    }
    flipped[g] <- ar
  }
  
  # Re-invert series and return as output.
  bckcasted <- rev(flipped)
  # print('Backcasted Values:')
  # for(h in ORDER:-1:1){
  #     print(paste(c('Year -',toString(h),': ',toString(bckcasted[h])),collapse=""))
  # }
  return(bckcasted)
}
