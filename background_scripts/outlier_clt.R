# outlier_clt.R ############################################################
# This subfunction determines the auto regressive outliers in the
# residuals that are greater than a given number of std devs using the
# central limit theorem.
# 99% of the observations lie within 2.58(std_res)
# 97.5% of the observations lie within 2.24(std_res)
# Ed recommends 3(std_res)
# Return dres, rmr, & type
outlier_clt <- function(input,fig2){
  
  # initialize variables
  type <- 0 # Type of outlier detected (1=pulse, 2=trend)
  lngth <- length(YEARS)
  a <- 9
  b <- 30
  if(b>(lngth/4)){
    b <- floor(lngth/4)
    print(paste(c('Maximum outlier detection length reduced to ',toString(b),' due to low ring #'),collapse=""))
  }
  lt <- a # Length of trend
  window <- 0
  nse <- 3.29
  # nse <- 1.96 # 95 pct CI
  # nse <- 2.58 # 99 pct CI
  # nse <- 3.29 # 99.9 pct CI
  dres <- array(0,lngth) # downweighted residuals
  mr <- array(0,c(lngth,(b-a+1))) # residuals mean in window
  rmr <- array(0,lngth)
  rmu <- 0
  rshat <- 0
  
  # initialize masked to ones.
  marker <- array(0,length(YEARS))
  masked <- array(0,length(YEARS))
  
  std_res <- (var(input))^0.5 # Calculate std dev of residuals
  
  if(a<=b){
    mam <- array(NA,(b-a+1))
    mim <- array(NA,(b-a+1))
    imax <- array(NA,(b-a+1))
    imin <- array(NA,(b-a+1))
    ymn <- array(NA,(b-a+1))
    varyh <- array(NA,(b-a+1))
    for(v in a:b){
      z <- v-a+1
      for(u in 1:(lngth-v)){
        window <- input[u:(u+v-1)]
        mr[u,z] <- mean(window,na.rm=TRUE)
      }
      # Uses Tukey's bi-weight mean instead (Hoaglin 1983, Meko's code)      
      bisqmean_out <- bisqmean_dan(mr[,z])
      ymn[z] <- bisqmean_out[[1]]
      varyh[z] <- bisqmean_out[[2]]
      se <- bisqmean_out[[6]]
      mam[z] <- max(mr[,z])
      imax[z] <- which(mr[,z]==mam[z],arr.ind=TRUE)[1]
      mim[z] <- min(mr[mr[,z]!=0,z])
      imin[z] <- which(mr[mr[,z]!=0,z]==mim[z],arr.ind=TRUE)[1]
    }    
    
    poso <- (mam-ymn) / varyh # Determines number of deviations from mean of means
    nego <- (ymn-mim) / varyh
    relmam <- max(poso,na.rm=TRUE)
    rimax <- which(poso==relmam,arr.ind=TRUE)
    relmim <- max(nego,na.rm=TRUE)
    rimin <- which(nego==relmim,arr.ind=TRUE)
    
    print(paste(c('Max departure = ',toString(relmam)),collapse=""))
    print(paste(c('Min departure = ',toString(relmim)),collapse=""))
    if((poso[rimax]>=nse)){#  & (relmam>=relmim)){ # Comment & for only + outliers
      type <- 2
      lt <- rimax+a-1 # length of trend
      dres[imax[rimax]:(imax[rimax]+lt-1)] <- ymn[rimax]
      print(paste(c('Release detected in ',toString(YEARS[imax[rimax]]),' for ',toString(lt),' years'),collapse=""))
      rmu <- ymn[rimax]
      rshat <- varyh[rimax]
      rmr <- mr[,rimax]
      # print(paste(c('rmu= ',toString(rmu)),collapse=""))
      # print(paste(c('rshat= ',toString(rshat)),collapse=""))
#     }else if(nego[rimin]>=nse){ # Comment elseif for only + outliers
#       type <- 3
#       lt <- rimin+a-1
#       dres[imin[rimin]:(imin[rimin]+lt-1)] <- ymn[rimin]
#       print(paste(c('Suppression detected in ',toString(YEARS[imin[rimin]]),' for ',toString(lt),' years'),collapse=""))
#       rmu <- ymn[rimin]
#       rshat <- varyh[rimin]
#       rmr <- mr[,rimin]
#       print(paste(c('rmu= ',toString(rmu)),collapse=""))
#       print(paste(c('rshat= ',toString(rshat)),collapse=""))
    }else{
      rmu <- ymn[1]
      rshat <- varyh[1]
      rmr <- mr[,1]
      print(paste(c('No trend outliers detected up to ',toString(b),' yrs'),collapse=""))
    }
  }
  
  #if(fig2==1){ # fig=1 if you want a figure as output
  #  windows()
  #  hist(rmr)
  #  h413 = findobj(gca,'Type','patch');
  #  set(h413,'FaceColor','k')
  #  # title('\bf Histogram of Running AR Residual Means')
  #  line([rmu-nse*rshat rmu+nse*rshat],[10 10],'Color',[.6 .6 .6])
  #  plot(rmu,10,'o','Color',[.6 .6 .6])
  #  ylabel('\bf Frequency')
  #  xlabel(['\bf' int2str(lt) '\bf Yr Residual Means'])
  # }
  return(list(dres,rmr,type))
}
