# v105pn.R
# This function extracts a single tree-ring time series from 
# an input RWL file and places it in a vector for time series analysis.
# The function power transforms and removes the mean to create transformed
# residuals.  The function then detrends with an iterative neg. exponential
# fit, or if that does not fit or fails to find a solution, then a linear 
# regression with either a positive or negative slope is fit to the data.  
# Using the maximum entropy model solution otherwise known as the Burg 
# method, the autoregressive model that is the best fit for the series is 
# determined.  Using the best fit model, the function searches for 
# autoregressive outliers iteratively.  These outliers may either be pulse 
# events (1 yr) or CSTs (> minimum no. of yrs).  After the first pass,
# the outliers are removed and the series is reconstituted.  The best ar 
# order is then redetermined and the function searches for additional 
# outliers.  The # of iterations is set by the user (8 should be enough).  
# This version uses a power transformation to minimize
# the heteroscedastic nature of my time series. 'fig' is a flag that 
# specifies whether you want a figure (=1) or not (=0).  Missing years are
# set to the average of neighboring rings.  The central limit theorem 
# is used to search the residuals for trend outliers.  This version also
# uses Dave Meko's biweight mean code and currently runs with a window of 
# 9 to 30 yrs.  Estimated values fpr missing rinngs are removed in the 
# output series.  This version uses a modified Hugershoff curve with a 
# potentially nonzero asymptote to detrend + and - disturbance events. 
# It also returns the transformed standardized series.
# 
# Updated input (9/25/2014): User defined working directory (direct)
# and input file (fileN).
# 
# Function written Sep 10, 2002.
# Function last revised Jun 6, 2014.
# Function converted to R Sep 25, 2014

# Test core
#core <- "TJO180"
#core <- "TJ190A"
#core <- "CGR03b"
core <- 29
fig <- 0
iter <- 8
#fileN <- "Bishop_Example.txt"
#fileN <- "CGR_all_5b_changed.rwl"
fileN <- "EasternMA_Presto_QUSP_Final.rwl"
direct <- "/Users/tessamandra/OneDrive - Harvard University/NortheastRegionalDisturbanceProject/Sorted regional data/Eastern MA/Presto"

v105pn = function(direct,fileN,core,fig,iter){
  
  # Load packages
  library(R.matlab)
  library(dplR)
  library(MASS)
  library(robustbase)
  library(minpack.lm)
  
  # Set working directory
  setwd(direct)
  
  # Load sub-functions
  source("bisqmean_dan.R")
  source("ar_order.R")
  source("outlier_clt.R")
  source("backcast.R")
  
  # Load tree-ring data (returns var *rings*)
  # Unneccesary output suppressed by 'capture.output' & assigned to *dummy.capt*
  # NOTE: Global env variables with same name will mask "attached" variables, if they exist
  dummy.capt <- capture.output(rings <- read.rwl(fileN,format="auto"))
  col_header <- names(rings)
  
  # Find pointer to start and end of series
  ind <- which(rings[,core]>0,arr.ind=TRUE)
  sos <- ind[1]
  eos <- tail(ind,1)
  
  # Assign years (as global variable) and raw widths to respective vectors.
  YEARS <<- as.numeric(rownames(rings[sos:eos,]))
  raw <- rings[sos:eos,core]
  print(paste(c('Core: ',core),collapse=""))
  nyrs <- length(YEARS)
  print(paste(c('Total no. of measured years: ',toString(nyrs)),collapse=""))
  print(paste(c('First year is ',toString(YEARS[1])),collapse=""))
  print(paste(c('Last year is ',toString(YEARS[nyrs])),collapse=""))
  
  # Estimate missing ring widths using mean of neighboring rings
  mss <- array(0,length(raw))
  
  if(any(raw==0)){    
    m1 <- which(raw==0,arr.ind=TRUE)
    print(paste(c('Missing rings at years ',toString(YEARS[m1])),collapse=""))
    for(nm in 1:length(m1)){
      ind <- which(raw[1:m1[nm]]>0,arr.ind=TRUE)
      prior <- mean(raw[tail(ind,1)])
      ind <- which(raw[m1[nm]:length(raw)]>0,arr.ind=TRUE)
      subs <- mean(raw[ind[1]+m1[nm]-1])
      mss[m1[nm]] <- mean(c(prior,subs))
    }
    raw <- raw+mss  # Replace missing values with imputed values
    #print(mss)
  }
  
  # Power transformation.
  fdiff <- array(0,c(length(YEARS),2))
  for(x in 1:((length(YEARS)-1))){ # Calculate 1st differences
    fdiff[x,1] <- raw[x+1]
    fdiff[x,2] <- abs(raw[x+1]-raw[x])
  }
  
  s <- 1
  for(q in 1:(length(YEARS)-1)){
    if(fdiff[q,1]!=0 & fdiff[q,2]!=0){
      if(s==1){
        nz_fdiff <- fdiff[q,1:2]  # non-zero ring widths
      }else{
        nz_fdiff <- rbind(nz_fdiff,fdiff[q,1:2])  # non-zero ring widths
      }
      s=s+1
    }
  }
  log_fdiff <- log(nz_fdiff)
  
  Y <- log_fdiff[,2]
  X <- log_fdiff[,1]
  bb <- lm(Y~X)
  optimal_line <- predict(bb)
  
  optimal_pwr <- 1-bb$coefficients[[2]]
  print(paste(c('Optimal Power = ',toString(optimal_pwr)),collapse=""))
  if(optimal_pwr <= 0.05){
    transformed <- log10(raw)
    tzero <- log10(0.001)
    print('Series was log10 transformed')
  }else if(optimal_pwr>1){
    optimal_pwr <- 1
    transformed <- raw^optimal_pwr # Don't need matrix operator in R for power function
    tzero <- 0.001^optimal_pwr
    print('Series was power transformed with power = 1')
  }else{
    transformed <- raw^optimal_pwr
    print(paste(c('Series was power transformed with power = ',toString(optimal_pwr)),collapse=""))
    tzero <- 0.001^optimal_pwr
  }
  
  transm <- mean(transformed)
  # print(paste(c('tzero = ',toString(tzero)),collapse=="))
  
  # tresids <- transformed-transm # transformed residuals cannot be used with
  # iterative age detrending since half of residuals are negative.
  
  
  options(warn=-1)
  # Nonlinear detrending option.
  # Function nls employs nonlinear least squares data fitting by the
  # Gauss-Newton Method.
  # Function nlsLM employs nonlinear least squares data fitting by the
  # Levenberg-Marquardt Method.
  crashed <- rep(0,nyrs)
  wlngth <- rep(0,nyrs)
  trendtype <- 0 # Neg exp = 1, neg linear reg = 2, or pos linear reg = 3 
  minyr <- 30 # minimum number of yrs to fit to nlinfit
  if(minyr>nyrs){
    print('Insufficient # of years to fit minimum nonlinear age trend.')
  }
  b <- array(0,c(nyrs,3))
  mse <- array(NA,nyrs)
  wExist <- array(NA,(nyrs-minyr+1))
  for(i in minyr:nyrs){
    last.warning <- NULL
    beta <- c(0.5,0.1,1)
    xyrs <- c(1:i) # set years from 1 to length of series
    Y <- transformed[1:i]
    X <- xyrs[1:i]
    # Stops code from crashing because of problems fitting exp curve
    nlin.mod <- try({
      nlin.mod <- nls(Y~beta1*exp(-beta2*X)+beta3,start=list(beta1=0.5,beta2=0.1,beta3=1))
      #nlin.mod <- nlsLM(Y~beta1*exp(-beta2*X)+beta3,start=list(beta1=0.5,beta2=0.1,beta3=1))
    },silent=TRUE)
    # Check for error/warning message
    wExist[i] <- is(nlin.mod,"try-error")
    if(is(nlin.mod,"try-error")){
      crashed[i] <- 2
      #print('Exponential fit function failed, reverting to linear regression.')
    }else{
      b[i,1:3] <- summary(nlin.mod)$coefficients[1:3]  # Model coefficients
      mse[i] <- mean((predict(nlin.mod)-Y)^2,na.rm=TRUE)  # Manually calculate MSE from residuals
      crashed[i]=1
      #print(paste(c('Variance of the error term: ',toString(mse[i])),collapse=""))
    }
  }
  options(warn=1)
  
  #print(mse)
  
  # Dissallow curve to be concave up and make sure nlinfit
  # converges by making b(2) sufficiently large.
  # constant b(3) must be >=0 in original mm
  i_c <- which(crashed==1 & b[,1]>=0 & b[,2]>0.001 & b[,3]>=tzero & wExist==FALSE,arr.ind=TRUE) # & b[,2]<0.5)
  
  if(length(i_c>0)){
    mmse <- min(mse[i_c])
    imse <- which.min(mse[i_c])
  }else{
    mmse <- 0
    i_c <- 0
    imse <- 1
  }
  
   if(fig==1){ # fig==1 if you want a figure as output
     windows()
     #par(mfrow=c(3,1))
     plot(x=YEARS,y=raw,col="black",lwd=2,lty=1,type="l",ylab="Ring width (mm)",xlab="Year")
     fig1atext <- paste(c("Optimal power = ",substring(toString(optimal_pwr),1,6)),collapse="")
     mtext(fig1atext,side=3,adj=0.5,line=-1.5)
     title(main=paste(c("CST Intervention Detection on Core ",core),collapse=""))
   }
  
  if(i_c[imse]>0){
    print(paste(c('Lowest error from fit = ',substring(toString(mmse),1,6)),collapse=""))
    print(paste(c('Best age trend fit from years ',toString(YEARS[1]),' to ',toString(YEARS[i_c[imse]])),collapse=""))
    print(paste(c('Best fit extends for ',toString(i_c[imse]),' years'),collapse=""))
    best <- b[i_c[imse],]  
    trendtype <- 1
    y_exp <- best[1]*exp(-best[2]*xyrs)+best[3]
    detrended <- transformed-y_exp
    
    print('Initial Age Detrending')
    print(paste(c('Y = ',substring(toString(best[1]),1,6),'*exp(-',substring(toString(best[2]),1,6),'*x)+',
                  substring(toString(best[3]),1,6)),collapse=""))
    
    if(fig==1){ # fig==1 if you want a figure as output
      windows()
      #par(mfrow=c(3,2))
      par(mar=c(0,0,0,0),oma=c(5,5,5,5))
      plot(x=YEARS,y=transformed,col="black",lwd=2,lty=1,type="l",ylab="",xlab="",xlim=c(min(YEARS),max(YEARS)))
      lines(YEARS,y_exp)
      mtext("Year",side=1,line=2.5,adj=0.5)
      mtext("Transformed Width",side=2,line=2.5,adj=0.5)
      fig1btext <- paste(c("Y = ",substring(toString(best[1]),1,6),"*exp(-",substring(toString(best[2]),1,6),
                           "*x)+",substring(toString(best[3]),1,6)),collapse="")
      mtext(fig1btext,side=3,adj=0.5,line=-1.5)
      par(new=TRUE)
      plot(YEARS[i_c],mse[i_c],axes=FALSE,xlim=c(min(YEARS),max(YEARS)),type="l",xlab="",ylab="",col="red")
      axis(4)
      mtext("Error Term Variance",side=4,line=2.5,adj=0.5)
      windows()
      plot(YEARS,detrended,type="l",lwd=2,xlab="Year",ylab="Detrended Width")
    }
  }else{
    trendtype=2;
    xyrs <- c(1:nyrs)
    # Linear detrending option used if neg. exponential curve dissallowed.
    mod1.lm <- lm(transformed~xyrs)
    b <- mod1.lm$coefficients
    sum.lm <- summary(mod1.lm)
    stats <- c(sum.lm$r.squared,sum.lm$fstatistic,sum.lm$coefficients[[8]],sum.lm$sigma)
    if(b[2]>=0){
      trendtype <- 3 
    }
    y_lin <- b[2]*xyrs + b[1]
    detrended <- transformed - y_lin
    print('Initial Age Detrending')
    print(paste(c('Y = ',substring(toString(b[2]),1,6),' * X + ',substring(toString(b[1]),1,6)),collapse=""))
    if(fig==1){ # fig=1 if you want a figure as output
      windows()
      #par(mfrow=c(3,2))
      #par(mar=c(0,0,0,0),oma=c(5,5,5,5))
      plot(x=YEARS,y=transformed,col="black",lwd=2,lty=1,type="l",ylab="",xlab="",xlim=c(min(YEARS),max(YEARS)),
           main="Transformed Series")
      lines(YEARS,y_lin)
      mtext("Year",side=1,line=2.5,adj=0.5)
      mtext("Transformed Width",side=2,line=2.5,adj=0.5)
      fig1btext <- paste(c("Y = ",substring(toString(b[2]),1,6),"* X +",substring(toString(b[1]),1,6)),collapse="")
      mtext(fig1btext,side=3,adj=0.5,line=-1.5)
      windows()
      plot(YEARS,detrended,type="l",lwd=2,xlab="Year",ylab="Detrended Width",main="Detrended Series")
    }
  }
  # Output age detrending info
  age <- c(as.character(col_header[core]),trendtype,YEARS[i_c[imse]])

  # Plot a histogram of the data to investigate its skew.
  if(fig==1){
    windows()
    hist(detrended,xlab="Detrended Width",main="Histogram of Detrended Ring Widths")
    box()
  }

  # Initialize arrays.
  next_iter <- 1 #Switch to determine whether next iteration is needed
  St <- detrended # St will be the iterated series (standardized)
  Atr<- array(NA,length(raw)) # Age trend re-expressed in raw units
  rline <- array(NA,length(raw)) # Just the slope of the intervention
  tline <- array(NA,c(length(raw),iter)) # Slope and constant of the intervention
  outs <- array(0,c(iter,5))
  
  for(q in 1:iter){ # Iterate AR model 'iter' times to remove all outliers
    if(next_iter==1){
      bckcasted <- 0
      ar_estimate <- 0
      residuals <- 0
      area_t <- 0
      iter_i <- St # Initial values of series for ith iteration.
      print(' ')
      print(paste(c('Statistics for AR model iteration ',toString(q),':'),collapse=""))
      
      # Calculate best AR model order and return in the following order:
      # residuals (white noise) and ar model estimates
      
      ar.mod <- ar_order(St)
      ar_white <- ar.mod[[1]]
      ar_model <- ar.mod[[2]]
      
      # Use new coefficients to prewhiten ORIGINAL series without
      # downweighted originals.
      
      # Backcast for pth order years of AR model.
      # bckcasted <- backcast(detrended)
      bckcasted <- backcast(St) #I think this is correct here.
      
      for(g in ORDER:length(bckcasted)){ # g = observation year
        ar <- 0 # ar model estimate for order i, year g
        for(k in 1:ORDER){ # kth parameter of order ORDER
          if((g-ORDER)>0){ # ensure obs yr > model order
            ar <- PARAM[k]*(bckcasted[g-k])+ar
          }
        }
        if((g-ORDER)>0){ # calculate model estimate and residuals
          if(detrended[g-ORDER]==0){ # Set missing rings to ar estimate value
            print(paste(c('Missing ring at year: ',toString(YEARS[g-ORDER])),collapse=""))
            bckcasted[g] <- ar
          }
          ar_estimate[g-ORDER] <- ar
          residuals[g-ORDER] <- bckcasted[g]-ar
        }
      }

      if(fig==1){ # fig=1 if you want a figure as output
        windows()
        plot(YEARS,St,type="l",lwd=2,main=paste(c(core," iteration ",toString(q)),collapse=""),
             ylab="Transformed Width",xlab="Year",ylim=c(min(St)*0.9,max(St)*1.1))
        lines(YEARS,ar_estimate)
        legend("topright",legend=c("Standardized Series","Autoregressive Estimate"),lwd=c(2,1),bty="n")
      }
      
      # Find release outliers
      outlier.out <- outlier_clt(residuals,fig)
      downres <- outlier.out[[1]]
      mres <- outlier.out[[2]]
      otype <- outlier.out[[3]]
      
      f <- which(downres!=0,arr.ind=TRUE)
      
      if(otype==1 & length(f)!=0){ # Pulse Outlier Detected
        St[f] <- ar_estimate[f]
      }else if(otype>1 & length(f)>1){ # Trend Outlier Detected
        w <- c(1:length(f))
        slope <- lm(St[f]~w)
        print(paste(c('Constant and slope = ',substring(toString(slope$coefficients[[1]]),1,6),", ",
                      substring(toString(slope$coefficients[[2]]),1,6)),collapse=""))
        
        # Fit Hugershoff curve to remainder of series
        lngthw <- c(min(f):length(St))
        lngthwf <- c((max(f)+1):length(St))
        lngthn <- c(1:length(lngthw))
        
        lngthn <- lngthn
        opts <- nls.control(maxiter=1000)
        X <- lngthn
        Y <- St[lngthw]
        nlin.mod <- try({
          #nlin.mod <- nls(Y~(beta1*(X^beta2)*exp(-beta3*X))+beta4,
          #                start=list(beta1=0.1,beta2=0.5,beta3=0.1,beta4=-0.1),control=opts)
          nlin.mod <- nlsLM(Y~(beta1*(X^beta2)*exp(-beta3*X))+beta4,
                        start=list(beta1=0.1,beta2=0.5,beta3=0.1,beta4=-0.1),control=opts)
        
          #targets <- data.frame(Y,X)
          #model.hug <- function(beta1,beta2,beta3,beta4){
          #  (beta1*(targets$X^beta2)*exp(-beta3*targets$X))+beta4
          #}
          #par <- list(beta1=0.1,beta2=0.5,beta3=0.1,beta4=-0.1)
          #var <- list(x="Y",mean="predicted",sd=1)
          #par_lo <- list(beta1=-25,beta2=-25,beta3=-25,beta4=-25)
          #par_hi <- list(beta1=25,beta2=25,beta3=25,beta4=25)
          #nlin.mod <- anneal(model=model.hug,
          #                   par=par,
          #                   var=var,source_data=targets,max_iter=1000,pdf=dnorm,dep_var="Y",
          #                   par_lo=par_lo,
          #                   par_hi=par_hi)
          
        },silent=TRUE)
        #print(nlin.mod)
        
        if(!is(nlin.mod,"try-error")){
          bw <- summary(nlin.mod)$coefficients[1:4]
          #bw <- unlist(nlin.mod$best_pars)
          print(paste(c('Hugershoff Parameters: ',toString(bw)),collapse=""))
          ar_est <- ar_estimate[f[1]] 
          rline[lngthw] <- -bw[1]*(lngthn^bw[2])*exp(-bw[3]*lngthn)-bw[4]
        }
        
        # If nlinfit returns error, then try again with diffent initial parameters.
        if(is(nlin.mod,"try-error")){
          rline[lngthw] <- 0
          print('Default initial parameters for Hugershoff curve failed')
          print('Fitting alternate, robust initial parameters [.1 .5 .1 .1]')
          inMat <- data.frame(Y,X)
          #stop("Adjust model in code to run robust bisquare model...")
          nlin.mod <- try({
            nlin.mod <- nlrob(Y~(beta1*(X^beta2)*exp(-beta3*X))+beta4,data=inMat,
                              start=list(beta1=0.1,beta2=0.5,beta3=0.1,beta4=-0.1),
                              psi=psi.bisquare,control=opts)
          },silent=TRUE)
          print(nlin.mod)
          if(!is(nlin.mod,"try-error")){
            print(paste(c('Hugershoff Parameters: ',toString(bw)),collapse=""))
            ar_est <- ar_estimate[f[1]]
            rline[lngthw] <- -bw[1]*(lngthn^bw[2])*exp(-bw[3]*lngthn)-bw[4]
          }
        }
        
        # If nlinfit returns error, then end outlier iterations and quit.
        if(is(nlin.mod,"try-error")){
          rline[lngthw] <- 0
          print('Unable to fit Hugershoff curve')
          ar_est <- 0
          next_iter <- 0
          outs[q,1:5] <- rep(0,5)
        }
        
        if(f[1]>1){
          St[lngthw] <- rline[lngthw]+St[lngthw]+ar_est
          tline[lngthw,q] <- -rline[lngthw]
        }else if(f[1]==1){ # If trend occurs in 1st yr of series
          St[lngthw] <- rline[lngthw]+St[lngthw]
          tline[lngthw,q] <- -rline[lngthw]
        }
        #print(c(YEARS[min(f)],YEARS[max(f)],slope$coefficients[[1]],slope$coefficients[[2]],otype))
        outs[q,1:5] <- c(YEARS[min(f)],YEARS[max(f)],slope$coefficients[[1]],slope$coefficients[[2]],otype)
      }
      
      if(length(f)==0){ # Determine whether any outliers...
        next_iter <- 0 # were detected on this iteration
      }
      
      if(q==iter & length(f)>0){
        print('Need to run additional iterations to resolve series!')
      }
      
      if(fig==1){ # fig=1 if you want a figure as output
        ymin=min(c(min(iter_i),min(St)))*1.1
        ymax=max(c(max(iter_i),max(St)))*1.1
        if(length(f)>0){ # Draw detrended regression line
          if(min(f)>1){
            windows()
            plot(c(YEARS[min(f)],YEARS[max(f)]),c(ar_est,ar_est),type="l",lwd=2,lty=2,ylab="Transformed Width",
                 xlab="Year",col="grey60",ylim=c(ymin,ymax))
            # lines(YEARS[f],St[f],lwd=2,lty=2,col="grey60")
          }else{
            windows()
            plot(c(YEARS[1],YEARS[max(f)]),rep(0,2),type="l",lwd=2,lty=2,ylab="Transformed Width",xlab="Year",
                 ylim=c(ymin,ymax),col="grey60")
          }
        }else{ # Draw same line, but set to first year of series
          windows()
          plot(YEARS,rep(0,length(YEARS)),type="l",lwd=2,lty=2,col="grey60",ylab="Transformed Width",xlab="Year",
               ylim=c(ymin,ymax))
        }
        lines(YEARS,iter_i,lwd=2)
        lines(YEARS,St)
        lines(YEARS,tline[,q],lwd=2,col="grey60")
        legend("topleft",legend=c("AR Estimate","Initial Series before detrend","Standardized Series","Hugershoff Fit"),
               lty=c(2,1,1,1),lwd=c(2,2,1,2),col=c("grey60","black","black","grey60"))
        windows()
        plot(YEARS,residuals,type="l",xlab="Year",ylab="Residuals",ylim=c(min(residuals)*1.1,max(residuals*1.1)),
             main=paste(c("Iteration ",toString(q),": AR Residuals and running residual mean for ",
                          toString(length(f))," years"),collapse=""))
        abline(h=0)
        lines(YEARS,mres,col="grey60",lwd=2)
      }
    }
  }
  
  
  if(fig==1){ # fig=1 if you want a figure as output
    # Shows final iterated series in transformed units
    windows()
    #transDt <- detrended - St # transformed outlier series
    # print(paste(c('transDtarea = ',toString(sum(transDt[f]))),collapse=""))
    # plot(YEARS,detrended,type="l")
    # lines(YEARS,transDt,lwd=2)
    plot(YEARS,detrended,type="l",lwd=2,main="Outlier Series",ylab="Transformed Width",xlab="Year",
         ylim=c(min(St),max(detrended)))
    lines(YEARS,St,lwd=2,col="grey60");
    legend("topleft",legend=c('Age-detrended series','Standardized series'),col=c("black","grey60"),
           lwd=c(2,2),bty="n")
    #     windows()
    #     plot(YEARS,St,lwd=2,main="Standardized Series",xlab="Years",ylab="Transformed Width")
  }
  
  # Shows final iterated series in original units (mm presumably)
  
  # St <- St+St_pos
  
  if(trendtype==1){ # negative exponential trend
    Stt <- y_exp+St # Size trend & first detrending
    # Stt <- Stt+transm # Add transformed mean back to series
    
    
    if(optimal_pwr<=0.05){
      Str <- 10^(Stt) # Size trend in original (raw) units
      Atr <- 10^(y_exp) # Age trend in original (raw) units
    }else{
      Stt[Stt<=0] <- 0 # Set neg values to zero
      Str <- (Stt)^(1/optimal_pwr)
      Atr <- (y_exp)^(1/optimal_pwr)
    }
  }else if(trendtype==2 | trendtype==3){ # linear regression trend
    Stt <- y_lin+St # Size trend & first detrending
    # Stt <- Stt+transm # Add transformed mean back to series
    
    if(optimal_pwr<=0.05){
      Str <- 10^(Stt) # Size trend in original (raw) units
      Atr <- 10^(y_lin) # Age trend in original (raw) units
    }else{
      Stt[Stt<=0] <- 0 # Set neg values to zero
      Str <- (Stt)^(1/optimal_pwr)
      Atr <- (y_lin)^(1/optimal_pwr)
    }
  }else{
    print('Error in trend type designation')
  }
  # outs[,6] <- mean(prline,na.rm=TRUE) # mean prior to intervention removal
  # outs[,7] <- mean(arline,na.rm=TRUE) # mean after intervention removal
  
  raw[mss>0] <- NA # Remove estimated values of missing rings
  Str[mss>0] <- NA # Remove estimated values of missing rings
  Dtr <- raw-Str # Remove estimated values of missing rings
  
  if(fig==1){ # fig=1 if you want a figure as output
    windows()
    plot(YEARS,Dtr,type="l",col="grey60",lwd=2,lty=2,xlab="Year",ylab="Ring Width (mm)",ylim=c(min(Dtr),max(raw)))
    lines(YEARS,Str,lwd=2,lty=2)
    lines(YEARS,raw,lwd=2)
    #lines(YEARS,Atr)
    legend("topleft",legend=c('Disturbance index','Standardized series','Original series'),
            col=c("grey60","black","black"),lty=c(2,2,1),lwd=c(2,2,2),bty="n")
   }
  
  return(list(YEARS,transformed,detrended,St,Str,Dtr,Atr,age,outs))
  
} # End of function
