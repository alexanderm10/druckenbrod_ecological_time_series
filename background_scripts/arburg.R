arburg = function(x,p){
  #   ARBURG   AR parameter estimation via Burg method.
  #   A = ARBURG(X,ORDER) returns the polynomial A corresponding to the AR
  #   parametric signal model estimate of vector X using Burg's method.
  #   ORDER is the model order of the AR system.
  #
  #   [A,E] = ARBURG(...) returns the final prediction error E (the variance
  #   estimate of the white noise input to the AR model).
  #
  #   [A,E,K] = ARBURG(...) returns the vector K of reflection 
  #   coefficients (parcor coefficients).
  #
  #   See also PBURG, ARMCOV, ARCOV, ARYULE, LPC, PRONY.
  
  #   Ref: S. Kay, MODERN SPECTRAL ESTIMATION,
  #              Prentice-Hall, 1988, Chapter 7
  #        S. Orfanidis, OPTIMUM SIGNAL PROCESSING, 2nd Ed.
  #              Macmillan, 1988, Chapter 5
  
  #   Author(s): D. Orofino and R. Losada
  #   Copyright 1988-2000 The MathWorks, Inc.
  #   $Revision: 1.10 $  $Date: 2000/06/09 22:03:11 $
  
  #error(nargchk(2,2,nargin))
  
  if(length(x)==0 | length(dim(x))>1){
    stop('X must be a vector with length greater than twice the model order.')
  }else if(length(p)==0 | p != round(p)){
    stop('Model order must be an integer.')
  }
  
  #if(issparse(x)){
  #  stop('Input signal cannot be sparse.')
  #}
  
  #x  = x(:);
  N  <- length(x)
  
  # Initialization
  ef <- x
  eb <- x
  a <- 1
  
  # Initial error
  E <- (x*x)/N
  
  k <- c()
  for(m in 1:p){
    # Calculate the next order reflection (parcor) coefficient
    efp <- ef[2:length(ef)]
    ebp <- eb[1:(length(eb)-1)]
    num <- -2*(t(ebp) %*% efp)
    den <- (t(efp) %*% efp)+(t(ebp) %*% ebp)
    
    k[m] <- num/den
    
    # Update the forward and backward prediction errors
    ef <- efp + k[m]*ebp
    eb <- ebp + k[m]*efp
    
    # Update the AR coeff.
    a <- c(a,0) + k[m]*c(0,Conj(rev(a)))
    
    # Update the prediction error
    E[m+1] = (1 - (k[m]*k[m]))*E[m]
  }

  varargout1 <- a[-1]
  varargout2 <- E[length(E)]
  varargout3 <- k

  return(list(varargout1,varargout2,varargout3))
}
