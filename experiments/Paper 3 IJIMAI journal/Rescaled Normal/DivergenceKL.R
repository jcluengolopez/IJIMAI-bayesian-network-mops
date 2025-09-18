

DivergenceKL = function(f,p,xMin,xMax)
{
  # it calculates the divergence Kullback-Leibler using the densities and numerical 
  # integration
  # INPUT: - f = original density function
  #        - p = coefficients of the fitted polynomial
  #        - (xMin,xMax) = vector with the interval where it is fitted
  # OUTPUT: - KL = divergence KL
  
  
  # delete the NAs
  p = p[!is.na(p)];
  
  # convert the coefficients of the polynomial in a function to integrate
  g = function(x,coef=p)
  {
    s = 0;
    for (i in 1:length(coef))
    {
      s = s + coef[i]*x^(i-1);
    }
    return(s)
  }
  
  # function to integrate
  #fToInt = function(x,f1=f,g1=g,coef=p)  f1(x)*log(f1(x)/g1(x,coef))
  fToInt1 = function(x,f1=f,g1=g,coef=p)  f1(x)*log(f1(x))
  fToInt2 = function(x,f1=f,g1=g,coef=p)  f1(x)*log(g1(x,coef))
  
  #KL = integrate(fToInt,xMin,xMax,subdivisions=1000L,rel.tol=1e-13,abs.tol=1e-16,
  #               stop.on.error=FALSE)$value
  KL1 = integrate(fToInt1,xMin,xMax,subdivisions=1000L,rel.tol=1e-13,abs.tol=1e-16,
                 stop.on.error=FALSE)$value
  KL2 = integrate(fToInt2,xMin,xMax,subdivisions=1000L,rel.tol=1e-13,abs.tol=1e-16,
                 stop.on.error=FALSE)$value
  
  
  return(KL1-KL2)
}
#########################################################






DivergenceKLtMoP = function(f,tMoP,xMin,xMax)
{
  # it calculates the divergence Kullback-Leibler of a tMoP density
  # INPUT: - f = original density function
  #        - tMoP = tMoP fit with polynomial and tails
  #        - (xMin,xMax) = vector with the interval where it is fitted
  # OUTPUT: - KL = divergence KL
  
  
  # extract the information of the tMoP
  poly = tMoP[[2]];        tails = tMoP[[3]];        interval = tMoP[[1]];
  
  # check if there is tail in the left
  if (tails[1]>0)   {   KL1 = DivergenceKL(f,tails[1],xMin,interval[1]);   }
  else    {   KL1 = 0   }
  
  # check if there is tail in the right
  if (tails[2]>0)   {   KL2 = DivergenceKL(f,tails[2],interval[2],xMax);   }
  else    {   KL2 = 0   }
  
  # divergence of the polynomial
  KL3 = DivergenceKL(f,poly,interval[1],interval[2]);
  
  return(KL1+KL2+KL3)
}
#########################################################




DivergenceKL2 = function(f,g,xMin,xMax)
{
  # it calculates the divergence Kullback-Leibler using the densities and numerical 
  # integration
  # INPUT: - f = original density function
  #        - g = estimated density function
  #        - (xMin,xMax) = vector with the interval where it is fitted
  # OUTPUT: - KL = divergence KL
  
  
  # function to integrate
  fToInt = function(x,f1=f,g1=g)  f1(x)*log(f1(x)/g1(x))
  
  KL = integrate(fToInt,xMin,xMax,subdivisions=1000L,rel.tol=1e-13,abs.tol=1e-16,
                 stop.on.error=FALSE)$value
  
  return(KL)
}
#########################################################






