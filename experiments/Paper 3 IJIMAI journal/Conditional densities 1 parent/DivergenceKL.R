
library('pracma')


DivKLpolyXY = function(fX,fYX,p,xMin,xMax,yMin,yMax)
{
  # it calculates the conditional divergence Kullback-Leibler from a polynomical fitting
  # of a variable X with parent Y using the densities and numerical integration
  # INPUT: - fX = original marginal density function of the parent   f(x)
  #        - fYX = original conditional density function   f(y|x)
  #        - p = coefficients of the fitted polynomial
  #        - (xMin,xMax) = vector with the definition interval of the parent X
  #        - (yMin,yMax) = vector with the definition interval of the child Y
  # OUTPUT: - KL = divergence KL
  
  
  # delete the NAs
  p = p[!is.na(p)];
  
  # convert the coefficients of the polynomial in a function to integrate
  g = function(y,coef=p)
  {
    s = 0;
    for (i in 1:length(coef))     s = s + coef[i]*y^(i-1)
    return(s)
  }
  
  # functions for both integrals (X is the parent, Y is the child)
  #  f1X(x) * f1YX(y,x) * log(f1YX(y,x) / g1(y,coef))
  fToInt1 = function(x,y,f1X=fX,f1YX=fYX)      f1X(x) * f1YX(y,x) * log(f1YX(y,x))
  fToInt2 = function(x,y,f1X=fX,f1YX=fYX,g1=g,coef=p)    f1X(x) * f1YX(y,x) * log(g1(y,coef))

  KL1 = integral2(fToInt1,xMin,xMax,yMin,yMax)$Q;
  KL2 = integral2(fToInt2,xMin,xMax,yMin,yMax)$Q;
  
  #cat('\n intervalos x=',xMin,xMax,'y=',yMin,yMax,'\n')
  #cat('\n kL = ',KL1-KL2)
  
  return(KL1-KL2)
}
#########################################################





DivKLpolyXYZ = function(fY,fZ,fXYZ,p,xMin,xMax,yMin,yMax,zMin,zMax)
{
  # it calculates the conditional divergence Kullback-Leibler from a polynomical fitting
  # of a variable X with parents Y and Z using the densities and numerical integration
  # INPUT: - fY, fZ = original marginal density function of the parents   f(y),  f(z)
  #        - fXYZ = original conditional density function   f(x|y,z)
  #        - p = coefficients of the fitted polynomial
  #        - (xMin,xMax) = vector with the definition interval of the child X
  #        - (yMin,yMax) = vector with the definition interval of the parent Y
  #        - (zMin,zMax) = vector with the definition interval of the parent Z
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
  
  # function of the first integral
  fToInt1 = function(x,y,z,f1Y=fY,f1Z=fZ,f1XYZ=fXYZ,g1=g,coef=p)  
    f1Y(y) * f1Z(z) * f1XYZ(x,y,z) * log(f1XYZ(x,y,z));
  
  # function of the first integral
  fToInt2 = function(x,y,z,f1Y=fY,f1Z=fZ,f1XYZ=fXYZ,g1=g,coef=p)  
    f1Y(y) * f1Z(z) * f1XYZ(x,y,z) * log(g1(x,coef));
  
  #cat('\n \n intervalo a integrar en X =',c(xMin,xMax))
  #cat('\n intervalo a integrar en Y =',c(yMin,yMax))
  #cat('\n intervalo a integrar en Z =',c(zMin,zMax))
  
  KL1 = integral3(fToInt1,xMin,xMax,yMin,yMax,zMin,zMax);
  KL2 = integral3(fToInt2,xMin,xMax,yMin,yMax,zMin,zMax);
  
  #cat('\n integral-1 = ',KL1, '  integral-2 = ',KL2)
  
  return(KL1-KL2)
}
#########################################################






DivKLtmopXY = function(fX,fYX,tMoP,yMin,yMax)
{
  # it calculates the divergence Kullback-Leibler of a tMoP density
  # INPUT: - fY = original marginal density function of the parent   f(x)
  #        - fYX = original conditional density function   f(y|x)
  #        - tMoP = tMoP fit with polynomial and tails in X
  #        - (yMin,yMax) = vector with the definition interval of the child Y
  # OUTPUT: - KL = divergence KL
  
  
  # cut points of the discretization of the parent X
  discretX = tMoP[[3]];
  
  # initialize the divergence
  KL = 0;
  
  # for each interval
  for (i in 1:nrow(tMoP[[2]]))
  {
    # extract the information of the tmop in each part of the discretization
    values = tMoP[[1]][,(3*i-2):(3*i)];
    interval = values[,1];        tails = values[,3];
    
    # extract the coefficients that are not NA (degree lower than max)
    poly = tMoP[[2]][i,];
    #poly = p[!is.na(p)];
    
    # calculate the divergence when there are tails
    if (tails[1]>0)     
      KL = KL + DivKLpolyXY(fX,fYX,tails[1],discretX[i],discretX[i+1],yMin,interval[1])
    
    if (tails[2]>0)
      KL = KL + DivKLpolyXY(fX,fYX,tails[2],discretX[i],discretX[i+1],interval[2],yMax)
    
    # divergence of the polynomial
    KL = KL + DivKLpolyXY(fX,fYX,poly,discretX[i],discretX[i+1],interval[1],interval[2])
  }
  
  return(KL)
}
#########################################################







DivKLDiscretXY = function(fX,fYX,prob,breaks)
{
  # it calculates the divergence Kullback-Leibler of a combination of uniform densities
  # INPUT: - fX = original marginal density function of the parent   f(x)
  #        - fYX = original conditional density function   f(y|x)
  #        - prob = matrix of probabilities
  #        - breaks = matrix with 2 rows including the breaks of each discretization
  # OUTPUT: - KL = divergence KL
  
  
  # initialize the divergence
  KL = 0;
  
  for (i in 1:nrow(prob))
  {
    for (j in 1:ncol(prob))
    {
      #cat('\n \n intervalo en padre',breaks[2,j],',',breaks[2,j+1])
      #cat('\n intervalo en hijo',breaks[1,i],',',breaks[1,i+1])
      KL = KL+DivKLpolyXY(fX,fYX,prob[i,j],breaks[2,j],breaks[2,j+1],breaks[1,i],breaks[1,i+1])
      
      #cat('\n kl = ',KL)
    }
  }
  
  return(KL)
}
#########################################################






DivKL1parent = function(fX,fYX,tMoP,yMin,yMax)
{
  # it calculates the divergence Kullback-Leibler of a tMoP density
  # INPUT: - fY = original marginal density function of the parent   f(x)
  #        - fYX = original conditional density function   f(y|x)
  #        - tMoP = tMoP fit with polynomial and tails in X
  #        - (yMin,yMax) = vector with the definition interval of the child Y
  # OUTPUT: - KL = divergence KL
  
  
  # cut points of the discretization of the parent X
  discretX = tMoP[[3]];
  
  # initialize the divergence
  KL = 0;
  
  # for each interval
  for (i in 1:nrow(tMoP[[2]]))
  {
    # extract the information of the tmop in each part of the discretization
    values = tMoP[[1]][,(3*i-2):(3*i)];
    interval = values[,1];        tails = values[,3];
    
    # extract the coefficients that are not NA (degree lower than max)
    poly = tMoP[[2]][i,];
    #poly = p[!is.na(p)];
    
    # calculate the divergence when there are tails
    if (tails[1]>0)     
      KL = KL + DivKLpolyXY(fX,fYX,tails[1],discretX[i],discretX[i+1],yMin,interval[1])
    
    if (tails[2]>0)
      KL = KL + DivKLpolyXY(fX,fYX,tails[2],discretX[i],discretX[i+1],interval[2],yMax)
    
    # divergence of the polynomial
    KL = KL + DivKLpolyXY(fX,fYX,poly,discretX[i],discretX[i+1],interval[1],interval[2])
  }
  
  return(KL)
}
#########################################################






# ponemos el directorio de este script y lo cambiamos para cargar AuxiliarFunctions.R
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../../../core")
source("AuxiliarFunctions.R");
# ponemos el directorio de este escript y lo cambiamos para cargar los datos
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


DivKL2parent = function(fY,fZ,fXYZ,tMoP,xMin,xMax,yMin,yMax,zMin,zMax)
{
  # it calculates the divergence Kullback-Leibler of a tMoP density
  # INPUT: - fY, fZ = original marginal density function of the parents   f(y),  f(z)
  #        - fXYZ = original conditional density function   f(x|y,z)
  #        - tMoP = tMoP fit with polynomial and tails in X
  #        - (xMin,xMax) = vector with the definition interval of the child X
  #        - (yMin,yMax) = vector with the definition interval of the parent Y
  #        - (zMin,zMax) = vector with the definition interval of the parent Z
  # OUTPUT: - KL = divergence KL
  
  
  # cut points of the discretization of the parent Y
  discretYZ = tMoP[[3]];
  
  # initialize the divergence
  KL = 0;
  
  # the number of different values of discr is the size of the discretization
  indices = as.numeric(levels(factor(discretYZ)));
  nIntervals = length(indices);
  
  for (i in 1:nIntervals)
  {
    # calculate the corresponding interval in each parent
    cutPoints = CutPointsDiscr(discretYZ,yMin,yMax,zMin,zMax,i);
    discretY = cutPoints[1,];          discretZ = cutPoints[2,];
    
    # extract the information of the tmop in each part of the discretization
    values = tMoP[[1]][,(3*i-2):(3*i)];
    interval = values[,1];        tails = values[,3];
    
    # extract the coefficients that are not NA (degree lower than max)
    p = tMoP[[2]][i,];
    poly = p[!is.na(p)];
    
    # calculate the divergence when there are tails
    if (tails[1]>0)     
      KL = KL + DivKLpolyXYZ(fY,fZ,fXYZ,tails[1],xMin,interval[1],discretY[1],discretY[2],
                             discretZ[1],discretZ[2]);
    
    if (tails[2]>0)
      KL = KL + DivKLpolyXYZ(fY,fZ,fXYZ,tails[2],interval[2],xMax,discretY[1],discretY[2],
                             discretZ[1],discretZ[2]);
    
    # divergence of the polynomial
    KL = KL + DivKLpolyXYZ(fY,fZ,fXYZ,poly,interval[1],interval[2],discretY[1],discretY[2],
                           discretZ[1],discretZ[2]);
  }
  
  return(KL)
}
#########################################################





DivKLDiscret = function(fX,fYX,prob,breaks)
{
  # it calculates the divergence Kullback-Leibler of a combination of uniform densities
  # INPUT: - fX = original marginal density function of the parent   f(x)
  #        - fYX = original conditional density function   f(y|x)
  #        - prob = matrix of probabilities
  #        - breaks = matrix with 2 rows including the breaks of each discretization
  # OUTPUT: - KL = divergence KL
  
  
  # initialize the divergence
  KL = 0;
  
  for (i in 1:nrow(prob))
  {
    for (j in 1:ncol(prob))
    {
      #cat('\n \n intervalo en padre',breaks[2,j],',',breaks[2,j+1])
      #cat('\n intervalo en hijo',breaks[1,i],',',breaks[1,i+1])
    KL = KL+DivKLpolyXY(fX,fYX,prob[i,j],breaks[2,j],breaks[2,j+1],breaks[1,i],breaks[1,i+1])
    
    #cat('\n kl = ',KL)
    }
  }
  
  return(KL)
}
#########################################################





