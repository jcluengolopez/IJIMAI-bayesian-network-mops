
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





DivKLpolyXYZ = function(fX,fY,fZXY,p,xMin,xMax,yMin,yMax,zMin,zMax)
{
  # it calculates the conditional divergence Kullback-Leibler from a polynomical fitting
  # of a variable X with parents Y and Z using the densities and numerical integration
  # INPUT: - fY, fZ = original marginal density function of the parents   f(y),  f(z)
  #        - fXYZ = original conditional density function   f(x|y,z)
  #        - p = coefficients of the fitted polynomial
  #        - (xMin,xMax) = vector with the definition interval of the parent X
  #        - (yMin,yMax) = vector with the definition interval of the parent Y
  #        - (zMin,zMax) = vector with the definition interval of the child Z
  # OUTPUT: - KL = divergence KL
  
  
  # delete the NAs
  p = p[!is.na(p)];
  
  # convert the coefficients of the polynomial in a function to integrate
  g = function(z,coef=p)
  {
    s = 0;
    for (i in 1:length(coef))      s = s + coef[i]*z^(i-1);
    return(s)
  }
  
  # function of the first integral
  fToInt1 = function(x,y,z,f1X=fX,f1Y=fY,f1ZXY=fZXY,g1=g,coef=p)  
    f1X(x) * f1Y(y) * f1ZXY(z,x,y) * log(f1ZXY(z,x,y));
  
  # function of the first integral
  fToInt2 = function(x,y,z,f1X=fX,f1Y=fY,f1ZXY=fZXY,g1=g,coef=p)  
    f1X(x) * f1Y(y) * f1ZXY(z,x,y) * log(g1(z,coef));
  
  #cat('\n \n intervalo a integrar en X =',c(xMin,xMax))
  #cat('\n intervalo a integrar en Y =',c(yMin,yMax))
  #cat('\n intervalo a integrar en Z =',c(zMin,zMax))
  
  KL1 = 0;
  KL2 = 0;
  KL1 = try({integral3(fToInt1,xMin,xMax,yMin,yMax,zMin,zMax)}, silent=TRUE)
  KL2 = try({integral3(fToInt2,xMin,xMax,yMin,yMax,zMin,zMax)}, silent=TRUE)
  
  #cat('\n KL1 =',KL1,'     KL2 =',KL2)
  #cat('\n tipo de KL',typeof(KL2))
  if (is.character(KL1) | is.character(KL2))   return(0)
  else   return(KL1-KL2)
  #KL1 = integral3(fToInt1,xMin,xMax,yMin,yMax,zMin,zMax);
  #KL2 = integral3(fToInt2,xMin,xMax,yMin,yMax,zMin,zMax);
  
  #return(KL1-KL2)
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






# ponemos el directorio de este script y lo cambiamos para cargar AuxiliarFunctions.R
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../../../core")
source("AuxiliarFunctions.R");
# ponemos el directorio de este escript y lo cambiamos para cargar los datos
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


DivKLtmopXYZ = function(fX,fY,fZXY,tMoP,xMin,xMax,yMin,yMax,zMin,zMax)
{
  # it calculates the divergence Kullback-Leibler of a tMoP density
  # INPUT: - fY, fZ = original marginal density function of the parents   f(y),  f(z)
  #        - fXYZ = original conditional density function   f(x|y,z)
  #        - tMoP = tMoP fit with polynomial and tails in X
  #        - (xMin,xMax) = vector with the definition interval of the parent X
  #        - (yMin,yMax) = vector with the definition interval of the parent Y
  #        - (zMin,zMax) = vector with the definition interval of the child Z
  # OUTPUT: - KL = divergence KL
  
  
  # cut points of the discretization of the parents X and Y
  discretXY = tMoP[[3]];
  
  # initialize the divergence
  KL = 0;
  
  # the number of different values of discr is the size of the discretization
  indices = as.numeric(levels(factor(discretXY)));
  nIntervals = length(indices);
  
  for (i in 1:nIntervals)
  {
    # calculate the corresponding interval in each parent
    cutPoints = CutPointsDiscr(discretXY,xMin,xMax,yMin,yMax,i);
    discretX = cutPoints[1,];          discretY = cutPoints[2,];
    
    # extract the information of the tmop in each part of the discretization
    values = tMoP[[1]][,(3*i-2):(3*i)];
    interval = values[,1];        tails = values[,3];
    
    # extract the coefficients that are not NA (degree lower than max)
    p = tMoP[[2]][i,];
    poly = p[!is.na(p)];
    
    #cat('\n \n intervalo en X =',discretX[1],',',discretX[2])
    #cat('\n intervalo en Y =',discretY[1],',',discretY[2])
    
    
    # calculate the divergence when there are tails
    if (tails[1]>0)     
      KL = KL + DivKLpolyXYZ(fX,fY,fZXY,tails[1],discretX[1],discretX[2],discretY[1],
                             discretY[2],zMin,interval[1]);
    
    if (tails[2]>0)
      KL = KL + DivKLpolyXYZ(fX,fY,fZXY,tails[2],discretX[1],discretX[2],discretY[1],
                             discretY[2],interval[2],zMax);
    
    # divergence of the polynomial
    KL = KL + DivKLpolyXYZ(fX,fY,fZXY,poly,discretX[1],discretX[2],discretY[1],
                           discretY[2],interval[1],interval[2]);
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








DivKLDiscretXYZ = function(fX,fY,fZXY,prob,breaks)
{
  # it calculates the divergence Kullback-Leibler of a combination of uniform densities
  # INPUT: - fX = original marginal density function of the parent   f(x)
  #        - fY = original marginal density function of the parent   f(y)
  #        - fZXY = original conditional density function   f(z|x,y)
  #        - prob = 3d matrix of probabilities
  #        - breaks = matrix with 3 rows including the breaks of each discretization (Z,X,Y)
  # OUTPUT: - KL = divergence KL
  
  
  # initialize the divergence
  KL = 0;
  
  # dimensions of the array of probabilities
  n = nrow(prob)
  
  for (i in 1:n)   # hijo Z
  {
    for (j in 1:n)   # padre X
    {
      for (k in 1:n)   # padre Y
      {
        KL = KL + DivKLpolyXYZ(fX,fY,fZXY,prob[i,j,k],breaks[2,j],breaks[2,j+1],
                               breaks[3,k],breaks[3,k+1],breaks[1,i],breaks[1,i+1])
      }
    }
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


DivKLtmopXYZ2 = function(fX,fYX,fZXY,tMoP,xMin,xMax,yMin,yMax,zMin,zMax)
{
  # it calculates the divergence Kullback-Leibler of a tMoP density
  # INPUT: - fY, fZ = original marginal density function of the parents   f(y),  f(z)
  #        - fXYZ = original conditional density function   f(x|y,z)
  #        - tMoP = tMoP fit with polynomial and tails in X
  #        - (xMin,xMax) = vector with the definition interval of the parent X
  #        - (yMin,yMax) = vector with the definition interval of the parent Y
  #        - (zMin,zMax) = vector with the definition interval of the child Z
  # OUTPUT: - KL = divergence KL
  
  
  # cut points of the discretization of the parents X and Y
  discretXY = tMoP[[3]];
  
  # initialize the divergence
  KL = 0;
  
  # the number of different values of discr is the size of the discretization
  indices = as.numeric(levels(factor(discretXY)));
  nIntervals = length(indices);
  
  for (i in 1:nIntervals)
  {
    # calculate the corresponding interval in each parent
    cutPoints = CutPointsDiscr(discretXY,xMin,xMax,yMin,yMax,i);
    discretX = cutPoints[1,];          discretY = cutPoints[2,];
    
    # extract the information of the tmop in each part of the discretization
    values = tMoP[[1]][,(3*i-2):(3*i)];
    interval = values[,1];        tails = values[,3];
    
    # extract the coefficients that are not NA (degree lower than max)
    p = tMoP[[2]][i,];
    poly = p[!is.na(p)];
    
    #cat('\n \n intervalo en X =',discretX[1],',',discretX[2])
    #cat('\n intervalo en Y =',discretY[1],',',discretY[2])
    
    
    # calculate the divergence when there are tails
    if (tails[1]>0)     
      KL = KL + DivKLpolyXYZ2(fX,fYX,fZXY,tails[1],discretX[1],discretX[2],discretY[1],
                             discretY[2],zMin,interval[1]);
    
    if (tails[2]>0)
      KL = KL + DivKLpolyXYZ2(fX,fYX,fZXY,tails[2],discretX[1],discretX[2],discretY[1],
                             discretY[2],interval[2],zMax);
    
    # divergence of the polynomial
    KL = KL + DivKLpolyXYZ2(fX,fYX,fZXY,poly,discretX[1],discretX[2],discretY[1],
                           discretY[2],interval[1],interval[2]);
  }
  
  return(KL)
}
#########################################################







DivKLDiscretXYZ2 = function(fX,fYX,fZXY,prob,breaks)
{
  # it calculates the divergence Kullback-Leibler of a combination of uniform densities
  # INPUT: - fX = original marginal density function of the parent   f(x)
  #        - fY = original marginal density function of the parent   f(y)
  #        - fZXY = original conditional density function   f(z|x,y)
  #        - prob = 3d matrix of probabilities
  #        - breaks = matrix with 3 rows including the breaks of each discretization (Z,X,Y)
  # OUTPUT: - KL = divergence KL
  
  
  # initialize the divergence
  KL = 0;
  
  # dimensions of the array of probabilities
  n = nrow(prob)
  
  for (i in 1:n)   # hijo Z
  {
    for (j in 1:n)   # padre X
    {
      for (k in 1:n)   # padre Y
      {
        KL = KL + DivKLpolyXYZ2(fX,fYX,fZXY,prob[i,j,k],breaks[2,j],breaks[2,j+1],
                                breaks[3,k],breaks[3,k+1],breaks[1,i],breaks[1,i+1])
      }
    }
  }
  
  
  return(KL)
}
#########################################################






DivKLpolyXYZ2 = function(fX,fYX,fZXY,p,xMin,xMax,yMin,yMax,zMin,zMax)
{
  # it calculates the conditional divergence Kullback-Leibler from a polynomical fitting
  # of a variable X with parents Y and Z using the densities and numerical integration
  # INPUT: - fY, fZ = original marginal density function of the parents   f(y),  f(z)
  #        - fXYZ = original conditional density function   f(x|y,z)
  #        - p = coefficients of the fitted polynomial
  #        - (xMin,xMax) = vector with the definition interval of the parent X
  #        - (yMin,yMax) = vector with the definition interval of the parent Y
  #        - (zMin,zMax) = vector with the definition interval of the child Z
  # OUTPUT: - KL = divergence KL
  
  
  # delete the NAs
  p = p[!is.na(p)];
  
  # convert the coefficients of the polynomial in a function to integrate
  g = function(z,coef=p)
  {
    s = 0;
    for (i in 1:length(coef))      s = s + coef[i]*z^(i-1);
    return(s)
  }
  
  # function of the first integral
  fToInt1 = function(x,y,z,f1X=fX,f1YX=fYX,f1ZXY=fZXY,g1=g,coef=p)  
    f1X(x) * f1YX(y,x) * f1ZXY(z,x,y) * log(f1ZXY(z,x,y));
  
  # function of the first integral
  fToInt2 = function(x,y,z,f1X=fX,f1YX=fYX,f1ZXY=fZXY,g1=g,coef=p)  
    f1X(x) * f1YX(y,x) * f1ZXY(z,x,y) * log(g1(z,coef));
  
  #cat('\n \n intervalo a integrar en X =',c(xMin,xMax))
  #cat('\n intervalo a integrar en Y =',c(yMin,yMax))
  #cat('\n intervalo a integrar en Z =',c(zMin,zMax))
  
  KL1 = 0;
  KL2 = 0;
  KL1 = try({integral3(fToInt1,xMin,xMax,yMin,yMax,zMin,zMax)}, silent=TRUE)
  KL2 = try({integral3(fToInt2,xMin,xMax,yMin,yMax,zMin,zMax)}, silent=TRUE)
  
  #cat('\n KL1 =',KL1,'     KL2 =',KL2)
  #cat('\n tipo de KL',typeof(KL2))
  if (is.character(KL1) | is.character(KL2))   return(0)
  else   return(KL1-KL2)
  #KL1 = integral3(fToInt1,xMin,xMax,yMin,yMax,zMin,zMax);
  #KL2 = integral3(fToInt2,xMin,xMax,yMin,yMax,zMin,zMax);
  
  #return(KL1-KL2)
}
#########################################################





