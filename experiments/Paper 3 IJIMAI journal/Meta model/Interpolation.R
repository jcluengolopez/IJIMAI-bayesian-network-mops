
# CARGAR DEPENDENCIAS AL PRINCIPIO (polinomios)
library(polynom)


Interpolate1param = function(param,tMoPs,xMin,xMax,d)
{
  # it makes the interpolation of a distribution with 1 parameter using a grid
  # INPUT: - (alpha,beta) = parameters of the distribution to estimate
  #        - tMoPs =list with all the tMoPs fitted in the grid
  #        - (xMin,xMax) = minimum and maximum value of the parameter
  #        - d = distance between points in the grid
  # OUTPUT: - newtMoP = tMoP estimated from the grid
  
  
  # calculate the equispaced points
  x = seq(xMin,xMax,d);
  
  # select the points of the grid (check if param is in the grid)
  if (length(x[x==param])>0)     nX = c(which(x==param),which(x==param))
  else     nX = c(max(which(x<param)),max(which(x<param))+1);
  
  # when the distribution is already done
  if (nX[1]==nX[2])     return(tMoPs[[nX[1]]])
  else
  {
    # calculate the distance from param to the points of the grid
    dist = c(param - x[nX[1]], x[nX[2]] - param);
    #distance1 = param - nX[1];
    #distance1 = nX[2] - param;
    #cat('\n distancias = ',dist, 'a los puntos ',nX);
    
    # mix the two closest tMoPs
    newtMoP = MixtMoP(tMoPs,nX,dist);
    
    return(newtMoP)
  }
}
#########################################################





Interpolate2param = function(alpha,beta,tMoPs,xMin,xMax,yMin,yMax,d)
{
  # it makes the interpolation of a distribution with 2 parameters using a grid
  # INPUT: - (alpha,beta) = parameters of the distribution to estimate
  #        - tMoPs =list with all the tMoPs fitted in the grid
  #        - (xMin,xMax) = minimum and maximum value of alpha
  #        - (xMin,xMax) = minimum and maximum value of beta
  #        - d = distance between points in the grid
  # OUTPUT: - newtMoP = tMoP estimated from the grid
  
  
  # vector with the parameters
  param = c(alpha,beta);
  
  # calculate the points of the grid
  x = seq(xMin,xMax,d);        y = seq(yMin,yMax,d);
  
  # select the indices the four points of the square where the new parameters are
  a = SelectGridPoints(xMin,xMax,yMin,yMax,x,y,alpha,beta);
  nX = a[1,];          nY = a[2,];
  #cat('\n nX = ',nX,'   nY = ',nY)
  
  # when the distribution is already done
  if (nX[1]==nX[2] & nY[1]==nY[2])
    return(tMoPs[[(nX[1]-1)*length(x) + nY[1]]])
  
  else
  {
    # four points and their indices in the grid
    i1 = c(nX[1],nY[1]);              i2 = c(nX[2],nY[1]);
    i3 = c(nX[1],nY[2]);              i4 = c(nX[2],nY[2]);
    p1 = c(x[nX[1]],y[nY[1]]);              p2 = c(x[nX[2]],y[nY[1]]);
    p3 = c(x[nX[1]],y[nY[2]]);              p4 = c(x[nX[2]],y[nY[2]]);
    
    # matrices with the indices and points
    ind = rbind(i1,i2,i3,i4);
    p = rbind(p1,p2,p3,p4);
    
    #cat('\n x = ',x);
    #cat('\n y = ',y);
    #cat('\n está entre los índices ',nX,' y ',nY);
    #cat('\n alpha = ',alpha,'     beta = ',beta);
    #cat('\n índices = \n');  print(ind);
    #cat('\n puntos = \n');  print(p);
    
    # calculate the distance from (alpha,beta) to the selected four points of the grid
    dist = c(sum((param-p1)^2), sum((param-p2)^2), sum((param-p3)^2), sum((param-p4)^2));
    
    # obtain the indices and the distances of the two nearest points
    points = c(which.min(dist),which(dist==sort(dist)[2]));
    distance1 = dist[points[1]];
    ind1 = ind[points[1],];
    distance2 = dist[points[2]];
    ind2 = ind[points[2],];
    
    #cat('\n distancias = ',dist);
    #cat('\n points = ',points)
    #cat('\n punto1 con índices ',ind1,' y distancia = ',distance1);
    #cat('\n punto2 con índices ',ind2,' y distancia = ',distance2);
    
    # position of the two selected points and the distance to each one
    n = c(length(x)*(ind1[1]-1) + ind1[2], length(x)*(ind2[1]-1) + ind2[2])
    d = c(distance1,distance2);
    
    # mix the two closest tMoPs
    newtMoP = MixtMoP(tMoPs,n,d);
    
    return(newtMoP)
  }
}
#########################################################






SelectGridPoints = function(xMin,xMax,yMin,yMax,x,y,alpha,beta)
{
  # it calculates the points of the grid between the parameters are
  # INPUT: - (xMin,xMax) = minimum and maximum value of alpha
  #        - (xMin,xMax) = minimum and maximum value of beta
  #        - (x,y) = vectos with the points of the grid in each dimension
  #        - (alpha,beta) = parameters of the distribution to estimate
  # OUTPUT: - nX = vector with the two indices of X between alpha is
  #         - nY = vector with the two indices of X between beta is
  
  
  # check if hte value is in the grid of X
  if (length(x[x==alpha])>0)     nX = c(which(x==alpha),which(x==alpha))
  else     nX = c(max(which(x<alpha)),max(which(x<alpha))+1);
  
  # check if hte value is in the grid of Y
  if (length(y[y==beta])>0)     nY = c(which(y==beta),which(y==beta))
  else     nY = c(max(which(y<beta)),max(which(y<beta))+1);
  
  return(rbind(nX,nY))
}
#########################################################





MixtMoP = function(tMoPs,n,dist)
{
  # it creates a new tMoP mixing the two closest tMoPs of the grid
  # INPUT: - tMoPs =list with all the tMoPs fitted in the grid
  #        - n = positions of the two closest tMoPs of the grid
  #        - dist = dinstance to the two closest tMoPs of the grid
  # OUTPUT: - newtMoP = tMoP estimated from the grid
  
  
  # define the interval where the tmop is defined
  xMin = -20;    xMax = 20;
  
  # select the corresponding tMoPs
  model1 = tMoPs[[n[1]]];
  model2 = tMoPs[[n[2]]];
  roots1 = model1[[1]];        poly1 = model1[[2]];        probTails1 = model1[[3]];
  roots2 = model2[[1]];        poly2 = model2[[2]];        probTails2 = model2[[3]];
  
  # the weights for the mean are inversely proportional to the distances
  factor = 1 / (1/dist[1] + 1/dist[2]);
  w1 = factor / dist[1];          w2 = factor / dist[2];
  
  # new tails
  probTailsNew = w1*probTails1 + w2*probTails2;
  
  # new roots (the roots where both polynomials are positive)
  rootsNew = c(max(roots1[1],roots2[1]), min(roots1[2],roots2[2]));
  
  #cat('\n \n \n INFORMACIÓN DE LOS MODELOS ORIGINALES');
  #cat('\n los pesos son ',w1,' y ',w2);
  #cat('\n raíces1 = ',roots1,' y raíces2 = ',roots2);
  #cat('\n prob colas de modelo 1 = ',probTails1);
  #cat('\n prob colas de modelo 2 = ',probTails2);
  #cat('\n polinomio 1 = ',poly1);
  #cat('\n polinomio 2 = ',poly2);
  
  # new polynomial
  polyNew = coef( w1 * as.polynomial(poly1) + w2 * as.polynomial(poly2));
  
  #cat('\n \n \n INFORMACIÓN DEL NUEVO MODELO');
  #cat('\n raíces = ',rootsNew);
  #cat('\n prob colas = ',probTailsNew);
  #cat('\n polinomio = ',polyNew);
  
  # calculate the integral of the polynomial
  polyInt = polynom::integral(as.polynomial(polyNew),rootsNew);
  
  # calculate the integral of the tails
  tailsInt = c(0,0)
  if (probTailsNew[1]>0)
    tailsInt[1] = polynom::integral(as.polynomial(probTailsNew[1]),c(xMin,rootsNew[1]));
  if (probTailsNew[2]>0)
    tailsInt[2] = polynom::integral(as.polynomial(probTailsNew[2]),c(rootsNew[2],xMax));
  
  # integral in the total domain
  totalInt = polyInt + sum(tailsInt);
  
  #cat('\n \n \n INTEGRALES')
  #cat('\n integral polinomio = ',polyInt);
  #cat('\n integral colas = ',tailsInt);
  #cat('\n integral completa = ',totalInt);
  
  # normalize
  polyNorm = polyNew / totalInt;
  
  # normalize the tails when they exist
  probTailsNorm = c(0,0);
  if (probTailsNew[1]>0)   probTailsNorm[1] = probTailsNew[1] / totalInt;
  if (probTailsNew[2]>0)   probTailsNorm[2] = probTailsNew[2] / totalInt;
  
  return(list(rootsNew,polyNorm,probTailsNorm));
}
#########################################################



