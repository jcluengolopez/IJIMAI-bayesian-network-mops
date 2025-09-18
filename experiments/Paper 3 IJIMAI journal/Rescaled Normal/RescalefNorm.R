
# CARGAR DEPENDENCIAS AL PRINCIPIO (polinomios)
library(polynom)


RescalefNorm = function(tmop01,mu,sigma,xMin01,xMax01,xMin,xMax)
{
  # it takes a tmop of a N(0,1) and re-scale it to a N(mu,sigma) in the specified interval
  # which should be at least the original times sigma
  # f(x) = g((x-mu)/sigma) = g(x/sigma - mu/sigma)
  # INPUT: - tmop01 = tmop of the N(0,1)
  #        - (mu,sigma) = median and sd of the new normal distribution
  #        - (xMin01,xMax01) = minimum and maximum of the tmop in N(0,1)
  #        - (xMin,xMax) = minimum and maximum of the tmop in N(mu,sigma)
  # OUTPUT: - tmop = tmop of the N(mu,sigma)
  
  
  # roots, polynom and tails
  roots01 = tmop01[[1]];        poly01 = tmop01[[2]];        tails01 = tmop01[[3]];
  
  # new polynomial
  poly = poly01 / sigma
  
  # change the deviation
  for (i in 1:length(poly))     poly[i] = poly[i] / sigma^(i-1)
  
  # rescale the polynom to mu/sigma
  poly = coef(change.origin(as.polynomial(poly),-mu))
  
  # change the roots multiplying by sigma and adding mu
  roots = roots01 * sigma + mu
  
  # rescale the tails using the length of its interval
  tails = c(0,0)
  if (tails01[1]>0)       tails[1] = tails01[1] * (roots01[1] - xMin01) / (roots[1] - xMin)
  if (tails01[2]>0)       tails[2] = tails01[2] * (xMax01 - roots01[2]) / (xMax - roots[2])
  
  return(list(roots,poly,tails))
}
#########################################################










