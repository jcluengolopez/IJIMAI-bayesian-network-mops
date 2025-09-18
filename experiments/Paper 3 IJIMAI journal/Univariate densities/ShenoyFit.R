
# CARGAR DEPENDENCIAS AL PRINCIPIO (polinomios y taylor)
library(polynom)
library(rootSolve)



ShenoyFit = function(f,maxDegree,xMin,xMax,nIntervals)
{
  # it makes the polynomial fit using Taylor series in each of the nIntervals. It uses
  # the maximum degree and decreases when there are roots
  # INPUT: - f = original density function
  #        - maxDegree = maximum degree for the polynomial fit
  #        - (xMin,xMax) = vector with the interval where it is fitted
  #        - nIntervals = number of intervals
  # OUTPUT: - polynomial = matrix with the coefficients of the normalized polynomials
  
  
  # a vector for the split points
  splitPoints = seq(xMin,xMax,length.out=nIntervals+1);
  
  # create the list for the taylor functions
  listTaylor = vector(mode='list', length=nIntervals);
  
  # for each interval
  for (i in 1:nIntervals)
  {
    # set the initial degree
    d = maxDegree;
    
    # the centre of the interval is the centre of the Taylor series
    centre = mean(splitPoints[i:(i+1)]);
    
    # create the Taylor series
    taylor = taylor2(f,n=d,a=centre);
    
    # evaluate the function in a vector inside the interval
    a = seq(splitPoints[i],splitPoints[i+1],length.out=1000);
    y = taylor(a);
    
    # negative values
    neg = y[y<0];
    
    # repeat while there are roots inside that interval
    while (length(neg)>0 & d>1)
    {
      # repeat the process decreasing the degree
      d = d - 1;
      
      # create the Taylor series
      taylor = taylor2(f,n=d,a=centre);
      
      # evaluate the function in a vector inside the interval
      a = seq(splitPoints[i],splitPoints[i+1],length.out=10000);
      y = taylor(a);
      
      # negative values
      neg = y[y<0];
    }
    
    if (d==1 & length(neg)>0)
    {
      taylor = function(x)   0*x + 1    # se ajusta una uniforme
      d = 0;
    }
    
    # normalize so that it integers the same as the original function in that interval
    i1 = integrate(taylor,splitPoints[i],splitPoints[i+1])$value;
    i2 = integrate(f,splitPoints[i],splitPoints[i+1])$value;
    
    # multiply the function by i2 and divide by i1
    taylorNorm = function(x)  x
    body(taylorNorm) = parse(text=paste(i2,'*(',as.expression(body(taylor)),')/',i1));
    
    listTaylor[[i]] = taylorNorm;
  }
  
  return(listTaylor);
}
#########################################################








taylor2 <- function(f, n, a) {
  
  ith_derivative <- as.expression(body(f))
  f_temp <- function(x) x
  series <- as.character(f(a))
  for (i in seq_len(n)) {
    ith_derivative <- body(f_temp) <- D(ith_derivative, "x")
    series <- paste0(series, "+", f_temp(a) / factorial(i), "*(x - ", a, ")^", i)
  }
  f_output <- function(x) x
  body(f_output) <- parse(text = series)
  f_output
}






