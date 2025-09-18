
# CARGAR DEPENDENCIAS AL PRINCIPIO (polinomios y kernel)
library(polynom)
library(KernSmooth)


PolynomialFit = function(data,maxDegree=10,dataMin,dataMax,domain,degree,
                          verbose=FALSE,plots=FALSE)
{
  # it makes the polynomial fit calculating the kernel points
  # INPUT: - data = original data
  #        - maxDegree = maximum degree for the polynomial fit
  #        - dataMin, dataMax = the minimum and maximum value in the variable
  #        - degree = 1 to fit from degree 1 to max and 2 to fit only the max and reduce
  #        - domain = 1 if the fit is made in the complete domain or 2 if the fit is made
  #        in the restricted domain and the queues extended
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - roots = minimum and maximum values where the polynomial is evaluated
  #         - polyNorm = coefficients of the normalized polynomial
  #         - probTails = probability of each tail of the fit


  # check if there are 0 or 1 point
  if (length(table(data))<2)
  {
    # the fit is a constant with no roots except the min and max
    polyNorm = 1/(dataMax-dataMin);
    roots = c(dataMin,dataMax);

    fit = list(roots,polyNorm,c(0,0));
  }
  else
  {
    # calculate the kernel points (x,y)
    p = KernelPoints(data);
    x = p[,1];      y = p[,2];

    # when domain is 1, extend the kernel points to the whole domain
    if (domain==1)
    {
      # add the minimum and maximum data
      x = c(dataMin,x,dataMax);
      y = c(min(y),y,min(y));
    }

    # when degree is 1 fit from 1 to the maximum
    if (degree==1)
    {
      fit = Fit1(data,x,y,dataMin,dataMax,maxDegree,domain);
    }
    else         # fit maximum degree and reduce
    {
      fit = Fit2(data,x,y,dataMin,dataMax,maxDegree,domain);
    }
  }

  return(fit);
}
#########################################################








Fit1 = function(data,x,y,dataMin,dataMax,maxDegree,domain,verbose=FALSE,plots=FALSE)
{
  # it makes the polynomial fit using from 1 to the maximum degree and choose the best fit
  # Two kind of root's selection are used (longest interval and higher proportion of kernel
  # points)
  # INPUT: - data = original data
  #        - (x,y) = coordinates of kernel points
  #        - dataMin, dataMax = the minimum and maximum value in the variable (not kernel)
  #        - maxDegree = maximum degree for the polynomial fit
  #        - domain = 1 if the fit is made in the complete domain or 2 if the fit is made
  #        in the restricted domain and the queues extended
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - roots = selected roots
  #         - poly = coefficients of the normalized polynomial
  #         - probTails = probability of each tail of the fit


  # make the first fit
  degree = maxDegree;
  fit = lm(y~poly(x,maxDegree,raw=TRUE));

  # select the first term of the polynomial which is not NA
  degreeNA = which(is.na(fit$coefficients)==TRUE);

  # check if there are NAs
  if (length(degreeNA)>0)
  {
    # reduce the degree to remove terms with NA
    degree = degreeNA[1] - 2;
  }

  # initialize the lists for the fit
  listPoly  = vector('list',degree);
  listRoots = listPoly;              listProbTail = listPoly;
  listBIC = vector('numeric',degree);

  # from degree 1 to the maximum
  for (d in 1:degree)
  {
    # make the fit
    fit = lm(y~poly(x,d,raw=TRUE));

    # select the coefficients
    poly = fit$coefficients;

    # calculate the roots using two different methods
    roots1 = CalcRoots1(polynom :: as.polynomial(poly),dataMin,dataMax);
    roots2 = CalcRoots2(polynom :: as.polynomial(poly),x,y,dataMin,dataMax);

    # calculate the error in each case
    error1 = CalcError(polynom :: as.polynomial(poly),x,y,roots1);
    error2 = CalcError(polynom :: as.polynomial(poly),x,y,roots2);

    # select the best root option
    if (is.na(error2)==TRUE)   {   roots = roots1;   }
    else if (is.na(error1)==TRUE)   {   roots = roots2;   }
    else if (error1 < error2)   {   roots = roots1;   }
    else   {   roots = roots2;   }

    # obtain the probability of each queue
    tails = ProbTail(data,roots,dataMin,dataMax);
    probTail = tails[[1]];
    probTailTotal = tails[[2]];

    # normalize the polynomial
    poly = poly / polynom :: integral(polynom :: as.polynomial(poly),roots);
    poly = poly * (1 - probTailTotal);

    # calculate the BIC of the fit (including the tails)
    values = cbind(roots,c(dataMin,dataMax),probTail);
    #bic = BICContinuousVar(poly,data,values,verbose=FALSE);
    bic = BIC0Parent(poly,data,values,1);

    # fill in the lists
    listPoly[[d]]  = poly;       listRoots[[d]] = roots;
    listBIC[d] = bic;            listProbTail[[d]] = probTail;
  }

  # select the best model
  best = which(listBIC==max(listBIC));
  bestPoly = listPoly[[best]];
  bestRoots = listRoots[[best]];
  bestProbTail = listProbTail[[best]];

  return(list(bestRoots,bestPoly,bestProbTail))
}
#########################################################









Fit2 = function(data,x,y,dataMin,dataMax,maxDegree,domain,verbose=FALSE,plots=FALSE)
{
  # it makes the polynomial fit using only the maximum degree and reduce this degree
  # when the whole fit (with queues) reduces the error. Two kind of root's
  # selection are used (longest interval and higher proportion of kernel points)
  # INPUT: - data = original data
  #        - (x,y) = coordinates of kernel points
  #        - dataMin, dataMax = the minimum and maximum value in the variable (not kernel)
  #        - maxDegree = maximum degree for the polynomial fit
  #        - domain = 1 if the fit is made in the complete domain or 2 if the fit is made
  #        in the restricted domain and the queues extended
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - roots = selected roots
  #         - poly = coefficients of the normalized polynomial
  #         - probTails = probability of each tail of the fit


  # make the first fit
  degree = maxDegree;
  fit = lm(y~poly(x,maxDegree,raw=TRUE));

  # select the first term of the polynomial which is no t NA
  degreeNA = which(is.na(fit$coefficients)==TRUE);

  # check if there are NAs
  if (length(degreeNA)>0)
  {
    # reduce the degree to remove terms with NA
    degree = degreeNA[1] - 2;

    # fit using the new degree
    fit = lm(y~poly(x,degree,raw=TRUE));
  }

  # select the coefficients
  poly = fit$coefficients;

  # calculate the roots using two different methods
  roots1 = CalcRoots1(polynom :: as.polynomial(poly),dataMin,dataMax);
  roots2 = CalcRoots2(polynom :: as.polynomial(poly),x,y,dataMin,dataMax);

  # calculate the error in each case
  error1 = CalcError(polynom :: as.polynomial(poly),x,y,roots1);
  error2 = CalcError(polynom :: as.polynomial(poly),x,y,roots2);

  # select the best root option
  if (is.na(error2)==TRUE)   {   roots = roots1;   }
  else if (is.na(error1)==TRUE)   {   roots = roots2;   }
  else if (error1 < error2)   {   roots = roots1;   }
  else   {   roots = roots2;   }

  # obtain the probability of each tail
  tails = ProbTail(data,roots,dataMin,dataMax);
  probTail = tails[[1]];
  probTailTotal = tails[[2]];

  # normalize the polynomial
  poly = poly / polynom :: integral(polynom :: as.polynomial(poly),roots);
  poly = poly * (1 - probTailTotal);

  # calculate the BIC of the fit (including the tails)
  values = cbind(roots,c(dataMin,dataMax),probTail);
  #bic = BICContinuousVar(poly,data,values,verbose=FALSE);
  bic = BIC0Parent(poly,data,values,1);

  # initialize the bic, roots, tails and fit
  newBic = bic;             newRoots = roots;          newPoly = poly;
  newProbTail = probTail;   newProbTailTotal = probTailTotal;

  # repeat the fit while the bic increases
  while (newBic>=bic & degree>1)
  {
    # upload  bic, roots, tails and fit
    bic = newBic;             roots = newRoots;       poly = newPoly;
    probTail = newProbTail;   probTailTotal = newProbTailTotal;

    # there is one coefficient more than the degree
    degree = degree - 1;

    # make the fit
    fit = lm(y~poly(x,degree,raw=TRUE));

    # select the coefficients
    newPoly = fit$coefficients;

    # calculate the roots using two different methods
    roots1 = CalcRoots1(polynom :: as.polynomial(newPoly),dataMin,dataMax);
    roots2 = CalcRoots2(polynom :: as.polynomial(newPoly),x,y,dataMin,dataMax);

    # calculate the error in each case
    error1 = CalcError(polynom :: as.polynomial(newPoly),x,y,roots1);
    error2 = CalcError(polynom :: as.polynomial(newPoly),x,y,roots2);

    # select the best root option
    if (is.na(error2)==TRUE)   {   newRoots = roots1;   }
    else if (is.na(error1)==TRUE)   {   newRoots = roots2;   }
    else if (error1 < error2)   {   newRoots = roots1;   }
    else   {   newRoots = roots2;   }

    # obtain the probability of each tail
    tails = ProbTail(data,newRoots,dataMin,dataMax);
    newProbTail = tails[[1]];
    newProbTailTotal = tails[[2]];

    # normalize the polynomial
    newPoly = newPoly / polynom :: integral(polynom :: as.polynomial(newPoly),newRoots);
    newPoly = newPoly * (1 - newProbTailTotal);

    # calculate the BIC of the fit (including the tails)
    values = cbind(newRoots,c(dataMin,dataMax),newProbTail);
    #newBic = BICContinuousVar(newPoly,data,values,verbose=FALSE);
    newBic = BIC0Parent(newPoly,data,values,1);
  }

  return(list(roots,poly,probTail))
}
#########################################################









ProbTail = function(data,roots,dataMin,dataMax)
{
  # it calculates the probability of the tails between the roots and the extremes
  # according to the distribution of x
  # INPUT: - data = original data
  #        - roots = two roots that defines the interval where the polynomial is fitted
  #        - dataMin, dataMax = the minimum and maximum value in the variable (not kernel)
  # OUTPUT:  - probTail = probability of both tails
  #          - probTailTotal = probability of joint tails


  # initialize the joint probability of both tails and each one
  probTailTotal = 0;        probTailLeft = 0;        probTailRight = 0;

  # check if the domain has to be adapted in the left part
  if (roots[1]>dataMin)
  {
    # check if there are points in that area
    if (length(data[data<roots[1]])==0)          # prob = 0.1%
    {
      # assign a probability of the tail
      probTailLeft = 0.001;

      # add that to the total probability
      probTailTotal = probTailTotal + probTailLeft;

      # divide by the length of the interval
      probTailLeft = probTailLeft / abs(roots[1] - dataMin);
    }
    else
    {
      # calculate the probability of the left tail
      probTailLeft = length(data[data<roots[1]]) / length(data);

      # add that to the total probability
      probTailTotal = probTailTotal + probTailLeft;

      # divide by the length of the interval
      probTailLeft = probTailLeft / abs(roots[1] - dataMin);
    }
  }

  # check if the domain has to be adapted in the right part
  if (roots[2]<dataMax)
  {
    # check if there are points in that area
    if (length(data[data>roots[2]])==0)          # prob = 0.1%
    {
      # assign a probability of the tail
      probTailRight = 0.001;

      # add that to the total probability
      probTailTotal = probTailTotal + probTailRight;

      # divide by the length of the interval
      probTailRight = probTailRight / abs(roots[2] - dataMax);
    }
    else
    {
      # calculate the probability of the left tail
      probTailRight = length(data[data>roots[2]]) / length(data);

      # add that to the total probability
      probTailTotal = probTailTotal + probTailRight;

      # divide by the length of the interval
      probTailRight = probTailRight / abs(roots[2] - dataMax);
    }
  }

  # probability of both tails
  tails = c(probTailLeft,probTailRight);

  return(list(tails,probTailTotal))
}
#########################################################








KernelPoints = function(data,nPoints=3000,nIntervals=10)
{
  # it calculates the kernel points and divide the interval in 10 subintervals to
  # dismiss points according to the number of points in each of them
  # INPUT: - data = original data
  #        - nPoints = number of kernel points
  #        - nIntervals = number of intervals to split and reduce the number of points
  # OUTPUT: - sol = (x,y) coordinates of kernel points


  # remove repeated data in the dataset
  d = as.numeric(as.data.frame(table(data),stringsAsFactors=F)[,1]);

  # calculate the bandwidth
  bandWidth = dpik(d);

  # duplicate the points in the left side
  dataLeft = min(data) - (data - min(data));

  # duplicate the points in the right side
  dataRight = max(data) + (max(data) - data);

  # extended set of points adding the duplication
  dataExpand = c(data,dataLeft,dataRight);

  # get the kernel points in the extended dataset
  kernelExpand = KernSmooth :: bkde(dataExpand,gridsize=nPoints*4,bandwidth=bandWidth);

  # select the kernel points inside the original interval
  select = which(kernelExpand$x>=min(data) & kernelExpand$x<=max(data));
  kernelX = kernelExpand$x[select];           kernelY = kernelExpand$y[select];

  # normalize the points
  kernelY = kernelY / sum(kernelY) * sum(kernelExpand$y);

  # points to separate the intervals
  cutPoints = seq(min(data),max(data),length.out=(nIntervals+1));

  # initialize x and y
  x = c();   y = c();

  # repeat the process for each interval
  for (i in 1:nIntervals)
  {
    # select the points of X that is inside the interval
    insideX = data[data>=cutPoints[i] & data<cutPoints[i+1]];
    insideX = c(cutPoints[i],insideX,cutPoints[i+1]);

    # percentage of points of X in this interval
    prob = length(insideX) / (length(data) + 2*nIntervals);

    # number of points that would correspond to each interval according to the percentage
    nPointsInterval = ceiling(prob*nPoints);

    # select the kernel points in the interval
    select = (kernelX>=min(insideX) & kernelX<max(insideX));
    intervalX = kernelX[select];          intervalY = kernelY[select];

    # check if there should be more points than the number of points in the interval
    if (nPointsInterval>(length(intervalX)-2))   # select all the points
    {
      x = append(x,intervalX);     y = append(y,intervalY);
    }
    else   # select some points
    {
      # add two more points to the sequence
      s = seq(1,length(intervalX),length.out=nPointsInterval+2);

      # round to get integer numbers and remove the first and last number
      sRound = floor(s)[2:(length(s)-1)];

      # add the points to the vectors
      x = append(x,intervalX[sRound]);     y = append(y,intervalY[sRound]);
    }
  }

  # return a matrix with one column for x and another for y
  sol = cbind(x,y);

  return(sol);
}
#########################################################








CalcRoots1 = function(poly,minX,maxX)
{
  # it calculates the roots of the polynomial and select the two roots that form the
  # longest interval (including the min and max of data)
  # INPUT:  - poly = polynomial
  #         - minX, maxX = m?n and max value of the original data
  # OUTPUT:    - range = the two most spread roots


  # calculate the roots of the polynomial and select only the real ones
  zeros = solve(poly);        # funci?n gen?rica, no necesita el paquete
  zerosReal = Re(zeros[Im(zeros)==0]);

  # select the zeros inside the range of data and add the min and max
  if (length(zerosReal)==0)   {   zerosRange = c(minX,maxX);   }   # no zeros
  else   {   zerosRange = c(minX,zerosReal[zerosReal>minX & zerosReal<maxX],maxX);   }

  # calculate the length of the intervals
  diff = diff(zerosRange);

  # calculate the middle point of each interval and evaluate the polynomial in them
  midX = zerosRange[-length(zerosRange)] + diff(zerosRange)/2;
  midY = predict(poly,midX);

  # discard the intervals with negative values
  diff[midY<0] = -1;

  # return the longest interval
  range = zerosRange[c(which(diff==max(diff)),which(diff==max(diff))+1)];

  return(range);
}
#########################################################








CalcRoots2 = function(poly,x,y,minX,maxX)
{
  # it calculates the roots of the polynomial and select the two roots with more area
  # between them (highest and positive)
  # INPUT:  - poly = polynomial
  #         - x, y = points the polynomial is adjusted
  #         - minX, maxX = m?n and max value of the original data
  # OUTPUT:    - range = the two roots with highest area between them


  # calculate the roots of the polynomial and select only the real ones
  zeros = solve(poly);         # funci?n gen?rica, no necesita el paquete
  zerosReal = Re(zeros[Im(zeros)==0]);

  # create a vector with the real roots and the min and max
  if (length(zerosReal)==0)   {   zerosRange = c(minX,maxX);   }   # no zeros
  else   {   zerosRange = c(minX,zerosReal[zerosReal>minX & zerosReal<maxX],maxX);   }

  # create a factor to split x using the roots
  factor = cut(x,zerosRange);

  # add the probabilities of Y in each interval
  prob = as.vector(tapply(y,factor,sum));

  # when there is no data in any interval, the probability is 0
  prob[is.na(prob)] = 0;

  # calculate the middle point of each interval and evaluate the polynomial in them
  midX = zerosRange[-length(zerosRange)] + diff(zerosRange)/2;
  midY = predict(poly,midX);

  # discard the intervals with negative values
  prob[midY<0] = -1;

  # return the interval with highest probability
  range = zerosRange[c(which(prob==max(prob)),which(prob==max(prob))+1)];

  return(range);
}
#########################################################







CalcError = function(poly,x,y,roots)
{
  # it calculates the error of the polynomial fit It uses the polynomial between
  # the roots and the probability for the tails (not dividing by the range)
  # INPUT:  - poly = polynomial
  #         - x, y = points the polynomial is adjusted
  #         - roots = the two roots where the polynomial is defined
  # OUTPUT:    - error = mean square error of the polynomial and the queues


  # split the points (points in the interval, at the right and left)
  insideX = x[x>=roots[1] & x<=roots[2]];
  insideY = y[x>=roots[1] & x<=roots[2]];
  leftY = y[x<roots[1]];          rightY = y[x>roots[2]];

  # evaluate the polynomial only between the two roots
  yPoly = predict(poly,insideX);

  # mean square error between the two roots
  errorInside = sum((insideY - yPoly) ^2) / length(insideY);

  # error in the left queue
  if (length(leftY)==0)   {   errorLeft = 0;   }     # zero points
  else
  {
    # calculate the height of the queue
    probLeft = sum(leftY) / length(leftY);

    # MSE of the queue
    errorLeft = sum((leftY - probLeft) ^2) / length(leftY);
  }

  # error in the right queue
  if (length(rightY)==0)   {   errorRight = 0;   }     # zero points
  else
  {
    # calculate the height of the queue
    probRight = sum(rightY) / length(rightY);

    # MSE of the queue
    errorRight = sum((rightY - probRight)^2) / length(rightY);
  }

  # add the error of each part
  error = errorInside + errorLeft + errorRight;

  return(error);
}
#########################################################
















