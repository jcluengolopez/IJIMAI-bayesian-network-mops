
# CARGAR DEPENDENCIAS AL PRINCIPIO (polinomios y kernel)
library(polynom)



BIC0Parent = function(prob,data,dataVal,varType,verbose=FALSE)
{
  # it calculates the BIC of a variable without parents (in continuous it counts coefficients
  # of the polynomial and one parameter per tail)
  # INPUT: - prob = the probability distribution of the variable
  #        - data = the observations of the variable
  #        - dataVal = all the values of the variable in the data frame
  #        - varType = 0 if it is discrete and 1 if it is continuous
  #        - verbose = show the 'cats' on screen
  # OUTPUT: - bic = BIC of the distribution
  
  
  # different process for continuous and discrete variable
  if (varType==0 | varType==2)     # discrete or discretized
  {
    # calculate the BIC using the formula
    bic = sum(log10(prob)) - (length(prob) - 1) / 2 * log10(length(data));
  }
  else   # continuous
  {
    # initialize the vector for the logarithms
    logData = rep(0,length(data));
    
    # extract the coefficients of the polynomial
    poly = prob[is.na(prob)==FALSE];
    
    # change the name of the variable
    values = dataVal;
    
    # part of the data corresponding to the polynomial
    dataCentral = (data >= values[1,1] & data <= values[2,1]);
    
    # predict the polynomial for the exact value
    polyVal = predict(polynom :: as.polynomial(poly),data[dataCentral]);
    
    #cat('\n raices = ',values[,1]);
    #cat('\n ');
    #print(values);
    #cat('\n polyVal = ',length(polyVal));
    
    # each part of the fit
    logData[data < values[1,1]] = log10(values[1,3]);
    logData[data > values[2,1]] = log10(values[2,3]);
    logData[dataCentral] = log10(polyVal);
    
    # change NAs by 0 (logarithms produce then when the polynomial is 0)
    logData[is.na(logData)] = 0;
    
    # calculate the number of tails
    q = values[,3] > 0;
    nTails = length(q[q==TRUE]);
    
    # calculate the BIC using the formula
    nCoef = length(poly) - 1 + nTails;
    bic = sum(logData) - nCoef / 2 * log10(length(data));
  }
  
  return(bic);
}
############################################################







# BICContinuousVar = function(poly,data,values,verbose=FALSE)
# {
#   # it calculates the BIC of a continuous variable
#   # INPUT: - poly = coefficients of the polynomial fit
#   #        - data = the observations of the variable
#   #        - values = matrix with roots, extremes and tails
#   #        - verbose = show the 'cats' on screen
#   # OUTPUT: - bic = BIC of the distribution
#   
#   
#   # part of the data corresponding to the polynomial
#   dataCentral = (data >= values[1,1] & data <= values[2,1]);
#   
#   # predict the polynomial for the exact value
#   polyVal = predict(polynom :: as.polynomial(poly),data[dataCentral]);
#   
#   # initialize the vector for the logarithms
#   logData = rep(0,length(data));
#   
#   # each part of the fit
#   logData[data < values[1,1]] = log10(values[1,3]);
#   logData[data > values[2,1]] = log10(values[2,3]);
#   logData[dataCentral] = log10(polyVal);
#   
#   
#   #logLeft = (data <= values[1,1]) * log10(values[1,3]);
#   #logRight = (data >= values[2,1]) * log10(values[2,3]);
#   #logPoly = (data > values[1,1] & data < values[2,1]) * log10(polyVal);
#   
#   # change NAs by 0 (logarithms produce then when the polynomial is 0)
#   logData[is.na(logData)] = 0;
#   #logLeft[is.na(logLeft)] = 0;
#   #logRight[is.na(logRight)] = 0;
#   #logPoly[is.na(logPoly)] = 0;
#   
#   # sum of logarithms
#   #logData = logLeft + logRight + logPoly;
#   
#   # calculate the number of tails
#   q = values[,3] > 0;
#   nTails = length(q[q==TRUE]);
#   
#   # calculate the BIC using the formula
#   nCoef = length(poly) - 1 + nTails;
#   bic = sum(logData) - nCoef / 2 * log10(length(data));
#   
#   #cat('\n \n Verosimilitud = ',sum(logData));
#   #cat('\n \n Penalizacion parametros = ',nCoef / 2 * log10(length(data)));
#   
#   
#   return(bic);
# }
############################################################









BIC1ParentD = function(prob,data,dataVal,C,classVal,varType,verbose=FALSE)
{
  # it calculates the BIC of a variable with 1 discrete parent
  # INPUT: - prob = the probability distribution of the variable (rows for dataVal)
  #        - data = the observations of the variable
  #        - dataVal = all the values of the variable in the data frame
  #        - C = the observations of the parent variable
  #        - classVal = the different values of the parent
  #        - varType = 0 if it is discrete and 1 if it is continuous
  #        - verbose = show the 'cats' on screen
  # OUTPUT: - bic = BIC of the distribution
  
  
  # initialize the BIC
  bic = 0;
  
  # different process for continuous and discrete variable
  if (varType==0 | varType==2)     # discrete or discretized
  {
    # for each value of the class and data
    for (i in 1:length(classVal))
    {
      # calculate the BIC when there are data
      if (length(data[C==classVal[i]])>0)
      {
        bic = bic + BIC0Parent(prob[,i],data[C==classVal[i]],dataVal,varType);
      }
    }
  }
  else     # continuous
  {
    # for each value of the parent
    for (i in 1:length(classVal))
    {
      #cat('\n Entra');
      # extract the matrix of values and the polynomial
      poly = prob[i,is.na(prob[i,])==FALSE];
      values = dataVal[,(3*i-2):(3*i)];
      
      #cat('\n se usa la fila ',i,' de prob = \n ');
      #print(prob);
      
      
      
      
      
      # add the BIC of the new factor when there are data
      if (length(data[C==classVal[i]])>0)
      {
        #cat('\n \n numero de datos en D = ',length(data[C==classVal[i]]));
        
        bic = bic + BIC0Parent(poly,data[C==classVal[i]],values,varType);
      }
    }
  }
  
  # subtract the number of parameters in classVal
  #bic = bic - (length(classVal) - 1) / 2 * log10(length(data));
  
  return(bic);
}
############################################################







BIC1ParentC = function(prob,data,dataVal,C,cutPoints,varType,verbose=FALSE)
{
  # it calculates the BIC of a variable with 1 continuous parent
  # INPUT: - prob = the probability distribution of the variable (rows for dataVal)
  #        - data = the observations of the variable
  #        - dataVal = all the values of the variable in the data frame
  #        - C = the observations of the parent variable
  #        - cutPoints = a vector with the cut points in the discretization of C
  #        - varType = 0 if it is discrete and 1 if it is continuous
  #        - verbose = show the 'cats' on screen
  # OUTPUT: - bic = BIC of the distribution
  
  
  # extract the used cutPoints of the discretization
  cutPoints = cutPoints[is.na(cutPoints)==FALSE];
  
  # initialize the BIC
  bic = 0;
  
  # different process for continuous and discrete variable
  if (varType==0 | varType==2)     # discrete or discretized
  {
    # for each interval of the discretization of the class (except the last interval)
    for (i in 1:(length(cutPoints) - 2))
    {
      # select the data for this interval
      dataClass = (data>=cutPoints[i] & data<cutPoints[i+1]);
      
      # calculate the BIC when there are data
      if (length(data[dataClass])>0)
      {
        bic = bic + BIC0Parent(prob[,i],data[dataClass],dataVal,varType);
      }
    }
    
    # for the last interval (to include the extremes)
    i = length(cutPoints) - 1;
    
    # select the data for this interval
    dataClass = (data>=cutPoints[i] & data<=cutPoints[i+1]);
    
    # calculate the BIC when there are data
    if (length(data[dataClass])>0)
    {
      bic = bic + BIC0Parent(prob[,i],data[dataClass],dataVal,varType);
    }
  }
  else     # continuous
  {
    # for each interval of the discretization of the class (except the last interval)
    for (i in 1:(length(cutPoints) - 2))
    {
      # extract the matrix of values, the polynomial and the data
      #cat('\n Entra 2');
      poly = prob[i,is.na(prob[i,])==FALSE];
      values = dataVal[,(3*i-2):(3*i)];
      dataClass = data[data>=cutPoints[i] & data<cutPoints[i+1]];
      
      #cat('\n Sale 2');
      
      
      
      # add the BIC of the new factor when there are data
      if (length(dataClass)>0)
      {
        #cat('\n \n numero de datos en C = ',length(data[dataClass]));
        
        bic = bic + BIC0Parent(poly,dataClass,values,varType);
      }
    }
    
    # for the last interval (to include the extremes)
    i = length(cutPoints) - 1;
    
    # extract the matrix of values, the polynomial and the data
    poly = prob[i,is.na(prob[i,])==FALSE];
    values = dataVal[,(3*i-2):(3*i)];
    dataClass = data[data>=cutPoints[i] & data<=cutPoints[i+1]];
    
    # add the BIC of the new factor when there are data
    if (length(dataClass)>0)
    {
      bic = bic + BIC0Parent(poly,dataClass,values,varType);
    }
  }
  
  # subtract the number of parameters in classVal
  #bic = bic - (length(cutPoints) - 1) / 2 * log10(length(data));
  
  return(bic);
}
############################################################









BIC2ParentsDD = function(prob,data,dataVal,parent,parentVal,C,classVal,varType,verbose=FALSE)
{
  # it calculates the BIC of a variable with 2 discrete parents
  # INPUT: - prob = the probability distribution of the variable
  #        - data = the observations of the variable
  #        - dataVal = all the values of the variable in the data frame
  #        - C, parent = the observations of both parents
  #        - classVal, parentVal = the different values of both parents
  #        - varType = 0 if it is discrete and 1 if it is continuous
  #        - verbose = show the 'cats' on screen
  # OUTPUT: - bic = BIC of the distribution
  
  
  # initialize the BIC
  bic = 0;
  
  # for each value of the second parent
  for (i in 1:length(parentVal))
  {
    # selected data
    select = (parent==parentVal[i]);
    
    # calculate the BIC with 1 parent for this part of the model when there are data
    if (length(data[select])>0)
    {
      # different process for continuous and discrete variable
      if (varType==0 | varType==2)     # discrete or discretized
      {
        bic = bic + BIC1ParentD(prob[,i,],data[select],dataVal,C[select],classVal,varType);
      }
      else     # continuous
      {
        # transpose the matrix prob to have classVal in rows
        bic = bic + BIC1ParentD(t(prob[i,,]),data[select],dataVal[(2*i-1):(2*i),],
                                C[select],classVal,varType);
      }
      
      #bic = bic + BIC1ParentD(prob[,,i],data[C==classVal[i]],dataVal,parent[C==classVal[i]],
      #                  parentVal,varType,verbose);
    }
  }
  
  # subtract the number of parameters in classVal
  #bic = bic - (length(classVal) - 1) / 2 * log10(length(data));
  
  return(bic);
}
############################################################







BIC2ParentsCD = function(prob,data,dataVal,parent,cutPoints,C,classVal,varType,verbose=FALSE)
{
  # it calculates the BIC of a variable with 1 discrete and 1 continuous parent
  # INPUT: - prob = the probability distribution of the variable
  #        - data = the observations of the variable
  #        - dataVal = all the values of the variable in the data frame
  #        - parent = the observations of the continuous parent
  #        - cutPoints = a vector with the cut points in the discretization of parent
  #        - C = the observations of the discrete parent
  #        - classVal = the different values of the first parent
  #        - varType = 0 if it is discrete and 1 if it is continuous
  #        - verbose = show the 'cats' on screen
  # OUTPUT: - bic = BIC of the distribution
  
  
  # initialize the BIC
  bic = 0;
  
  # for each value of the discrete parent
  for (i in 1:length(classVal))
  {
    # calculate the BIC with 1 parent for this part of the model when there are data
    if (length(data[C==classVal[i]])>0)
    {
      # different process for continuous and discrete variable
      if (varType==0 | varType==2)     # discrete or discretized
      {
        bic = bic + BIC1ParentC(prob[,,i],data[C==classVal[i]],dataVal,
                                parent[C==classVal[i]],cutPoints[i,],varType,verbose);
      }
      else     # continuous
      {
        bic = bic + BIC1ParentC(prob[,,i],data[C==classVal[i]],dataVal[(2*i-1):(2*i),],
                                parent[C==classVal[i]],cutPoints[i,],varType,verbose);
      }
    }
  }
  
  # subtract the number of parameters in classVal
  #bic = bic - (length(classVal) - 1) / 2 * log10(length(data));
  
  
  
  return(bic);
}
############################################################




  





BIC2ParentsCC = function(prob,data,dataVal,C,parent,discr,varType,verbose=FALSE)
{
  # it calculates the BIC of a variable with 2 continuous parents
  # INPUT: - prob = the probability distribution of the variable  (columns for dataVal)
  #        - data = the observations of the variable
  #        - dataVal = all the values of the variable in the data frame
  #        - C, parent = the observations of both parents
  #        - discr = matrix whith the discretization of C and parent
  #        - varType = 0 if it is discrete and 1 if it is continuous
  #        - verbose = show the 'cats' on screen
  # OUTPUT: - bic = BIC of the distribution
  
  
  # the number of different values of discr is the size of the discretization
  indices = as.numeric(levels(factor(discr)));
  nIntervals = length(indices);
  
  # discretize the class and parent according to the matrix discr
  parentDisc = DiscretizeParent(discr,parent,C);
  
  # number of cut points in the discretization of the parents
  sizeDiscr = SizeDiscr(discr);
  
  # different process for continuous and discrete variable
  if (varType==0 | varType==2)     # discrete or discretized
  {
    # calculate the BIC using the discretization of both parents as a new variable
    # transpose the matrix prob to have dataVal in rows
    bic = BIC1ParentD(t(prob),data,dataVal,parentDisc,indices,varType);
  }
  
  else     # continuous
  {
    # calculate the BIC using the discretization of both parents as a new variable
    bic = BIC1ParentD(prob,data,dataVal,parentDisc,indices,varType);
  }
  
  # upload the BIC adding the rest of split points
  #bic = bic - (sum(sizeDiscr) - 1 - nIntervals) / 2 * log10(length(data));
  
  return(bic);
}
############################################################








BICNetwork = function(links,BICpartX,BICpartXY,BICpartXYZ)
{
  # it calculates the BIC of a network using the BIC of each distribution
  # INPUT: - links = matrix of the links to build the model
  #        - BICpartX,BICpartXY,BICpartXYZ = previously calculated partial BIC
  # OUTPUT: - bic = BIC of the whole network
  
  
  # initialize the bic
  bic = 0;
  
  # obtain the bic of each variable (marginal or conditioned)
  for (i in 1:ncol(links))
  {
    # number of parents
    nParents = sum(links[,i]);
    
    # depending on the number of parents
    if (nParents==0)   {   bic = bic + BICpartX[i];   }
    else if (nParents==1)     # 1 parent
    {
      # obtain the position of the parent
      posParent = which(links[,i]==1);
      
      # obtain the position in the list of distributions
      posfXY = i + (posParent - 1) * ncol(links);
      
      bic = bic + BICpartXY[posfXY];
    }
    else     # 2 parents
    {
      # obtain the position of both parents
      posParent1 = which(links[,i]==1)[1];
      posParent2 = which(links[,i]==1)[2];
    
      # obtain the position in the list of distributions
      posfXYZ = i + (posParent1 - 1) * ncol(links) + (posParent2 - 1) * ncol(links)^2;
    
      bic = bic + BICpartXYZ[posfXYZ];
    }
  }
  
  return(bic);
}
############################################################


















