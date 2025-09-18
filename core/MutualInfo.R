
# CARGAR DEPENDENCIAS AL PRINCIPIO (discretization)
library(infotheo);
library(discretization);
library(classInt);



MutualInfo = function(data,varType,verbose=FALSE)
{
  # calculate the mutual information of each variable and the class using the package infotheo
  # INPUT: - data = a data-frame with the data of study
  #        - varType = a vector with 0 if the variable is discrete, 1 if it's
  #       continuous and 2 if it has to be discretizated. Program calculates it
  #        - verbose = show the 'cats' on screen
  # OUTPUT: - MI: vector with the mutual information of each variables and the class
  
  
  # number of variables
  nVar = dim(data)[2];
  
  # make the discretization of all the continuous variables (discretized are already done)
  disc = DiscretizationMI(data,varType,verbose);
  data = disc[[1]];
  dataVal = disc[[2]];
  
  # infotheo needs discrete values
  for (i in 1:length(varType))   {   data[,i] = as.factor(data[,i]);   }
  
  # create the vector to fill in with the MI
  MI = vector(mode='numeric',length=nVar-1);
  
  # for each observation variable
  for (i in 1:(nVar-1))
  {
    MI[i] = infotheo :: condinformation(data[,i],data[,nVar],method='emp');
  }
  
  return(MI);
}
############################################################








MutualInfoClass = function(data,varType,verbose)
{
  # calculate the mutual information of each pair of variables conditioned to the class
  # using the package infotheo
  # INPUT: - data = a data-frame with the data of study
  #        - varType = a vector with 0 if the variable is discrete, 1 if it's
  #       continuous and 2 if it has to be discretizated. Program calculates it
  # OUTPUT: - MIC = matrix with the mutual information of each pair of variables
  
  
  # number of variables
  nVar = dim(data)[2];
  
  # make the discretization of all the continuous variables (discretizable are done
  # in the main program)
  disc = DiscretizationMI(data,varType,verbose);
  data = disc[[1]];
  dataVal = disc[[2]];
  
  # infotheo needs discrete values
  for (i in 1:length(varType))   {   data[,i] = as.factor(data[,i]);   }
  
  # create the matrix to fill in with the MIC
  MIC = matrix(-9999,nVar-1,nVar-1);
  
  # for each pair of variables
  for (i in 1:(nVar-2))
  {
    for (j in (i+1):(nVar-1))
    {
      MIC[j,i] = infotheo :: condinformation(data[,i],data[,j],data[,nVar],method='emp');
    }
  }
  
  return(MIC);
}
############################################################







MutualInfoGeneral = function(data,varType,verbose=FALSE)
{
  # calculate the mutual information of each pair of variables and also the mutual 
  # information of each pair conditioned to a third variable using the package infotheo
  # INPUT: - data = a data-frame with the data of study
  #        - varType = a vector with 0 if the variable is discrete, 1 if it's
  #       continuous and 2 if it has to be discretizated. Program calculates it
  #        - verbose = show the 'cats' on screen
  # OUTPUT: - MI: vector with the mutual information of each variables and the class
  #         - MIC = matrix with the mutual information of each pair of variables
  
  
  # number of variables
  nVar = ncol(data);
  
  # make the discretization of all the continuous variables (discretized are already done)
  disc = DiscretizationMI2(data,varType,verbose);
  data = disc[[1]];
  dataVal = disc[[2]];
  
  # infotheo needs discrete values
  for (i in 1:length(varType))   {   data[,i] = as.factor(data[,i]);   }
  
  # create the matrix to fill in with the MI initialited to 0
  MI = matrix(-9999,nrow=nVar,ncol=nVar);
  
  # for each pair of variables
  for (i in 1:(nVar-1))
  {
    for (j in (i+1):nVar)
    {
      MI[j,i] = infotheo :: condinformation(data[,i],data[,j],method='emp');
    }
  }
  
  
  
  # create the 3d-matrix to fill in with the MIC
  MIC = array(-9999,dim=c(nVar,nVar,nVar));
  
  aux = 1:nVar;
  
  # for each pair of variables and each variable to make the condition
  for (k in 1:nVar)
  {
    for (i in aux[aux!=k & aux<nVar])   # avoid k and the last variable
    {
      for (j in aux[aux!=k & aux>i])   # avoid k and start from i+1
      {
        MIC[j,i,k] = infotheo :: condinformation(data[,i],data[,j],data[,k],method='emp');
      }
    }
  }
  
  sol = list(MI,MIC)
  return(sol);
}
############################################################









DiscretizationMI = function(data,varType,verbose=FALSE)
{
  # it makes the discretization of all continuous variables in the dataframe and
  # calculates the values of each variable. If the class is continuous (regression), 
  # it also discretizes it using an equal frequency method
  # INPUT: - data = a data-frame with the data of study
  #        - varType = a vector with 0 if the variable is discrete, 1 if it's
  #       continuous and 2 if it has to be discretizated. Program calculates it
  # OUTPUT: - data = a data frame with the discretized data
  #         - dataVal = all the values of the variable in the data frame
  
  
  nVar = dim(data)[2];
  
  # do the discretization of the class using equal width method (only for mdlp)
  if (varType[nVar]==1)
  {
    data[,nVar] = arules :: discretize(data[,nVar],method='frequency',labels=FALSE);
    # from now on the class is discrete
    varType[nVar] = 0;
  }
  
  # change the continuous variables to pseudo-discrete and pseudo-discrete to discrete
  aux = varType;
  varType[aux==1] = 2;    varType[aux==2] = 0;
  
  # change the class by the discretized one
  discret = Discretization(data,varType,2);
  
  # the class isn't discretized
  data[,1:(nVar-1)] = discret[[1]][,1:(nVar-1)];
  cut.points = discret[[2]];
  
  # calculate the values of each variable
  dataVal = vector(mode='list',length=nVar);
  for (i in 1:nVar)
  {
    dataVal[[i]] = rle(sort(data[,i]))$values;
  }
  
  # return the discretized data and the values
  sol = list(data,dataVal)
  return(sol);
}
############################################################









DiscretizationMI2 = function(data,varType,verbose=FALSE)
{
  # it makes the discretization of all continuous variables in the dataframe and
  # calculates the values of each variable with an equal with method
  # INPUT: - data = a data-frame with the data of study
  #        - varType = a vector with 0 if the variable is discrete, 1 if it's
  #        continuous and 2 if it has to be discretized. Program calculates it
  # OUTPUT: - data = a data frame with the discretized data
  #         - dataVal = all the values of the variable in the data frame
  
  
  nVar = ncol(data);
  
  # maximum and minimum of each variable
  dataMin = as.numeric(sapply(data,min)) - 0.000000001;
  dataMax = as.numeric(sapply(data,max)) + 0.000000001;
  
  # index of continuous and pseudo-discrete variables
  indexDisc = which(varType>0);
  nDisc = length(indexDisc);
  
  # check if there're discretizable variables
  if (nDisc==0)      # nothing to discretize
  {
    indexDisc = 0;        cutPoints = 0;
  }
  else
  {
    # initialize the split points
    splitPoints = vector(mode='list',length=nDisc);
    
    for (i in 1:nDisc)
    {
      # calculate the cut points for the discretization
      cutPoints = classInt :: classIntervals(data[,indexDisc[i]],style='equal')$brks;
      
      # change the first point
      cutPoints[1] = cutPoints[1] - 0.00001;
      
      # make the discretization
      data[,indexDisc[i]] = cut(data[,indexDisc[i]],breaks=cutPoints,labels=FALSE);
      
      # quit the first and last cut points of the discretization
      if (length(cutPoints) > 2)
      {
        splitPoints[[i]] = cutPoints[2:(length(cutPoints)-1)];
      }
    }
    
    # add the minimum and maximum to the split points and check if there are variables 
    # with only a factor in the discretization
    data2 = data;
    for (i in 1:nDisc)
    {
      if (length(which(splitPoints[[i]]=='All')))
      {
        # discretize using the median (1 < median <= 2)
        data2[,indexDisc[i]] = rep(2,nrow(data));
        
        # median of the data
        m = median(data[,indexDisc[i]]);
        
        # check if the median is equal to any of the extremes
        if (dataMin[indexDisc[i]]==m)   {   m = m + 0.00001;   }
        else if (dataMax[indexDisc[i]]==m)   {   m = m - 0.00001;   }
        
        data2[data[,indexDisc[i]] < m, indexDisc[i]] = 1;
        
        # the split points are the median, the maximum and the minimum
        splitPoints[[i]] = c(dataMin[indexDisc[i]],m,dataMax[indexDisc[i]]);
      }
      else
      {
        # the vector of the split points includes minimum and maximum
        splitPoints[[i]] = c(dataMin[indexDisc[i]],splitPoints[[i]],dataMax[indexDisc[i]]);
      }
    }
    
    # rename the dataset
    data = data2;
    
    # create a list for the split points of all the variables
    cutPoints = vector(mode='list',nVar);
    cutPoints[indexDisc] = splitPoints;
  }
  
  sol = list(data,cutPoints);
  return(sol);
}
############################################################






