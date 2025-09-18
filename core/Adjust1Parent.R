

Adjust1ParentDD = function(data,dataVal,C,classVal,verbose=FALSE,plots=FALSE)
{
  # it calculates the matrix of a discrete variable with a discrete parent
  # INPUT: - data = the observations of the discrete variable
  #        - dataVal = all the values of the variable in the data frame
  #        - C = the observations of the parent variable
  #        - classVal = the different values of the parent
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - dataMatrix = a matrix with the probab. of data conditioned to C
  #         - bic = BIC of the final distribution
  
  
  # create an empty matrix for the probabilities
  dataMatrix = matrix(nrow=length(dataVal),ncol=length(classVal));
  
  # add the names of the probabilities
  rownames(dataMatrix) = paste('X',dataVal,sep='_');
  colnames(dataMatrix) = paste('C',classVal,sep='_');
  
  # calculate the probability for each value of the parent
  for (i in 1:length(classVal))
  {
    # select each value of the parent and the correspondent values of data
    x = data[C==classVal[i]];
    
    # calculate the probabilities of x
    prob = DiscreteProb(x,dataVal,verbose);
    
    # put the probabilities in the table
    dataMatrix[,i] = prob[,2];
    
    if (length(x)==0)   {   cat('\n \n sin datos');   }
  }
  
  # calculate the BIC of the distribution
  bic = BIC1ParentD(dataMatrix,data,dataVal,C,classVal,varType=0);
  
  return(list(dataMatrix,bic));
}
############################################################






Adjust1ParentCD = function(data,dataMin,dataMax,C,classVal,n,maxDegree,verbose=FALSE,
                          plots=FALSE,domain=2,degree=2)
{
  # it makes the polynomial adjust of a continuous variable with a discrete parent
  # INPUT: - data = the observations of the continuous variable
  #        - dataMin, dataMax = the minimum and maximum value in the variable
  #        - C = the observations of the parent
  #        - classVal = the different values of the parent
  #        - n = position of the variable in the dataframe
  #        - maxDegree = maximum degree for the polynomial adjust
  #        - degree = 1 to adjust from degree 1 to max and 2 to adjust only the max and reduce
  #        - domain (optional) = 1 if the adjust is made in the complete domain or 2 if
  #       the adjust is made in the restricted domain and the queues extended
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - values = matrix with two rows and 3*C.val columns. It contains the
  #         min and max root, the min and max value and the probab. of the queues
  #         - polynomial = matrix with the coefficients of the polynomial  for 
  #         each value of the class
  #         - probTails = sum of the probabilities of both queues
  #         - bic = BIC of the final distribution
  
  
  # create the matrix for the polynomials and a matrix for their roots and the
  # probabilities of the values that are out of the range
  polynomial = matrix(nrow=length(classVal),ncol=maxDegree+1);
  values = matrix(0,nrow=2,ncol=length(classVal)*3);
  probTailsTotal = rep(0,length(classVal));
  
  # assign names to the three parts of values
  r = rep(c('roots','extremes','probQueue'),length(classVal));
  colnames(values) = paste(r,rep(classVal,each=3),sep='_');
  rownames(values) = c('lower','upper');
  
  # assign names to the rows of the matrix of polynomials
  rownames(polynomial) = paste('C',classVal,sep='_');
  colnames(polynomial) = paste('x',0:maxDegree,sep='^');
  
  # repeat the adjust for each value of the parent
  for (i in 1:length(classVal))
  {
    # select each value of the parent and the correspondent values of data
    x = data[C==classVal[i]]
    
    # fit polynomial and tails
    fit = PolynomialFit(x,maxDegree,dataMin,dataMax,domain,degree,verbose=FALSE,plots=FALSE)
    
    #cat('\n \n Estado de la clase ',i);
    #cat('\n Datos a ajustar entre',min(x),' y',max(x));
    #cat('\n M?nimo de los datos =',dataMin,' y m?ximo de los datos =',dataMax);
    #cat('\n Ra?ces del polinomio =',p[[1]]);
    
    if (plots==TRUE)   {
      title(paste('Variable ',n,' en la clase ',i));   }
    
    # the roots are in positions 1,4,7,.. the min and max of x in 2,5,8,...
    # and the probabilities of the tails in 3,6,9...
    values[,3*i-2] = fit[[1]];
    values[,3*i-1] = c(dataMin,dataMax);
    values[,3*i] = fit[[3]];
    
    # add the total probability of both tails
    probTailsLeft = (values[1,3*i-2] - dataMin) * values[1,3*i];
    probTailsRight = (dataMax - values[2,3*i-2]) * values[2,3*i];
    probTailsTotal[i] = probTailsLeft + probTailsRight;
    
    # the polynomials
    polynomial[i,1:length(fit[[2]])] = fit[[2]]; 
  }
  
  # calculate the BIC of the distribution
  bic = BIC1ParentD(polynomial,data,values,C,classVal,varType=1);
  
  # return the list made of roots and polynomials
  sol = list(values,polynomial,probTailsTotal,bic);
  return(sol);
}
############################################################







Adjust1ParentDC = function(data,dataVal,C,classMin,classMax,verbose=FALSE,plots=FALSE)
{
  # It calculates the probabilities of a discrete variable with one continuous
  # parent (the class) dividing it in 16 pieces
  # INPUT: - data = the observations of the discrete variable
  #        - dataVal = all the values of the variable in the data frame
  #        - C = the observations of the class
  #        - classMin, classMax = the minimum and maximum value in the class
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - prob = a matrix with the probabilities of data conditioned to C
  #         - cutPoints = a vector with the cut points in the discretization of C
  #         - bic = BIC of the final distribution
  
  
  # initialize the vector with the indices of the discretization
  discr = vector(mode='numeric',length=16);
  discr[1:8] = 1;   discr[9:16] = 2;
  
  # cut points for the discretization of the class in two intervals
  mid = (classMin + classMax) / 2;
  cutPoints = c(classMin,mid,classMax);
  
  # initialize the discretization of the class
  CDisc = vector(mode='numeric',length(C));
  
  # do the discretization of C
  CDisc[C>=cutPoints[1] & C<cutPoints[2]] = 1;
  CDisc[C>=cutPoints[2] & C<=cutPoints[3]] = 2;
  
  # the class is a discrete variable from now on
  # add all the values in case they are not represented in the discretization
  ClassVal = 1:(length(cutPoints)-1);
  
  # the function that calculates the conditional probabilities and the values and BIC
  p = Adjust1ParentDD(data,dataVal,CDisc,ClassVal,verbose,plots);
  prob = p[[1]];      bic = p[[2]];
  
  # calculate the BIC of the first discretization
  #bic = BIC1Parent(prob,data,dataVal,CDisc,ClassVal,varType=0,verbose);
  
  # initialize the decision of dividing the class
  divide = TRUE;
  
  # continue subdividing
  while (divide==TRUE)
  {
    # different values of the division (ordered as they appear in the vector)
    indices = as.numeric(levels(factor(discr,levels=unique(discr))));
    
    # number of intervals
    nIntervals = max(indices);
    
    # create a vector for the bic and the cut points of each division
    newBic = vector('numeric',nIntervals);
    newCutPoints = vector('list',nIntervals);
    p = vector('list',nIntervals);
    
    # for each interval of the discretization
    for (j in 1:nIntervals)
    {
      # insert the new subdivision in the vector of cutPoints
      newPoint = (cutPoints[j] + cutPoints[j+1]) / 2;
      newCutPoints[[j]] = sort(append(cutPoints,newPoint));
      
      # discretize the class again and calculate the number of observations
      nData = vector(mode='numeric',nIntervals+1);
      for (k in 1:(nIntervals+1))
      {
        CDisc[C>=newCutPoints[[j]][k] & C<=newCutPoints[[j]][k+1]] = k;
        nData[k] = length(C[C>=cutPoints[k] & C<=cutPoints[k+1]]);
      }
      
      # when there are less than 5% of data in one interval, that division is dismissed
      if (max(nData)<(length(C)/20))   {   newBic[j] = NA;   }
      else
      {
        # add all the values in case they are not represented in the discretization
        newClassVal = 1:(nIntervals+1);
        
        # the function that calculates the conditional probabilities and the values and BIC
        aux = Adjust1ParentDD(data,dataVal,CDisc,newClassVal,verbose,plots);
        p[[j]] = aux[[1]];        newBic[j] = aux[[2]];
        
        # calculate the BIC of the new discretization
        #newBic[j] = BIC1Parent(p[[j]],data,dataVal,CDisc,newClassVal,varType=0,verbose);
        #cat('\n newBIC = ',newBic);
      }
    }
    
    # select the best division
    best = which(newBic==max(newBic,na.rm=TRUE))[1];
    
    # check if the division increases the bic
    if (newBic[best]<=bic)   {   divide = FALSE;   }
    else
    {
      # select the cut points
      cutPoints = newCutPoints[[best]];
      
      # select the indices to divide (the second half of the indices to divide
      # is changed by a new number, which is the max+1)
      ind = which(discr==indices[best]);
      discr[ceiling((ind[1]+ind[length(ind)])/2):ind[length(ind)]] = max(discr) + 1;
      
      # select the bic
      bic = newBic[best];
      
      # select the matrix of probabilities
      prob = p[[best]];
      
      # calculate the level of the new discretization in each interval
      levelDiscr = vector(mode='numeric',nIntervals+1);
      for (k in 1:(nIntervals+1))   {   levelDiscr[k] = length(discr[discr==k]);   }
      
      # check if the level is 1, in that case, break
      if (min(levelDiscr)==1)   {   divide = FALSE;   }
    }
  }
  
  sol = list(prob,cutPoints,bic);
  return(sol);
}
############################################################






Adjust1ParentCC = function(data,dataMin,dataMax,C,classMin,classMax,n,maxDegree,
                           verbose=FALSE,plots=FALSE)
{
  # it makes the polynomial adjust of a continuous variable with a continuous parent
  # INPUT: - data = the observations of the continuous variable
  #        - dataMin = the minimum value in the variable
  #        - dataMax = the maximum value in the variable
  #        - C = the observations of the parent
  #        - classMin = the minimum value in the parent
  #        - classMax = the maximum value in the parent
  #        - n = position of the variable in the dataframe
  #        - maxDegree = maximum degree for the polynomial fit
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - values = matrix with two rows and 3*C.val columns. It contains the
  #         min and max root, the min and max value and the probab. of the queues
  #         - poly = matrix with the coefficients of the polynomial  for each value
  #         of the parent
  #         - cutPoints = a vector with the cut points in the discretization of C
  #         - bic = BIC of the final distribution
  
  
  # when there are no observations to adjust
  if (length(data)==0)
  {
    # the adjust is an uniform that integers 1
    poly = 1 / (dataMax - dataMin);
    
    # the minimum and maximum values, roots and probabilities
    valMin = c(dataMin,dataMin,0);
    valMax = c(dataMax,dataMax,0);
    
    # matrix with the values of the adjust
    values = rbind(valMin,valMax);
    
    rownames(values) = c('lower','upper');
    colnames(values) = c('roots','extremes','probQueue');
    
    # the cutPoints of the discretization of the parent are the min and max values
    cutPoints = c(classMin,classMax)
    
    return(list(values,poly,cutPoints));
  }
  
  # initialize the vector with the indices of the discretization
  discr = vector(mode='numeric',length=16);
  discr[1:8] = 1;   discr[9:16] = 2;
  
  # cut points for the discretization of the class in two intervals
  mid = (classMin + classMax) / 2;
  cutPoints = c(classMin,mid,classMax);
  
  # initialize the matrices for the polynomials of each step
  F2 = matrix(nrow=2,ncol=maxDegree+1);          F4 = matrix(nrow=4,ncol=maxDegree+1);
  F8 = matrix(nrow=8,ncol=maxDegree+1);          F16 = matrix(nrow=16,ncol=maxDegree+1);
  
  # initialize the matrices for the roots and prob of the queues of each step
  A2 = matrix(nrow=2,ncol=3*2);         A4 = matrix(nrow=2,ncol=3*4);
  A8 = matrix(nrow=2,ncol=3*8);         A16 = matrix(nrow=2,ncol=3*16);
  
  # initialize the discretization of the class
  CDisc = vector(mode='numeric',length(C));
  
  # do the discretization of C
  CDisc[C>=cutPoints[1] & C<cutPoints[2]] = 1;
  CDisc[C>=cutPoints[2] & C<=cutPoints[3]] = 2;
  
  # the class is a discrete variable from now on
  # add all the values in case they are not represented in the discretization
  classVal = 1:(length(cutPoints)-1);
  
  # the function that calculates the conditional probabilities and the values
  prob = Adjust1ParentCD(data,dataMin,dataMax,CDisc,classVal,n,maxDegree,verbose,plots);
  poly = prob[[2]];   values = prob[[1]];
  
  # fill in the first level matrices
  F2 = prob[[2]];    A2 = prob[[1]];
  
  # calculate the BIC of the first discretization
  bic = BIC1ParentD(F2,data,A2,CDisc,classVal,varType=1,verbose);
  
  # initialize the decision of dividing the class
  divide = TRUE;
  
  # continue subdividing
  while (divide==TRUE)
  {
    # different values of the division (ordered as they appear in the vector)
    indices = as.numeric(levels(factor(discr,levels=unique(discr))));
    
    # number of intervals
    nIntervals = max(indices);
    
    # create a vector for the bic, cut points, level and adjust of each division
    newBic = vector('numeric',nIntervals);
    newCutPoints = vector('list',nIntervals);
    newPoly = vector('list',nIntervals);
    newValues = vector('list',nIntervals);
    
    # create a matrix with discr in each row to modify in each subdivision
    newDiscr = matrix(discr,nrow=nIntervals,ncol=length(discr),byrow=TRUE);
    
    # for each interval of the discretization
    for (j in 1:nIntervals)
    {
      # insert the new subdivision in the vector of cutPoints
      newPoint = (cutPoints[j] + cutPoints[j+1]) / 2;
      newCutPoints[[j]] = sort(append(cutPoints,newPoint));
      
      # discretize the class again and calculate the number of observations
      nData = vector(mode='numeric',nIntervals+1);
      for (k in 1:(nIntervals+1))
      {
        CDisc[C>=newCutPoints[[j]][k] & C<=newCutPoints[[j]][k+1]] = k;
        nData[k] = length(CDisc[CDisc==k]);
      }
      
      # when there are less than 5% of data in one interval, that division is dismissed
      if (max(nData)<(length(C)/20))   {   newBic[j] = NA;   }
      else
      {
        # add all the values in case they are not represented in the discretization
        newClassVal = 1:(nIntervals+1);
        
        # select the indices to divide (the second half of the indices to divide
        # is changed by a new number, which is max+1)
        ind = which(discr==indices[j]);
        newDiscr[j,ceiling((ind[1]+ind[length(ind)])/2):ind[length(ind)]] = max(discr) + 1;
        
        # different values of the new division (ordered as they appear in the vector)
        newIndices = as.numeric(levels(factor(newDiscr[j,],levels=unique(newDiscr[j,]))));
        
        # create the matrices for the adjust dividing the current interval
        polyj = matrix(nrow=nIntervals+1,ncol=maxDegree+1);
        valuesj = matrix(nrow=2,ncol=3*(nIntervals+1));
        
        # for each interval in the new division
        for (k in 1:(nIntervals+1))
        {
          # calculate the level of the discretization and the beginning of the interval
          x = newDiscr[j,];
          levelDiscr = length(x[x==newIndices[k]]);
          start = which(x==newIndices[k])[1];
          
          # data to adjust
          d = data[CDisc==k];      r = rep(0,length(d));
          
          # check if a previous adjust can be used again
          if (levelDiscr==8 && start==1)     # first of 2 pieces of 8 elements
          {
            # the adjust was previously done
            polyj[k,] = F2[1,];       valuesj[,(3*k-2):(3*k)] = A2[,1:3];
          }
          else if (levelDiscr==8 && start==9)     # second of 2 pieces of 8 elements
          {
            # the adjust was previously done
            polyj[k,] = F2[2,];       valuesj[,(3*k-2):(3*k)] = A2[,4:6];   
          }
          else if (levelDiscr==4 && start==1)     # first of 4 pieces of 4 elements
          {
            # when the adjust was previously done
            if (is.na(F4[1,1])==FALSE)
            {
              polyj[k,] = F4[1,];       valuesj[,(3*k-2):(3*k)] = A4[,1:3];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F4[1,] = prob[[2]];       A4[,1:3] = prob[[1]];
            }
          }
          else if (levelDiscr==4 && start==5)     # second of 4 pieces of 4 elements
          {
            # when the adjust was previously done
            if (is.na(F4[2,1])==FALSE)
            {
              polyj[k,] = F4[2,];       valuesj[,(3*k-2):(3*k)] = A4[,4:6];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F4[2,] = prob[[2]];       A4[,4:6] = prob[[1]];
            }
          }
          else if (levelDiscr==4 && start==9)     # third of 4 pieces of 4 elements
          {
            # when the adjust was previously done
            if (is.na(F4[3,1])==FALSE)
            {
              polyj[k,] = F4[3,];       valuesj[,(3*k-2):(3*k)] = A4[,7:9];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F4[3,] = prob[[2]];       A4[,7:9] = prob[[1]];
            }
          }
          else if (levelDiscr==4 && start==13)     # fourth of 4 pieces of 4 elements
          {
            # when the adjust was previously done
            if (is.na(F4[4,1])==FALSE)
            {
              polyj[k,] = F4[4,];       valuesj[,(3*k-2):(3*k)] = A4[,10:12];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F4[4,] = prob[[2]];       A4[,10:12] = prob[[1]];
            }
          }
          else if (levelDiscr==2 && start==1)     # first of 8 pieces of 2 elements
          {
            # when the adjust was previously done
            if (is.na(F8[1,1])==FALSE)
            {
              polyj[k,] = F8[1,];       valuesj[,(3*k-2):(3*k)] = A8[,1:3];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F8[1,] = prob[[2]];       A8[,1:3] = prob[[1]];
            }
          }
          else if (levelDiscr==2 && start==3)     # second of 8 pieces of 2 elements
          {
            # when the adjust was previously done
            if (is.na(F8[2,1])==FALSE)
            {
              polyj[k,] = F8[2,];       valuesj[,(3*k-2):(3*k)] = A8[,4:6];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F8[2,] = prob[[2]];       A8[,4:6] = prob[[1]];
            }
          }
          else if (levelDiscr==2 && start==5)     # third of 8 pieces of 2 elements
          {
            # when the adjust was previously done
            if (is.na(F8[3,1])==FALSE)
            {
              polyj[k,] = F8[3,];       valuesj[,(3*k-2):(3*k)] = A8[,7:9];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F8[3,] = prob[[2]];       A8[,7:9] = prob[[1]];
            }
          }
          else if (levelDiscr==2 && start==7)     # fourth of 8 pieces of 2 elements
          {
            # when the adjust was previously done
            if (is.na(F8[4,1])==FALSE)
            {
              polyj[k,] = F8[4,];       valuesj[,(3*k-2):(3*k)] = A8[,10:12];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F8[4,] = prob[[2]];       A8[,10:12] = prob[[1]];
            }
          }
          else if (levelDiscr==2 && start==9)     # fifth of 8 pieces of 2 elements
          {
            # when the adjust was previously done
            if (is.na(F8[5,1])==FALSE)
            {
              polyj[k,] = F8[5,];       valuesj[,(3*k-2):(3*k)] = A8[,13:15];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F8[5,] = prob[[2]];       A8[,13:15] = prob[[1]];
            }
          }
          else if (levelDiscr==2 && start==11)     # sixth of 8 pieces of 2 elements
          {
            # when the adjust was previously done
            if (is.na(F8[6,1])==FALSE)
            {
              polyj[k,] = F8[6,];       valuesj[,(3*k-2):(3*k)] = A8[,16:18];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F8[6,] = prob[[2]];       A8[,16:18] = prob[[1]];
            }
          }
          else if (levelDiscr==2 && start==13)     # seventh of 8 pieces of 2 elements
          {
            # when the adjust was previously done
            if (is.na(F8[7,1])==FALSE)
            {
              polyj[k,] = F8[7,];       valuesj[,(3*k-2):(3*k)] = A8[,19:21];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F8[7,] = prob[[2]];       A8[,19:21] = prob[[1]];
            }
          }
          else if (levelDiscr==2 && start==15)     # eighth of 8 pieces of 2 elements
          {
            # when the adjust was previously done
            if (is.na(F8[8,1])==FALSE)
            {
              polyj[k,] = F8[8,];       valuesj[,(3*k-2):(3*k)] = A8[,22:24];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F8[8,] = prob[[2]];       A8[,22:24] = prob[[1]];
            }
          }
          else if (levelDiscr==1 && start==1)     # first of 16 pieces of 1 element
          {
            # when the adjust was previously done
            if (is.na(F16[1,1])==FALSE)
            {
              polyj[k,] = F16[1,];       valuesj[,(3*k-2):(3*k)] = A16[,1:3];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F16[1,] = prob[[2]];       A16[,1:3] = prob[[1]];
            }
          }
          else if (levelDiscr==1 && start==2)     # second of 16 pieces of 1 element
          {
            # when the adjust was previously done
            if (is.na(F16[2,1])==FALSE)
            {
              polyj[k,] = F16[2,];       valuesj[,(3*k-2):(3*k)] = A16[,4:6];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F16[2,] = prob[[2]];       A16[,4:6] = prob[[1]];
            }
          }
          else if (levelDiscr==1 && start==3)     # third of 16 pieces of 1 element
          {
            # when the adjust was previously done
            if (is.na(F16[3,1])==FALSE)
            {
              polyj[k,] = F16[3,];       valuesj[,(3*k-2):(3*k)] = A16[,7:9];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F16[3,] = prob[[2]];       A16[,7:9] = prob[[1]];
            }
          }
          else if (levelDiscr==1 && start==4)     # fourth of 16 pieces of 1 element
          {
            # when the adjust was previously done
            if (is.na(F16[4,1])==FALSE)
            {
              polyj[k,] = F16[4,];       valuesj[,(3*k-2):(3*k)] = A16[,10:12];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F16[4,] = prob[[2]];       A16[,10:12] = prob[[1]];
            }
          }
          else if (levelDiscr==1 && start==5)     # fifth of 16 pieces of 1 element
          {
            # when the adjust was previously done
            if (is.na(F16[5,1])==FALSE)
            {
              polyj[k,] = F16[5,];       valuesj[,(3*k-2):(3*k)] = A16[,13:15];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F16[5,] = prob[[2]];       A16[,13:15] = prob[[1]];
            }
          }
          else if (levelDiscr==1 && start==6)     # sixth of 16 pieces of 1 element
          {
            # when the adjust was previously done
            if (is.na(F16[6,1])==FALSE)
            {
              polyj[k,] = F16[6,];       valuesj[,(3*k-2):(3*k)] = A16[,16:18];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F16[6,] = prob[[2]];       A16[,16:18] = prob[[1]];
            }
          }
          else if (levelDiscr==1 && start==7)     # seventh of 16 pieces of 1 element
          {
            # when the adjust was previously done
            if (is.na(F16[7,1])==FALSE)
            {
              polyj[k,] = F16[7,];       valuesj[,(3*k-2):(3*k)] = A16[,19:21];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F16[7,] = prob[[2]];       A16[,19:21] = prob[[1]];
            }
          }
          else if (levelDiscr==1 && start==8)     # eighth of 16 pieces of 1 element
          {
            # when the adjust was previously done
            if (is.na(F16[8,1])==FALSE)
            {
              polyj[k,] = F16[8,];       valuesj[,(3*k-2):(3*k)] = A16[,22:24];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F16[8,] = prob[[2]];       A16[,22:24] = prob[[1]];
            }
          }
          else if (levelDiscr==1 && start==9)     # ninth of 16 pieces of 1 element
          {
            # when the adjust was previously done
            if (is.na(F16[9,1])==FALSE)
            {
              polyj[k,] = F16[9,];       valuesj[,(3*k-2):(3*k)] = A16[,25:27];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F16[9,] = prob[[2]];       A16[,25:27] = prob[[1]];
            }
          }
          else if (levelDiscr==1 && start==10)     # 10th of 16 pieces of 1 element
          {
            # when the adjust was previously done
            if (is.na(F16[10,1])==FALSE)
            {
              polyj[k,] = F16[10,];       valuesj[,(3*k-2):(3*k)] = A16[,28:30];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F16[10,] = prob[[2]];       A16[,28:30] = prob[[1]];
            }
          }
          else if (levelDiscr==1 && start==11)     # 11th of 16 pieces of 1 element
          {
            # when the adjust was previously done
            if (is.na(F16[11,1])==FALSE)
            {
              polyj[k,] = F16[11,];       valuesj[,(3*k-2):(3*k)] = A16[,31:33];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F16[11,] = prob[[2]];       A16[,31:33] = prob[[1]];
            }
          }
          else if (levelDiscr==1 && start==12)     # 12th of 16 pieces of 1 element
          {
            # when the adjust was previously done
            if (is.na(F16[12,1])==FALSE)
            {
              polyj[k,] = F16[12,];       valuesj[,(3*k-2):(3*k)] = A16[,34:36];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F16[12,] = prob[[2]];       A16[,34:36] = prob[[1]];
            }
          }
          else if (levelDiscr==1 && start==13)     # 13th of 16 pieces of 1 element
          {
            # when the adjust was previously done
            if (is.na(F16[13,1])==FALSE)
            {
              polyj[k,] = F16[13,];       valuesj[,(3*k-2):(3*k)] = A16[,37:39];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F16[13,] = prob[[2]];       A16[,37:39] = prob[[1]];
            }
          }
          else if (levelDiscr==1 && start==14)     # 14th of 16 pieces of 1 element
          {
            # when the adjust was previously done
            if (is.na(F16[14,1])==FALSE)
            {
              polyj[k,] = F16[14,];       valuesj[,(3*k-2):(3*k)] = A16[,40:42];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F16[14,] = prob[[2]];       A16[,40:42] = prob[[1]];
            }
          }
          else if (levelDiscr==1 && start==15)     # 15th of 16 pieces of 1 element
          {
            # when the adjust was previously done
            if (is.na(F16[15,1])==FALSE)
            {
              polyj[k,] = F16[15,];       valuesj[,(3*k-2):(3*k)] = A16[,43:45];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F16[15,] = prob[[2]];       A16[,43:45] = prob[[1]];
            }
          }
          else if (levelDiscr==1 && start==16)     # 16th of 16 pieces of 1 element
          {
            # when the adjust was previously done
            if (is.na(F16[16,1])==FALSE)
            {
              polyj[k,] = F16[16,];       valuesj[,(3*k-2):(3*k)] = A16[,46:48];   
            }
            else    # the adjust is not done
            {
              prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
              polyj[k,] = prob[[2]];       valuesj[,(3*k-2):(3*k)] = prob[[1]];
              F16[16,] = prob[[2]];        A16[,46:48] = prob[[1]];
            }
          }
        }
        
        # calculate the BIC of the new discretization
        newBic[j] = BIC1ParentD(polyj,data,valuesj,CDisc,newClassVal,varType=1,verbose);
        
        # fill in the list with polynomials and values of different discretizations
        newPoly[[j]] = polyj;     newValues[[j]] = valuesj;
      }
    }
    
    # select the best division
    best = which(newBic==max(newBic,na.rm=TRUE))[1];
    
    # check if the division increases the bic
    if (newBic[best]<=bic)   {   divide = FALSE;   }
    else
    {
      # select the cut points and the positions in the discretization
      cutPoints = newCutPoints[[best]];     discr = newDiscr[best,];
      
      # select the best bic
      bic = newBic[best];
      
      # select the best polynomial and values
      poly = newPoly[[best]];     values = newValues[[best]]
      
      # names for the rows and columns of the matrices
      rownames(values) = c('lower','upper');
      colnames(values) = rep(c('roots','extremes','probQueue'),nIntervals+1);
      rownames(poly) = paste('int',1:(nIntervals+1));
      colnames(poly) = paste('x',0:maxDegree,sep='^');
      
      # calculate the level of the new discretization in each interval
      levelDiscr = vector(mode='numeric',nIntervals+1);
      
      for (k in 1:(nIntervals+1))   {   levelDiscr[k] = length(discr[discr==k]);   }
      
      # check if the level is 1, in that case, break
      if (min(levelDiscr)==1)   {   divide = FALSE;   }
    }
  }
  
  sol = list(values,poly,cutPoints,bic);
  return(sol);
}
############################################################












