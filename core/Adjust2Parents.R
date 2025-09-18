

Adjust2ParentsDDD = function(data,dataVal,parent,parentVal,C,classVal,verbose=FALSE,
                             plots=FALSE)
{
  # It calculates the probabilities of a discrete variable with two discrete parents
  # the class and another discrete variable
  # INPUT: - data = the observations of the discrete variable
  #        - dataVal = all the values of the variable in the data frame
  #        - parent = the observation of the parent
  #        - parentVal = all the values of the parent
  #        - C = the observations of the class
  #        - classVal = the different values of the class
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - dataMatrix = a matrix 3D with the probabilities of data conditioned to 
  #         both parents
  #         - bic = BIC of the final distribution
  
  
  # create an empty 3d-matrix for the probabilities
  dataMatrix = array(dim=c(length(dataVal),length(parentVal),length(classVal)));
  
  # assign names to the three dimensions of the array
  dimnames(dataMatrix)[[1]] = paste('X',dataVal,sep='_');
  dimnames(dataMatrix)[[2]] = paste('p1',parentVal,sep='_');
  dimnames(dataMatrix)[[3]] = paste('C',classVal,sep='_');
  
  # calculate the probability for each value of the parent and the class
  for (i in 1:length(classVal))
  {
    for (j in 1:length(parentVal))
    {
      # filter the data by each value of the parent and the class
      x = data[C==classVal[i] & parent==parentVal[j]];
      
      # calculate the probabilities of x
      prob = DiscreteProb(x,dataVal);
      
      # put the probabilities in the table
      dataMatrix[,j,i] = prob[,2];
    }
  }
  
  # calculate the BIC of the distribution
  bic = BIC2ParentsDD(dataMatrix,data,dataVal,parent,parentVal,C,classVal,varType=0,verbose);
  
  return(list(dataMatrix,bic));
}
############################################################





Adjust2ParentsCDD = function(data,dataMin,dataMax,parent,parentVal,C,classVal,n,
                             maxDegree,verbose=FALSE,plots=FALSE)
{
  # It makes the polynomial adjust of a continuous variable with two discrete parents
  # the class and another discrete variable
  # INPUT: - data = the observations of the continuous variable
  #        - dataMin, dataMax = the minimum and maximum value in the variable
  #        - parent = the observation of the parent
  #        - parentVal = all the values of the parent
  #        - C = the observations of the parent
  #        - classVal = the different values of the parent
  #        - n = position of the variable in the dataframe
  #        - maxDegree = maximum degree for the polynomial fit
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - values = matrix with 2*P rows and 3*C.val columns. It contains the
  #         min and max root, the min and max value and the probab. of the tails
  #         - polynomial = matrix with the coefficients of the polynomial  for 
  #         each value of the class
  #         - probQueues = sum of the probabilities of both queues
  #         - bic = BIC of the final distribution
  
  
  # create the list for the polynomials and a matrix for their roots and the
  # probabilities of the values that are out of the range
  polynomial = array(dim=c(length(parentVal), maxDegree+1, length(classVal)));
  values = matrix(nrow=length(parentVal)*2, ncol=length(classVal)*3);
  probQueues = matrix(nrow=length(parentVal), ncol=length(classVal));
  
  # assign names to the three dimensions of polynomial
  dimnames(polynomial)[[1]] = paste('p1',parentVal,sep='_');
  dimnames(polynomial)[[2]] = paste('x',0:maxDegree,sep='^');
  dimnames(polynomial)[[3]] = paste('C',classVal,sep='_');
  
  # assign names to the three dimensions of values
  r1 = rep(c('lower','upper'),length(parentVal));
  r2 = rep(c('roots','extremes','probQueue'),length(classVal));
  rownames(values) = paste(r1,rep(parentVal,each=2),sep='_');
  colnames(values) = paste(r2,rep(classVal,each=3),sep='_');
  
  # adjust polynomials for each parent
  for (i in 1:length(parentVal))
  {
    # filter the data and the class by each value of the parent
    x = data[parent==parentVal[i]];
    class = C[parent==parentVal[i]];
    
    # adjust 
    p = Adjust1ParentCD(x,dataMin,dataMax,class,classVal,n,maxDegree,verbose,plots);
    
    # put the polynomials and values in the matrices
    values[(2*i-1):(2*i),] = p[[1]];
    polynomial[i,,] = t(p[[2]]);
    probQueues[i,] = p[[3]];
  }
  
  # calculate the BIC of the distribution
  bic = BIC2ParentsDD(polynomial,data,values,parent,parentVal,C,classVal,varType=1,verbose);
    
  # return the list made of roots, polynomials and BIC
  sol = list(values,polynomial,probQueues,bic);
  return(sol);
}
############################################################








Adjust2ParentsDCD = function(data,dataVal,parent,parentMin,parentMax,C,classVal,
                             verbose=FALSE,plots=FALSE)
{
  # It calculates the probabilities of a discrete variable with one discrete and one
  # continuous parent
  # INPUT: - data = the observations of the discrete variable
  #        - dataVal = all the values of the variable in the data frame
  #        - parent = the observation of the continuous parent
  #        - parentMin, parentMax = the minimum and maximum value in the parent
  #        - C = the observations of the discrete class
  #        - classVal = the different values of the discrete class
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - dataMatrix = a matrix 3D with the probabilities of data conditioned to 
  #         both parents
  #         - cutPoints = matrix with the cut points of each discretization
  #         - bic = BIC of the final distribution
  
  
  # create an empty 3d-matrix for the probabilities
  dataMatrix = array(dim=c(length(dataVal), 16, length(classVal)));
  
  # assign names to the three dimensions of the array of probabilities
  dimnames(dataMatrix)[[1]] = paste('X',dataVal,sep='_');
  dimnames(dataMatrix)[[2]] = paste('p1',1:16,sep='_');
  dimnames(dataMatrix)[[3]] = paste('p2',classVal,sep='_');
  
  # create an empty matrix for the cutPoints of the discretizations
  cutPoints = matrix(nrow=length(classVal),ncol=16);
  
  # create an empty vector for the number of intervals in each discretization
  nIntervals = vector(mode='numeric',length(classVal));
  
  # initialize the BIC
  bic = 0;
  
  # calculate the probability for each value of the class
  for (i in 1:length(classVal))
  {
    # filter the data by each value of the class
    x = data[C==classVal[i]];        p = parent[C==classVal[i]];
    
    # calculate the probabilities of x conditioned to the continuous parent
    prob = Adjust1ParentDC(x,dataVal,p,parentMin,parentMax,verbose,plots);
    
    # the number of intervals
    nIntervals[i] = length(prob[[2]]) - 1;
    
    # put the probabilities in the table
    dataMatrix[,1:nIntervals[i],i] = prob[[1]];
    
    # put the cutPoints in the matrix
    cutPoints[i,1:(nIntervals[i]+1)] = prob[[2]];
    
    # update the BIC checking if there are some data
    if (length(x)>0)   {   bic = bic + prob[[3]];   }
    #bic = bic + prob[[3]];
  }
  
  # reduce the matrices to the maximum number of intervals in the discretizations
  n = max(nIntervals);
  dataMatrix = array(dataMatrix[,1:n,],dim=c(dim(dataMatrix)[1],n,dim(dataMatrix)[3]));
  cutPoints = cutPoints[,1:(n+1)];
  
  # subtract the number of parameters in classVal
  #bic = bic - (length(classVal) - 1) / 2 * log10(length(data));
  
  sol = list(dataMatrix,cutPoints,bic);
  return(sol);
}
############################################################






Adjust2ParentsCCD = function(data,dataMin,dataMax,parent,parentMin,parentMax,C,classVal,
                             n,maxDegree,verbose=FALSE,plots=FALSE)
{
  # It calculates the probabilities of a continuous variable with one discrete and one
  # continuous parent (parent is continuous and C is discrete)
  # INPUT: - data = the observations of the continuous variable
  #        - dataMin, dataMax = the minimum and maximum value in the variable
  #        - parent = the observation of the continuous parent
  #        - parentMin, parentMax = the minimum and maximum value in the parent
  #        - C = the observations of the discrete class
  #        - classVal = the different values of the discrete class
  #        - n = position of the variable in the dataframe
  #        - maxDegree = maximum degree for the polynomial fit
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - values = matrix with 2*P and 3*C.val columns. It contains the
  #         min and max root, the min and max value and the probab. of the tails
  #         - polynomial = matrix with the coefficients of the polynomial  for 
  #         each value of the class and the discretization of the other parent
  #         - cutPoints = matrix with the cut points of each discretization
  #         - bic = BIC of the final distribution
  
  
  # create an empty 3d-matrix for the probabilities
  polynomial = array(dim=c(16, maxDegree+1, length(classVal)));
  values = matrix(nrow=length(classVal)*2, ncol=3*16);
  
  # assign names to the three dimensions of polynomial
  dimnames(polynomial)[[1]] = paste('p1',1:16,sep='_');
  dimnames(polynomial)[[2]] = paste('x',0:maxDegree,sep='^');
  dimnames(polynomial)[[3]] = paste('p2',classVal,sep='_');
  
  # assign names to the three dimensions of values
  r1 = rep(c('lower','upper'),length(classVal));
  r2 = rep(c('roots','extremes','probQueue'),16);
  rownames(values) = paste(r1,rep(classVal,each=2),sep='_');
  colnames(values) = paste(r2,rep(1:16,each=3),sep='_');
  
  # create an empty matrix for the cutPoints of the discretizations
  cutPoints = matrix(nrow=length(classVal), ncol=16);
  
  # create an empty vector for the number of intervals in each discretization
  nIntervals = vector(mode='numeric',length(classVal));
  
  # initialize the BIC
  bic = 0;
  
  if (verbose==TRUE)   {
    cat('\n En Ajuste2ParentsCCD');
    cat('\n C = ',C);
    cat('\n classVal = ',classVal);
    cat('\n dataMin = ',dataMin,'y dataMax = ',dataMax);
  }
  
  # calculate the probability for each value of the class
  for (i in 1:length(classVal))
  {
    # filter the data by each value of the class
    x = data[C==classVal[i]];        p = parent[C==classVal[i]];
    
    # adjust the polynomials for x conditioned to the continuous parent
    prob = Adjust1ParentCC(x,dataMin,dataMax,p,parentMin,parentMax,n,maxDegree,verbose,
                           plots);
    nIntervals[i] = length(prob[[3]]) - 1;
    
    if (verbose==TRUE)   {
      cat('\n datos = ',x);  
      cat('\n values 2 = ',dim(prob[[1]]));
      cat('\n polynomial 2 = ',dim(prob[[2]]));
      cat('\n cutPoints 2 = ',length(prob[[3]]));
      cat('\n prob4 = ',prob[[4]]);
    }
    
    # the values and polynomials
    values[(2*i-1):(2*i),1:(3*nIntervals[i])] = prob[[1]];
    polynomial[1:nIntervals[i],,i] = prob[[2]];
    
    # the cut points of the discretization of the continuous parent
    cutPoints[i,1:(nIntervals[i]+1)] = prob[[3]];
    
    # update the BIC checking if there are some data
    #cat('\n \n \n va a fallar en el BIC de 2 padres CCD')
    if (length(x)>0)   {   bic = bic + prob[[4]];   }
  }
  
  # subtract the number of parameters in classVal
  #bic = bic - (length(classVal) - 1) / 2 * log10(length(data));
  
  sol = list(values,polynomial,cutPoints,bic);
  return(sol);
}
############################################################









Adjust2ParentsDCC = function(data,dataVal,parent,parentMin,parentMax,C,classMin,
                             classMax,n,maxDegree,verbose=FALSE,plots=FALSE)
{
  # It calculates the probabilities of a discrete variable with two continuous parents
  # INPUT: - data = the observations of the continuous variable
  #        - dataVal = all the values of the variable in the data frame
  #        - parent = the observation of the continuous parent
  #        - parentMin, parentMax = the minimum and maximum value in the parent
  #        - C = the observations of the discrete class
  #        - classMin, classMax = the minimum and maximum value in the class
  #        - n = position of the variable in the dataframe
  #        - maxDegree = maximum degree for the polynomial adjust
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - prob = a matrix with the probabilities of data conditioned to C
  #         - cutPoints = matrix with the cut points of each discretization
  #         - bic = BIC of the final distribution
  
  
  # initialize the matrix with the indices of both discretizations
  discr = matrix(nrow=8,ncol=8);
  discr[1:4,1:4] = 1;   discr[1:4,5:8] = 2;
  discr[5:8,1:4] = 3;   discr[5:8,5:8] = 4;
  
  # cut points for the discretization of the class in two intervals
  midClass = (classMin + classMax) / 2;
  cutPointsClass = c(classMin,midClass,classMax);
  
  # cut points for the discretization of the other parent in two intervals
  midParent = (parentMin + parentMax) / 2;
  cutPointsParent = c(parentMin,midParent,parentMax);
  
  # initialize the matrix for the conditional probabilities
  prob = matrix(nrow=4,ncol=length(dataVal));
  
  # do the discretization of C
  CDisc = vector(mode='numeric',length(C));
  CDisc[C>=cutPointsClass[1] & C<cutPointsClass[2]] = 1;
  CDisc[C>=cutPointsClass[2] & C<=cutPointsClass[3]] = 2;
  
  # do the discretization of the other parent
  parentDisc = vector(mode='numeric',length(C));
  parentDisc[parent>=cutPointsParent[1] & parent<cutPointsParent[2]] = 1;
  parentDisc[parent>=cutPointsParent[2] & parent<=cutPointsParent[3]] = 2;
  
  # both parents are discrete variables from now on
  # add all the values in case they are not represented in the discretization
  classVal = 1:(length(cutPointsClass)-1);
  parentVal = 1:(length(cutPointsParent)-1);
  
  # the function that calculates the conditional probabilities, which are matrices
  # with the info by rows (row 1...n for 1...n in matrix discr)
  p = Adjust2ParentsDDD(data,dataVal,parentDisc,parentVal,CDisc,classVal,verbose,plots);
  
  prob[1,] = p[[1]][,1,1];               prob[2,] = p[[1]][,1,2];
  prob[3,] = p[[1]][,2,1];               prob[4,] = p[[1]][,2,2];
  
  # calculate the BIC of the first discretization
  bic = BIC2ParentsCC(prob,data,dataVal,C,parent,discr,varType=0,verbose);
  
  # initialize the decision of dividing the class
  divide = TRUE;
  
  # continue subdividing
  while (divide==TRUE)
  {
    # different values of the division (ordered as they appear in the vector)
    indices = as.numeric(levels(factor(discr)));
    
    # number of intervals
    nIntervals = max(indices);
    
    # create a vector for the bic, cut points, level and adjust of each division
    newBicP = vector('numeric',nIntervals);              newBicC = newBicP;
    newProbP = vector('list',nIntervals);                newProbC = newProbP;
    
    # create an array with discr in each row to modify in each subdivision
    newDiscrP = replicate(nIntervals,discr);               newDiscrC = newDiscrP;
    
    # for each interval of the discretization (horizontal divisions)
    for (j in 1:nIntervals)
    {
      # add all the values in case they are not represented in the discretization
      newClassVal = 1:(nIntervals+1);
      
      # select the indices to divide (the second half of the indices to divide
      # is changed by a new number, which is max+1)
      ind = PartMatrix(discr,indices[j]);
      
      # when the length of the interval is 1, don't continue with the discretization
      if (length(ind[[1]]) > 1)
      {
        # create the new discretization including the new number
        x = ind[[1]];
        rows = ceiling((x[1] + x[length(x)])/2) : x[length(x)];
        newDiscrP[rows,ind[[2]],j] = max(discr) + 1;
        
        #cat('\n nueva discretizacion horizontal \n');   ShowMatrix(newDiscrP[,,j]);
        
        # different values of the new division (ordered as they appear in the vector)
        newIndices = as.numeric(levels(factor(newDiscrP[,,j])));
        
        # discretize the parents again and calculate the number of observations
        parentDisc = DiscretizeParent(newDiscrP[,,j],parent,C);
        nData = as.data.frame(table(parentDisc))[,2];
        
        # when there are less than 5% of data in one interval, that division is dismissed
        if (max(nData)<(length(parent)/20))   {   newBicP[j] = NA;   }
        else
        {
          # create the matrix for the adjust dividing the current interval
          probj = matrix(nrow=nIntervals+1,ncol=length(dataVal));
          
          # do the adjust for each interval in the new division
          for (k in 1:(nIntervals+1))
          {
            # data to adjust
            d = data[parentDisc==k];
            r = rep(0,length(d));
            
            # the function that calculates the conditional probabilities and the values
            aux = Adjust1ParentDD(d,dataVal,r,0,verbose,plots);
            probj[k,] = aux[[1]];
          }
          
          # calculate the BIC of the new discretization
          newBicP[j] = BIC2ParentsCC(probj,data,dataVal,C,parent,newDiscrP[,,j],varType=0,
                                   verbose);
          
          # fill in the list with the conditional probabilities
          newProbP[[j]] = probj;
        }
      }
      
      else     {   newBicP[j] = NA;   }
    }
    
    
    #cat('\n \n \n TERMINA LAS PARTICIONES HORIZONTALES');
    
    # for each interval of the discretization (vertical divisions)
    for (j in 1:nIntervals)
    {
      # add all the values in case they are not represented in the discretization
      newClassVal = 1:(nIntervals+1);
      
      # select the indices to divide (the second half of the indices to divide
      # is changed by a new number, which is max+1)
      ind = PartMatrix(discr,indices[j]);
      
      # when the length of the interval is 1, don't continue with the discretization
      if (length(ind[[2]]) > 1)
      {
        # create the new discretization including the new number
        y = ind[[2]];
        columns = ceiling((y[1] + y[length(y)])/2) : y[length(y)];
        newDiscrC[ind[[1]],columns,j] = max(discr) + 1;
        
        #cat('\n nueva discretizacion vertical \n');   ShowMatrix(newDiscrC[,,j]);
        
        # different values of the new division (ordered as they appear in the vector)
        newIndices = as.numeric(levels(factor(newDiscrC[,,j])));
        
        # discretize the parents again and calculate the number of observations
        parentDisc = DiscretizeParent(newDiscrC[,,j],parent,C);
        nData = as.data.frame(table(parentDisc))[,2];
        
        # when there are less than 5% of data in one interval, that division is dismissed
        if (max(nData)<(length(parent)/20))   {   newBicC[j] = NA;   }
        else
        {
          # create the matrix for the adjust dividing the current interval
          probj = matrix(nrow=nIntervals+1,ncol=length(dataVal));
          
          # do the adjust for each interval in the new division
          for (k in 1:(nIntervals+1))
          {
            # data to adjust
            d = data[parentDisc==k];
            r = rep(0,length(d));
            
            # the function that calculates the conditional probabilities and the values
            aux = Adjust1ParentDD(d,dataVal,r,0,verbose,plots);
            probj[k,] = aux[[1]];
          }
          
          # calculate the BIC of the new discretization
          newBicC[j] = BIC2ParentsCC(probj,data,dataVal,C,parent,newDiscrC[,,j],varType=0,
                                   verbose);
          
          # fill in the list with the conditional probabilities
          newProbC[[j]] = probj;
        }
      }
      
      else     {   newBicC[j] = NA;   }
    }
    
    #cat('\n \n \n TERMINA LAS PARTICIONES VERTICALES');
    
    #cat('\n bic antiguo = ',bic);
    #cat('\n bic horizontales = ',newBicP);
    #cat('\n bic verticales = ',newBicC);
    
    # select the best division among the class and the parent
    # there can be NAs when the interval has length 1 and is not divided
    bestP = which(newBicP==max(newBicP,na.rm=TRUE))[1];
    bestC = which(newBicC==max(newBicC,na.rm=TRUE))[1];
    
    # decide which division is the best one
    if (newBicP[bestP]>newBicC[bestC])
    {
      # check if the division increases the bic
      if (newBicP[bestP]<=bic)   {   divide = FALSE;   }
      else
      {
        # select the new discretization
        discr = newDiscrP[,,bestP];
        
        #cat('\n Mejor discretizacion horizontal = ');    ShowMatrix(discr);
        # select the best bic
        bic = newBicP[bestP];
        
        # select the best conditional probabilities
        prob = newProbP[[bestP]];
        
        # calculate the level of the new discretization in each interval
        levelDiscr = matrix(nrow=nIntervals+1,ncol=2);
        for (k in 1:(nIntervals+1))
        {
          # calculate the level of the discretization and the beginning of the interval
          part = PartMatrix(discr,k);
          levelDiscr[k,1] = length(part[[1]]);        levelDiscr[k,2] = length(part[[2]]);
        }
        
        # when the size of one part is 1x1, stop divisions
        if (min(apply(levelDiscr,1,prod))==1)   {   divide = FALSE;   }
      }
    }
    else
    {
      # check if the division increases the bic
      if (newBicC[bestC]<=bic)   {   divide = FALSE;   }
      else
      {
        # select the new discretization
        discr = newDiscrC[,,bestC];
        
        #cat('\n Mejor discretizacion vertical = ');    ShowMatrix(discr);
        # select the best bic
        bic = newBicC[bestC];
        
        # select the best conditional probabilities
        prob = newProbC[[bestC]];
        
        # calculate the level of the new discretization in each interval
        levelDiscr = matrix(nrow=nIntervals+1,ncol=2);
        for (k in 1:(nIntervals+1))
        {
          # calculate the level of the discretization and the beginning of the interval
          part = PartMatrix(discr,k);
          levelDiscr[k,1] = length(part[[1]]);      levelDiscr[k,2] = length(part[[2]]);
        }
        
        # when the size of one part is 1x1, stop divisions
        if (min(apply(levelDiscr,1,prod))==1)   {   divide = FALSE;   }
      }
    }
    
    #cat('\n nueva discretizacion = ');   ShowMatrix(discr);
  }
  
  # assign names to the matrix of probabilities
  rownames(prob) = paste('int',1:dim(prob)[1],sep='_');
  colnames(prob) = paste('X',dataVal,sep='_');
  
  return(list(prob,discr,bic));
}
############################################################









Adjust2ParentsCCC = function(data,dataMin,dataMax,parent,parentMin,parentMax,C,classMin,
                             classMax,n,maxDegree,verbose=FALSE,plots=FALSE)
{
  # It calculates the probabilities of a continuous variable with two continuous parents
  # INPUT: - data = the observations of the continuous variable
  #        - dataMin, dataMax = the minimum and maximum value in the variable
  #        - parent = the observation of the continuous parent
  #        - parentMin, parentMax = the minimum and maximum value in the parent
  #        - C = the observations of the discrete class
  #        - classMin, classMax = the minimum and maximum value in the class
  #        - n = position of the variable in the dataframe
  #        - maxDegree = maximum degree for the polynomial fit
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - values = matrix with 2*P and 3*C.val columns. It contains the
  #         min and max root, the min and max value and the probab. of the queues
  #         - polynomial = matrix with the coefficients of the polynomial for 
  #         each value of the discretization of both parents
  #         - cutPoints = matrix with the cut points of each discretization
  #         - bic = BIC of the final distribution
  
  
  # initialize the matrix with the indices of both discretizations
  discr = matrix(nrow=8,ncol=8);
  discr[1:4,1:4] = 1;   discr[1:4,5:8] = 2;
  discr[5:8,1:4] = 3;   discr[5:8,5:8] = 4;
  
  cat('\n discretización al principio: ')
  ShowMatrix(discr)
  
  
  
  # cut points for the discretization of the class in two intervals
  midClass = (classMin + classMax) / 2;
  cutPointsClass = c(classMin,midClass,classMax);
  
  # cut points for the discretization of the other parent in two intervals
  midParent = (parentMin + parentMax) / 2;
  cutPointsParent = c(parentMin,midParent,parentMax);
  
  # initialize the matrices for the polynomials of each step
  F22 = array(dim=c(2, 2, maxDegree+1));        F24 = array(dim=c(2, 4, maxDegree+1));
  F42 = array(dim=c(4, 2, maxDegree+1));        F44 = array(dim=c(4, 4, maxDegree+1));
  F28 = array(dim=c(2, 8, maxDegree+1));        F82 = array(dim=c(8, 2, maxDegree+1));
  F48 = array(dim=c(4, 8, maxDegree+1));        F84 = array(dim=c(8, 4, maxDegree+1));
  F88 = array(dim=c(8, 8, maxDegree+1));
  
  # initialize the matrices for the roots and prob of the queues of each step
  A22 = matrix(nrow=2*2,ncol=2*3);         A24 = matrix(nrow=2*2,ncol=4*3);
  A42 = matrix(nrow=4*2,ncol=2*3);         A44 = matrix(nrow=4*2,ncol=4*3);
  A28 = matrix(nrow=2*2,ncol=8*3);         A82 = matrix(nrow=8*2,ncol=2*3);
  A48 = matrix(nrow=4*2,ncol=8*3);         A84 = matrix(nrow=8*2,ncol=4*3);
  A88 = matrix(nrow=8*2,ncol=8*3);
  
  # initialize the matrix for polynomials and values
  poly = matrix(nrow=4,ncol=maxDegree+1);
  values = matrix(nrow=2,ncol=3*4);
  
  # do the discretization of C
  CDisc = vector(mode='numeric',length(C));
  CDisc[C>=cutPointsClass[1] & C<cutPointsClass[2]] = 1;
  CDisc[C>=cutPointsClass[2] & C<=cutPointsClass[3]] = 2;
  
  # do the discretization of the other parent
  parentDisc = vector(mode='numeric',length(C));
  parentDisc[parent>=cutPointsParent[1] & parent<cutPointsParent[2]] = 1;
  parentDisc[parent>=cutPointsParent[2] & parent<=cutPointsParent[3]] = 2;
  
  # both parents are discrete variables from now on
  # add all the values in case they are not represented in the discretization
  classVal = 1:(length(cutPointsClass)-1);
  parentVal = 1:(length(cutPointsParent)-1);
  
  # the function that calculates the conditional probabilities and the values, which
  # are matrices with the info by rows (row 1...n for 1...n in matrix discr)
  prob = Adjust2ParentsCDD(data,dataMin,dataMax,parentDisc,parentVal,CDisc,classVal,n,
                           maxDegree,verbose,plots);
  poly[1,] = prob[[2]][1,,1];              poly[2,] = prob[[2]][1,,2];
  poly[3,] = prob[[2]][2,,1];              poly[4,] = prob[[2]][2,,2];
  values[,1:3] = prob[[1]][1:2,1:3];             values[,4:6] = prob[[1]][1:2,4:6];
  values[,7:9] = prob[[1]][3:4,1:3];             values[,10:12] = prob[[1]][3:4,4:6];
  
  # fill in the first level matrices
  F22 = aperm(prob[[2]],c(1,3,2));            A22 = prob[[1]];
  
  # calculate the BIC of the first discretization
  bic = BIC2ParentsCC(poly,data,values,C,parent,discr,varType=1,verbose);
  
  # initialize the decision of dividing the class
  divide = TRUE;
  
  # continue subdividing
  while (divide==TRUE)
  {
    # different values of the division (ordered as they appear in the vector)
    indices = as.numeric(levels(factor(discr)));
    
    # number of intervals
    nIntervals = max(indices);
    
    # create a vector for the bic, cut points, level and adjust of each division
    newBicP = vector('numeric',nIntervals);              newBicC = newBicP;
    newPolyP = vector('list',nIntervals);                newPolyC = newPolyP;
    newValuesP = vector('list',nIntervals);              newValuesC = newValuesP;
    
    # create an array with discr in each row to modify in each subdivision
    newDiscrP = replicate(nIntervals,discr);               newDiscrC = newDiscrP;
    
    # for each interval of the discretization (horizontal divisions)
    for (j in 1:nIntervals)
    {
      # add all the values in case they are not represented in the discretization
      newClassVal = 1:(nIntervals+1);
      
      # select the indices to divide (the second half of the indices to divide
      # is changed by a new number, which is max+1)
      ind = PartMatrix(discr,indices[j]);
      
      # when the length of the interval is 1, don't continue with the discretization
      if (length(ind[[1]]) > 1)
      {
        # create the new discretization including the new number
        x = ind[[1]];
        rows = ceiling((x[1] + x[length(x)])/2) : x[length(x)];
        newDiscrP[rows,ind[[2]],j] = max(discr) + 1;
        
        cat('\n nueva discretizacion horizontal \n');   ShowMatrix(newDiscrP[,,j]);
        
        # different values of the new division (ordered as they appear in the vector)
        newIndices = as.numeric(levels(factor(newDiscrP[,,j])));
        
        # discretize the parents again and calculate the number of observations
        parentDisc = DiscretizeParent(newDiscrP[,,j],parent,C);
        nData = as.data.frame(table(parentDisc))[,2];
        
        # when there are less than 5% of data in one interval, that division is dismissed
        if (max(nData)<(length(parent)/20))   {   newBicP[j] = NA;   }
        else
        {
          # create the matrices for the adjust dividing the current interval
          polyj = matrix(nrow=nIntervals+1,ncol=maxDegree+1);
          valuesj = matrix(nrow=2,ncol=3*(nIntervals+1));
          
          # do the adjust for each interval in the new division
          for (k in 1:(nIntervals+1))
          {
            #cat('\n \n Intervalo horizontal  ',k);
            # calculate the level of the discretization and the beginning of the interval
            levelDiscr = PartMatrix(newDiscrP[,,j],newIndices[k]);
            levelP = length(levelDiscr[[1]]);        levelC = length(levelDiscr[[2]]);
            
            # first index of the submatrix
            startP = levelDiscr[[1]][1];             startC = levelDiscr[[2]][1];
            
            # data to fit
            d = data[parentDisc==k];
            
            # check if a previous adjust can be used again
            if (levelP==4 && levelC==4)             # 2x2 matrix (blocks 4x4)
            {
              # change the matrices according to the previous info or the new adjusts
              sol2x2 = Cases2x2(F22,A22,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,
                                startC,d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol2x2[[1]];         valuesj[,(3*k-2):(3*k)] = sol2x2[[2]];
              #cat('\n cambia un 2x2');
            }
            
            else if (levelP==4 && levelC==2)             # 2x4 matrix (blocks 4x2)
            {
              # change the matrices according to the previous info or the new adjusts
              sol2x4 = Cases2x4(F24,A24,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,startC,
                                d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol2x4[[1]];         valuesj[,(3*k-2):(3*k)] = sol2x4[[2]];
              F24 = sol2x4[[3]];               A24 = sol2x4[[4]];
              #cat('\n cambia un 2x4');
            }
            
            else if (levelP==2 && levelC==4)             # 4x2 matrix (blocks 2x4)
            {
              # change the matrices according to the previous info or the new adjusts
              sol4x2 = Cases4x2(F42,A42,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,startC,
                                d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol4x2[[1]];         valuesj[,(3*k-2):(3*k)] = sol4x2[[2]];
              F42 = sol4x2[[3]];               A42 = sol4x2[[4]];
              #cat('\n cambia un 4x2');
            }
            
            else if (levelP==2 && levelC==2)             # 4x4 matrix (blocks 2x2)
            {
              # change the matrices according to the previous info or the new adjusts
              sol4x4 = Cases4x4(F44,A44,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,startC,
                                d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol4x4[[1]];         valuesj[,(3*k-2):(3*k)] = sol4x4[[2]];
              F44 = sol4x4[[3]];               A44 = sol4x4[[4]];
              #cat('\n cambia un 4x4');
            }
            
            else if (levelP==4 && levelC==1)             # 2x8 matrix (blocks 4x1)
            {
              # change the matrices according to the previous info or the new adjusts
              sol2x8 = Cases2x8(F28,A28,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,startC,
                                d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol2x8[[1]];         valuesj[,(3*k-2):(3*k)] = sol2x8[[2]];
              F28 = sol2x8[[3]];               A28 = sol2x8[[4]];
              #cat('\n cambia un 2x8');
            }
            
            else if (levelP==1 && levelC==4)             # 8x2 matrix (blocks 1x4)
            {
              # change the matrices according to the previous info or the new adjusts
              sol8x2 = Cases8x2(F82,A82,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,startC,
                                d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol8x2[[1]];         valuesj[,(3*k-2):(3*k)] = sol8x2[[2]];
              F82 = sol8x2[[3]];               A82 = sol8x2[[4]];
              #cat('\n cambia un 8x2');
            }
            
            else if (levelP==2 && levelC==1)             # 4x8 matrix (blocks 2x1)
            {
              # change the matrices according to the previous info or the new adjusts
              sol4x8 = Cases4x8(F48,A48,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,startC,
                                d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol4x8[[1]];         valuesj[,(3*k-2):(3*k)] = sol4x8[[2]];
              F48 = sol4x8[[3]];               A48 = sol4x8[[4]];
              #cat('\n cambia un 4x8');
            }
            
            else if (levelP==1 && levelC==2)             # 8x4 matrix (blocks 1x2)
            {
              # change the matrices according to the previous info or the new adjusts
              sol8x4 = Cases8x4(F84,A84,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,startC,
                                d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol8x4[[1]];         valuesj[,(3*k-2):(3*k)] = sol8x4[[2]];
              F84 = sol8x4[[3]];               A84 = sol8x4[[4]];
              #cat('\n cambia un 8x4');
            }
            
            else if (levelP==1 && levelC==1)             # 8x8 matrix (blocks 1x1)
            {
              # change the matrices according to the previous info or the new adjusts
              sol8x8 = Cases8x8(F88,A88,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,startC,
                                d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol8x8[[1]];         valuesj[,(3*k-2):(3*k)] = sol8x8[[2]];
              F88 = sol8x8[[3]];               A88 = sol8x8[[4]];
              #cat('\n cambia un 8x8');
            }
            
          }
          
          # calculate the BIC of the new discretization
          newBicP[j] = BIC2ParentsCC(polyj,data,valuesj,C,parent,newDiscrP[,,j],varType=1,
                                     verbose);
          
          # fill in the list with polynomials and values of different discretizations
          newPolyP[[j]] = polyj;     newValuesP[[j]] = valuesj;
          
        }
      }
      else     {   newBicP[j] = NA;   }
    }
    
    cat('\n \n \n TERMINA LAS PARTICIONES HORIZONTALES');
    
    # for each interval of the discretization (vertical divisions)
    for (j in 1:nIntervals)
    {
      # add all the values in case they are not represented in the discretization
      newClassVal = 1:(nIntervals+1);
      
      # select the indices to divide (the second half of the indices to divide
      # is changed by a new number, which is max+1)
      ind = PartMatrix(discr,indices[j]);
      
      # when the length of the interval is 1, don't continue with the discretization
      if (length(ind[[2]]) > 1)
      {
        # create the new discretization including the new number
        y = ind[[2]];
        columns = ceiling((y[1] + y[length(y)])/2) : y[length(y)];
        newDiscrC[ind[[1]],columns,j] = max(discr) + 1;
        
        cat('\n nueva discretizacion vertical \n');   ShowMatrix(newDiscrC[,,j]);
        
        # different values of the new division (ordered as they appear in the vector)
        newIndices = as.numeric(levels(factor(newDiscrC[,,j])));
        
        # discretize the parents again and calculate the number of observations
        parentDisc = DiscretizeParent(newDiscrC[,,j],parent,C);
        nData = as.data.frame(table(parentDisc))[,2]
        
        # when there are less than 5% of data in one interval, that division is dismissed
        if (max(nData)<(length(parent)/20))   {   newBicC[j] = NA;   }
        else
        {
          # create the matrices for the adjust dividing the current interval
          polyj = matrix(nrow=nIntervals+1,ncol=maxDegree+1);
          valuesj = matrix(nrow=2,ncol=3*(nIntervals+1));
          
          # do the adjust for each interval in the new division
          for (k in 1:(nIntervals+1))
          {
            # calculate the level of the discretization and the beginning of the interval
            levelDiscr = PartMatrix(newDiscrC[,,j],newIndices[k]);
            levelP = length(levelDiscr[[1]]);        levelC = length(levelDiscr[[2]]);
            
            # first index of the submatrix
            startP = levelDiscr[[1]][1];             startC = levelDiscr[[2]][1];
            
            # data to fit
            d = data[parentDisc==k];
            
            # check if a previous adjust can be used again
            if (levelP==4 && levelC==4)             # 2x2 matrix (blocks 4x4)
            {
              # change the matrices according to the previous info or the new adjusts
              sol2x2 = Cases2x2(F22,A22,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,
                                startC,d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol2x2[[1]];         valuesj[,(3*k-2):(3*k)] = sol2x2[[2]];
              #cat('\n cambia un 2x2');
            }
            
            else if (levelP==4 && levelC==2)             # 2x4 matrix (blocks 4x2)
            {
              # change the matrices according to the previous info or the new adjusts
              sol2x4 = Cases2x4(F24,A24,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,startC,
                                d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol2x4[[1]];         valuesj[,(3*k-2):(3*k)] = sol2x4[[2]];
              F24 = sol2x4[[3]];               A24 = sol2x4[[4]];
              #cat('\n cambia un 2x4');
            }
            
            else if (levelP==2 && levelC==4)             # 4x2 matrix (blocks 2x4)
            {
              # change the matrices according to the previous info or the new adjusts
              sol4x2 = Cases4x2(F42,A42,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,startC,
                                d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol4x2[[1]];         valuesj[,(3*k-2):(3*k)] = sol4x2[[2]];
              F42 = sol4x2[[3]];               A42 = sol4x2[[4]];
              #cat('\n cambia un 4x2');
            }
            
            else if (levelP==2 && levelC==2)             # 4x4 matrix (blocks 2x2)
            {
              # change the matrices according to the previous info or the new adjusts
              sol4x4 = Cases4x4(F44,A44,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,startC,
                                d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol4x4[[1]];         valuesj[,(3*k-2):(3*k)] = sol4x4[[2]];
              F44 = sol4x4[[3]];               A44 = sol4x4[[4]];
              #cat('\n cambia un 4x4');
            }
            
            else if (levelP==4 && levelC==1)             # 2x8 matrix (blocks 4x1)
            {
              # change the matrices according to the previous info or the new adjusts
              sol2x8 = Cases2x8(F28,A28,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,startC,
                                d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol2x8[[1]];         valuesj[,(3*k-2):(3*k)] = sol2x8[[2]];
              F28 = sol2x8[[3]];               A28 = sol2x8[[4]];
              #cat('\n cambia un 2x8');
            }
            
            else if (levelP==1 && levelC==4)             # 8x2 matrix (blocks 1x4)
            {
              # change the matrices according to the previous info or the new adjusts
              sol8x2 = Cases8x2(F82,A82,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,startC,
                                d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol8x2[[1]];         valuesj[,(3*k-2):(3*k)] = sol8x2[[2]];
              F82 = sol8x2[[3]];               A82 = sol8x2[[4]];
              #cat('\n cambia un 8x2');
            }
            
            else if (levelP==2 && levelC==1)             # 4x8 matrix (blocks 2x1)
            {
              # change the matrices according to the previous info or the new adjusts
              sol4x8 = Cases4x8(F48,A48,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,startC,
                                d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol4x8[[1]];         valuesj[,(3*k-2):(3*k)] = sol4x8[[2]];
              F48 = sol4x8[[3]];               A48 = sol4x8[[4]];
              #cat('\n cambia un 4x8');
            }
            
            else if (levelP==1 && levelC==2)             # 8x4 matrix (blocks 1x2)
            {
              # change the matrices according to the previous info or the new adjusts
              sol8x4 = Cases8x4(F84,A84,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,startC,
                                d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol8x4[[1]];         valuesj[,(3*k-2):(3*k)] = sol8x4[[2]];
              F84 = sol8x4[[3]];               A84 = sol8x4[[4]];
              #cat('\n cambia un 8x4');
            }
            
            else if (levelP==1 && levelC==1)             # 8x8 matrix (blocks 1x1)
            {
              # change the matrices according to the previous info or the new adjusts
              sol8x8 = Cases8x8(F88,A88,polyj[k,],valuesj[,(3*k-2):(3*k)],startP,startC,
                                d,dataMin,dataMax,n,k,maxDegree,verbose,plots);
              
              # update the matrices
              polyj[k,] = sol8x8[[1]];         valuesj[,(3*k-2):(3*k)] = sol8x8[[2]];
              F88 = sol8x8[[3]];               A88 = sol8x8[[4]];
              #cat('\n cambia un 8x8');
            }
          }
          
          # calculate the BIC of the new discretization
          newBicC[j] = BIC2ParentsCC(polyj,data,valuesj,C,parent,newDiscrC[,,j],varType=1,
                                     verbose);
          
          # fill in the list with polynomials and values of different discretizations
          newPolyC[[j]] = polyj;     newValuesC[[j]] = valuesj;
        }
      }
      
      else     {   newBicC[j] = NA;   }
      
    }
    
    cat('\n \n \n TERMINA LAS PARTICIONES VERTICALES');
    
    cat('\n bic antiguo = ',bic);
    cat('\n bic horizontales = ',newBicP);
    cat('\n bic verticales = ',newBicC);
    
    # select the best division among the class and the parent
    # there can be NAs when the interval has length 1 and is not divided
    bestP = which(newBicP==max(newBicP,na.rm=TRUE))[1];
    bestC = which(newBicC==max(newBicC,na.rm=TRUE))[1];
    
    # choose the best division
    if (newBicP[bestP]>newBicC[bestC])
    {
      # check if the division increases the bic
      if (newBicP[bestP]<=bic)   {   divide = FALSE;   }
      else
      {
        # select the new discretization
        discr = newDiscrP[,,bestP];
        
        cat('\n Mejor discretizacion horizontal = ');    ShowMatrix(discr);
        
        
        
        
        # select the best bic
        bic = newBicP[bestP];
        
        # select the best polynomial and values
        poly = newPolyP[[bestP]];     values = newValuesP[[bestP]];
        
        # calculate the level of the new discretization in each interval
        levelDiscr = matrix(nrow=nIntervals+1,ncol=2);
        for (k in 1:(nIntervals+1))
        {
          # calculate the level of the discretization and the beginning of the interval
          part = PartMatrix(discr,k);
          levelDiscr[k,1] = length(part[[1]]);        levelDiscr[k,2] = length(part[[2]]);
        }
        
        # when the size of one part is 1x1, stop divisions
        if (min(apply(levelDiscr,1,prod))==1)   {   divide = FALSE;   }
      }
    }
    else
    {
      # check if the division increases the bic
      if (newBicC[bestC]<=bic)   {   divide = FALSE;   }
      else
      {
        # select the new discretization
        discr = newDiscrC[,,bestC];
        
        cat('\n Mejor discretizacion vertical = ');    ShowMatrix(discr);
        
        
        
        # select the best bic
        bic = newBicC[bestC];
        
        # select the best polynomial and values
        poly = newPolyC[[bestC]];     values = newValuesC[[bestC]];
        
        # calculate the level of the new discretization in each interval
        levelDiscr = matrix(nrow=nIntervals+1,ncol=2);
        for (k in 1:(nIntervals+1))
        {
          # calculate the level of the discretization and the beginning of the interval
          part = PartMatrix(discr,k);
          levelDiscr[k,1] = length(part[[1]]);        levelDiscr[k,2] = length(part[[2]]);
        }
        
        # when the size of one part is 1x1, stop divisions
        if (min(apply(levelDiscr,1,prod))==1)   {   divide = FALSE;   }
      }
    }
    
    cat('\n nueva discretizacion = ');   ShowMatrix(discr);
  }
  
  # assign names to the three dimensions of values
  r = rep(c('roots','extremes','probQueue'),dim(poly)[1]);
  colnames(values) = paste(r,rep(1:dim(poly)[1],each=3),sep='_');
  rownames(values) = c('lower','upper');
                           
  # assign names to the matrix of polynomials
  rownames(poly) = paste('int',1:dim(poly)[1],sep='_');
  colnames(poly) = paste('x',0:maxDegree,sep='^');
  
  return(list(values,poly,discr,bic));
}
############################################################









