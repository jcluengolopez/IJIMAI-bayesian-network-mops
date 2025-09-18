
# CARGAR DEPENDENCIAS AL PRINCIPIO (discretizaci?n)
library(discretization);
library(classInt);


DiscOrCont = function(data,discretize)
{
  # It decides if each variable of the data is discrete or continuous
  # and returns a vector composed by 0 = discr, 1 = cont or 2 = discretized
  # INPUT:  - data = a data-frame with the data of study
  #         - discretize = it decides the type of discretization. 
  #       0=never, 1=between 20 and 5% of the values, 2=more than 20 values
  # OUTPUT:  - x = a vector with 0 (discrete), 1 (continuous) or 2 (discretized)
  
  
  # number of variables
  nVar = ncol(data);
  
  # change to 0 if the variable is a character and numeric variables depend on
  # the value of the discretize parameter
  if (discretize==0)   # no discretization
  {
    # create a 1-vector (the class is discretizated in regression)
    x = rep(1,nVar);
    for (i in 1:nVar)
    {
      if (is.character(data[,i]) | length(table(data[,i]))<=20)   { x[i] = 0; }
    }
  }
  else if (discretize==1)   # discretize pseudo-discrete variables
  {
    # create a 1-vector (the class is discretized in regression)
    x = rep(1,nVar);
    
    for (i in 1:nVar)
    {
      # check if it is a dscrete variable
      if (is.character(data[,i]) | length(table(data[,i]))<=20)   {   x[i] = 0;   }
      
      # when there are less values than 5% of the data 
      else if (length(table(data[,i]))<=0.05*nrow(data))   {   x[i] = 2;   }
    }
  }
  else if (discretize==2)   # discretize every non-discrete variable
  {
    # create a 2-vector
    x = rep(2,nVar);
    
    for (i in 1:nVar)
    {
      # check if it is a dscrete variable
      if (is.character(data[,i]) | length(table(data[,i]))<=20)   { x[i] = 0; }
    }
  }
  
  # return the vector with 0s, 1s and 2s
  return(x);
}
############################################################






Discretization = function(data,varType,discretizeType,dataMin,dataMax,model)
{
  # it does the discretization of the variables whose varType is 2 and substitute 
  # the discretized data in the dataframe. It makes two different discretizations
  # INPUT:  - data = dataframe with the data to discretize
  #         - varType = a vector with 0 if the variable is discrete, 1
  #       if it's continuous and 2 if it has to be discretized. Program calculates it
  #        - discretizeType = it decides the type of discretization. 0=dynamic, 1=equal width
  #        - dataMin,dataMax = minimum and maximum value of each variable
  #        - model (optional) = 'regression' or 'classification'
  # OUTPUT:  - data = dataframe with the discretized data
  #          - splitPoints = the cut points of the discretization (it can be NULL)
  
  
  # number of variables
  nVar = ncol(data);
  
  # if there's not model, calculate it
  if (missing(model) & varType[nVar]==0)   {   model = 'classification';   }
  else if (missing(model) & varType[nVar]>0)   {   model = 'regression';   }
  
  # if there're not dataMin and dataMax, calculate them
  if (missing(dataMin))   {   dataMin = as.numeric(sapply(data,min)) - 0.000000001;   }
  if (missing(dataMax))   {   dataMax = as.numeric(sapply(data,max)) + 0.000000001;   }
  
  # index of the discretizable variables including the class
  indexDisc = which(varType==2);
  nDisc = length(indexDisc);
  
  # index of the discretizable variables not including the class
  indexDiscNoClass = which(varType[-nVar]==2);
  nDiscNoClass = length(indexDiscNoClass);
  
  # check if there're discretizable variables
  if (nDisc==0)      # nothing to discretize
  {
    indexDisc = 0;        cutPoints = 0;
  }
  else
  {
    if (discretizeType==0)      # dynamic discretization
    {
      # check if there are observation variables to discretize
      if (nDiscNoClass>0)     # variables to discretize (classification or regression)
      {
        # select the variables to discretize
        disc = data[,c(indexDiscNoClass,nVar)];
        
        # discretize the data set
        dataDisc = discretization :: mdlp(disc);
        
        # copy the values and the splitPoints of the discretization
        data[,indexDiscNoClass] = dataDisc$Disc.data[,1:nDiscNoClass];
        splitPoints = dataDisc$cutp;
        
        # check if the class has to be discretized
        if (varType[nVar]==2)
        {
          # calculate the cut points for the discretization of the class
          cutPoints = classInt :: classIntervals(data[,nVar],style ='equal')$brks;
          
          # change the first point
          cutPoints[1] = cutPoints[1] - 0.00001;
          
          # make the discretization of the class using the equal width method
          data[,nVar] = cut(data[,nVar],breaks=cutPoints,labels=FALSE);
          
          # quit the first and last cut points of the discretization
          if (length(cutPoints) > 2)
          {
            splitPoints[[nDisc]] = cutPoints[-c(1,length(cutPoints))];
          }
          else   {   splitPoints[[nDisc]]=='All';   }
        }
      }
      else     # no variable to discretize (only the class in regression)
      {
        # calculate the cut points for the discretization of the class
        cutPoints = classInt :: classIntervals(data[,nVar],style ='equal')$brks;
        
        # change the first point
        cutPoints[1] = cutPoints[1] - 0.00001;
        
        # make the discretization of the class using the equal width method
        data[,nVar] = cut(data[,nVar],breaks=cutPoints,labels=FALSE);
        
        # quit the first and last cut points of the discretization
        if (length(cutPoints) > 2)
        {
          splitPoints = list(cutPoints[-c(1,length(cutPoints))]);
        }
        else   {   splitPoints[[nDisc]] = list('All');   }
      }
    }
    else     # equal width discretization (for classification and regression)
    {
      # initialize the split points
      splitPoints = vector(mode='list',length=nDisc);
      
      for (i in 1:nDisc)
      {
        # calculate the cut points for the discretization
        cutPoints = classInt :: classIntervals(data[,indexDisc[i]],style ='equal')$brks;
        #cutPoints = classInt :: classIntervals(data[,indexDisc[i]],n=3,style ='equal')$brks;
        
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






DiscreteProb = function(data,dataVal,verbose=FALSE,plots=FALSE)
{
  # it calculates the probability of each value in a discrete variable
  # INPUT: - data = the observations of the discrete variable
  #        - dataVal = the different values of the discrete variable
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - dataTable = a 2 column matrix with the different values of the 
  #         variable in the first column and each probability in the second column
  
  
  # add the different values of the variable to the observations
  extendedData = append(data,dataVal);
  
  # different values of the extended data and the frequency
  dataTable = as.data.frame(table(extendedData));
  
  # subtract one extra data to each value in order to exclude the previous one added
  # apply Laplace's correction
  dataTable[,2] = dataTable[,2] / (length(data) + length(dataVal));
  
  # return the probabilities of each value of the data
  return(dataTable);
}
############################################################








PositionMatrix = function(dataX,dataY,minX,maxX,minY,maxY)
{
  # it obtains the position in a 8x8 matrix of the observation (dataX, dataY)
  # INPUT - dataX, dataY = the observations of two variables
  #       - minx, maxX = the min and max values of variable X
  #       - minY, maxY = the min and max values of variable Y
  # OUTPUT: - pos = a vector with the row and column in the matrix
  
  
  # range in both variables
  rangeX = maxX - minX;              rangeY = maxY - minY;
  
  # obtain the position
  posX = ceiling(dataX/rangeX);             posY = ceiling(dataY/rangeY);
  
  return(c(posX,posY));
}
############################################################





PartMatrix = function(matrix,val)
{
  # it obtains the indices of the submatrix filled in with the value val
  # INPUT - matrix = a matrix with 1...n values
  #       - val = a number between 1 and n
  # OUTPUT: - rows = a vector with the rows with that value
  #         - columns = a vector with the columns with that value
  
  
  # pairs of sub indices with that value
  subInd = which(matrix==val,arr.ind=TRUE);
  
  # rows that contain the value
  rows = as.numeric(levels(as.data.frame(table(subInd[,1]))[,1]));
  
  # columns that contain the value
  columns = as.numeric(levels(as.data.frame(table(subInd[,2]))[,1]));
  
  return(list(rows,columns));
}
############################################################





CutPointsDiscr = function(discr,parentMin,parentMax,classMin,classMax,val)
{
  # it create the cut points for a discretization taking in account the part of
  # matrix where val is
  # INPUT  - discr = a matrix with 1...n values
  #        - parentMin, parentMax = the minimum and maximum value in the parent
  #        - classMin, classMax = the minimum and maximum value in the class
  #        - val = a number between 1 and n
  # OUTPUT: - rows = a vector with the rows with that value
  #         - columns = a vector with the columns with that value
  
  
  # select the indices of val
  ind = PartMatrix(discr,val);
  rows = ind[[1]];              columns = ind[[2]];
  
  # beginning of each variable
  begParent = parentMin + (parentMax - parentMin) * (rows[1] - 1) / 8;
  begClass = classMin + (classMax - classMin) * (columns[1] - 1) / 8;
  
  # end of each variable
  endParent = parentMin + (parentMax - parentMin) * rows[length(rows)] / 8;
  endClass = classMin + (classMax - classMin) * columns[length(columns)] / 8;
  
  # return the beginning and end as a 2x2 matrix with parent in row1 and C in row 2
  x = c(begParent,begClass,endParent,endClass);
  
  return(matrix(x,nrow=2));
}
############################################################







DiscretizeParent = function(discr,parent,C)
{
  # it discretizes both parents (parent and class) that fits the discretization
  # made in discr
  # INPUT  - discr = a matrix with 1...n values
  #        - parent = the observation of the continuous parent
  #        - C = the observations of the discrete class
  # OUTPUT  - newDisc = vector with the combined discretization for both parents
  
  
  # the number of different values of discr will be the size of the discretization
  indices = as.numeric(levels(factor(discr)));
  nIntervals = length(indices);
  
  # initialize the discretizations
  newDisc = vector(mode='numeric',length(parent));
  
  # for each interval of the new discretization
  for (k in 1:nIntervals)
  {
    # cut points to discretize the parent and the class
    cutPoints = CutPointsDiscr(discr,min(parent),max(parent),min(C),max(C),k);
    
    # discretization of the parent (the class will be the same)
    newDisc[(parent>=cutPoints[1,1] | round(parent,6)>=round(cutPoints[1,1],6)) & 
            (parent<=cutPoints[1,2] | round(parent,6)<=round(cutPoints[1,2],6)) & 
            (C>=cutPoints[2,1] | round(C,6)>=round(cutPoints[2,1],6)) & 
            (C<=cutPoints[2,2] | round(C,6)<=round(cutPoints[2,2],6))] = k;
  }
  
  a = which(newDisc==0);
  if (length(a)>0)
  {
    cat('\n \n \n a = ',a);
    cat('\n parent(a) = ',parent[a]);
    cat('\n c(a) = ',C[a]);
    cat('\n cutPoints en DiscretizeParent = \n ');    ShowMatrix(cutPoints);
    cat('\n min(parent) = ',min(parent),'y max(parent) =',max(parent));
    cat('\n min(C) = ',min(C),'y max(C) =',max(C));
  }
  
  return(newDisc);
}
############################################################






SizeDiscr = function(discr)
{
  # it calculates the number of cut points in the discretization of two continuous
  # variables
  # INPUT:  - discr = a matrix with 1...n values
  # OUTPUT: - sizeDiscr = vector with the number of cut points in each variable
  
  
  # create the matrices for the differences
  diffRow = diff(discr);            diffCol = t(diff(t(discr)));
  
  # add matrices by rows and by columns to calculate the number of changes
  sum1 = abs(apply(diffRow, 1, sum));              sum2 = abs(apply(diffCol, 2, sum));
  
  # extract the number of changes
  changes1 = length(sum1[sum1>0]);            changes2 = length(sum2[sum2>0]);
  
  sizeDiscr = c(changes1,changes2);
  
  return(sizeDiscr);
}
############################################################








SortVariablesMIC = function(MIC,select,verbose)
{
  # it sorts the variables according to the mutual information conditioned to the class
  # excluding the variables that are already in the model
  # INPUT: - MIC = matrix with the mutual information of each pair of variables
  #        - select = variables included in the model
  # OUTPUT: - ranking = ranking of the variables sorted by its MIC
  
  
  # number of variables
  nVar = dim(MIC)[1];
  
  # order the mutual info by decreasing and removing zeroes (not necessary)
  order = sort.int(t(MIC),decreasing=TRUE,index.return=TRUE)$ix;
  
  # select only the values with mutual information >0
  n = (nVar^2-nVar)/2;   # number of elements under the diagonal
  order = order[1:n];
  
  # get the position in the matrix of each mutual information
  nRow = ceiling(order/nVar);               nCol = order - (nRow-1) * nVar;
  pos = rbind(nRow,nCol);
  
  # remove the elements that include only one of the selected variables (not 0, not 2)
  pos2 = pos[,apply(pos,2,function(x) length(x[x %in% select])) == 1];
  
  # select the variables related to the ones included in select
  var = apply(pos2,2,function(x) x[which(!x %in% select)]);
  
  # remove the repeated variables to get the ranking of variables
  ranking = unique(var);
  
  return(ranking)
}
############################################################





DiscretizationEqWid = function(data,nInt)
{
  # it makes an equal width discretization of a continuous dataset
  # INPUT:  - data = dataframe with the data to discretize
  #         - nInt = number of intervals to discretiza
  # OUTPUT:  - dataDisc = dataframe with the discretized data
  #          - breaks = the cut points of the discretization
  
  
  nVar = ncol(data);
  
  breaks = matrix(nrow=nVar,ncol=nInt+1)
  dataDisc = data
  
  # for each variable
  for (i in 1:nVar)
  {
    # calculate the breaks of the discretization
    breaks[i,] = classIntervals(data[,i],nInt,style='equal')$brks
    breaks[i,1] = breaks[i,1] - 0.00001;
    
    # discretize using the breaks
    dataDisc[,i] = cut(data[,i],breaks[i,],labels=FALSE)
  }
  
  return(list(dataDisc,breaks))
}
############################################################






