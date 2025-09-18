

PredictDataClassificationEF = function(prob,data,verbose=FALSE,plots=FALSE)
{
  # it predicts the discrete class of one single observation data using the model (efficient)
  # INPUT: - prob = the TAN model
  #        - data = set of observations without the class
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - prediction = the selected value of the class for each observation
  
  
  # change the structure of the model
  #prob = StructureModelIverse(prob2);
  
  # type and number of variables
  varType = prob[[1]];               nVar = length(varType) - 1;
  
  # matrix with the arcs of the tree
  links = prob[[2]][-(nVar+1),-(nVar+1)];
  
  # sum the columns of links to get the number of parents of each variable
  nParents = apply(links,2,sum);
  
  # distribution of the class
  classVal = prob[[length(prob)-1]];
  nClass = length(classVal);
  classProb = prob[[length(prob)]];
  
  # number of observations
  nObservations = nrow(data);
  
  
  #cat('\n classVal = ',classVal);
  #cat('\n nParents = ',nParents);
  #cat('\n nObservations = ',nObservations);
  
  # function with the probability of the class for one observation
  probC = function(x)   {   classProb[which(classVal==x)];   }
  
  
  # function with the probability of the variable with one parent for one observation
  probXC = function(x,class,posX)
  {
    # INPUT: - x = observations of the concrete variable
    #        - class = one of the factors of the class
    #        - posX = position of the variable X in the dataframe
    
    
    # the probabilities depends on the type of variable
    if (varType[posX]==1)     # continuous
    {
      # extract the roots, extremes and probabilities of the queues from the matrix values
      values = prob[[3*posX]][,(3*class-2):(3*class)];
      roots = values[,1];        extremes = values[,2];        probQueues = values[,3];
      
      # extract the coefficients that are not NA (degree lower than max)
      p = prob[[3*posX+1]][class,];
      poly = p[!is.na(p)];
      
      # predict the polynomial for the exact value
      polyVal = predict(as.polynomial(poly),x);
      
      # probability in each part
      probability1 = (x>=extremes[1] & x<roots[1]) * probQueues[1];
      probability2 = (x>=roots[1] & x<=roots[2]) * polyVal;
      probability3 = (x>roots[2] & x<=extremes[2]) * probQueues[2];
      
      # the sum of the three probabilities
      density = probability1 + probability2 + probability3;
    }
    else     # discrete or discretized
    {
      # the values of the variable depends on its type
      if (varType[posX]==0)   {   dataVal = prob[[3*posX]];   }
      else          # discretized
      {
        # the used part of the model
        splitPoints = prob[[3*posX]];    # split points of the variable
        dataVal = 1:(length(splitPoints)-1);    # the rows are 1, 2, 3, 4...
      }
      
      # extract its conditional probabilities
      dataMatrix = prob[[3*posX+1]];
      
      # obtain the position in the matrix of probabilities
      density = vapply(x, function(i) {
        if (i %in% dataVal)   {   density = dataMatrix[dataVal == i,class];   }
        else   {   density = 0;   }
        return(density)   },  1);
    }
    
    return(density);
  }
  
  
  
  
  # function with the probability of the variable with two parents for one observation
  probXYC = function(x,parent,class,posX,posParent)
  {
    # INPUT: - x = observations of the concrete variable
    #        - parent = observations of the parent
    #        - class = one of the factors of the class
    #        - posX, posParent = position of the variable X and its parent in the dataframe
    
    
    # the model depends on the type of parent
    if (varType[posParent]==1)   # 2 parents (the class and another continuous)
    {
      # the model depends on the type of variable
      if (varType[posX]==1)  # when the variable is continuous
      {
        # the used part of the model (with the rows for the class)
        values = prob[[3*posX]][(2*class-1):(2*class),];
        cutPoints = prob[[3*posX+1]][class,];
        polynom = prob[[3*posX+2]];
        
        # extract the split points in the discretization of the continuous parent
        splitPoints = cutPoints[!is.na(cutPoints)];
        
        # position in the discretization of the continuous parent
        interval = vapply (parent, function(i) {   which(i<=splitPoints[-1])[1]    }, 1);
        
        # extract the roots and probabilities of the queues from the matrix values
        rootsMin = values[1,3*interval-2];        rootsMax = values[2,3*interval-2];
        extremesMin = values[1,3*interval-1];       extremesMax = values[2,3*interval-1];
        probQueuesMin = values[1,3*interval];       probQueuesMax = values[2,3*interval];
        
        # predict the value of the polynomial for each observation
        polyVal = vapply(1:length(x), function(i) {
          p = polynom[interval[i],,class]
          p = p[!is.na(p)]
          return(predict(as.polynomial(p),x[i]))   },  1);
        
        # probability in each part
        probability1 = (x>=extremesMin & x<rootsMin) * probQueuesMin;
        probability2 = (x>=rootsMin & x<=rootsMax) * polyVal;
        probability3 = (x>rootsMax & x<=extremesMax) * probQueuesMax;
        
        # the sum of the three probabilities
        density = probability1 + probability2 + probability3;
        
        
        
        
        #cat('\n \n \n \n Variable ',posX,'con padre',posParent);
        #cat('\n x = ',x,'   padre = ',parent);
        #cat('\n splitPoints = ',splitPoints);
        #cat('\n interval = ',interval);
        #cat('\n roots = ',c(rootsMin,rootsMax),'y extremes = ',c(extremesMin,extremesMax));
        #cat('\n probQueues = ',c(probQueuesMin,probQueuesMax));
        #cat('\n polyVal = ',polyVal);
        #cat('\n probabilities = ',c(probability1,probability2,probability3));
        #cat('\n density = ',density);
      }
      
      
      else     # discrete or discretized
      {
        # the values of the variable depends on its type
        if (varType[posX]==0)   {   dataVal = prob[[3*posX]];   }
        else          # discretized
        {
          dataVal = 1:(length(prob[[3*posX]])-1);    # the values are 1, 2, 3, 4...
        }
        
        # the used part of the model (with the rows for the class)
        cutPoints = prob[[3*posX+1]][class,];
        dataMatrix = prob[[3*posX+2]];
        
        # extract the split points in the discretization of the continuous parent
        splitPoints = cutPoints[!is.na(cutPoints)];
        
        # position in the discretization of the continuous parent
        interval = vapply (parent, function(i) {   which(i<=splitPoints[-1])[1]    }, 1);
        
        # obtain the position in the matrix of probabilities
        density = vapply(1:length(x), function(i) {
          if (x[i] %in% dataVal)   {   density = dataMatrix[dataVal==x[i],interval[i],class];  }
          else   {   density = 0;   }
          return(density)   },  1);
      }
    }
    
    
    
    else    # 2 parents (the discrete class and another discrete or discretized)
    {
      # extract the position of the observation of the discrete parent
      if (varType[posParent]==0)   {   parentVal = prob[[3*posParent]];   }
      else          # discretized
      {
        parentVal = 1:(length(prob[[3*posParent]])-1);    # the values are 1, 2, 3, 4...
      }
      
      # position of the exact value of the parent
      interval = vapply (parent, function(i) {   which(i==parentVal)    }, 1);
      
      # the model depends on the type of variable
      if (varType[posX]==1)  # when the variable is continuous
      {
        # the used part of the model (with the rows for the class and the other parent)
        values = prob[[3*posX]][,(3*class-2):(3*class)];
        polynom = prob[[3*posX+1]];
        
        # extract the roots and probabilities of the queues from the matrix values
        rootsMin = values[2*interval-1,1];        rootsMax = values[2*interval,1];
        extremesMin = values[2*interval-1,2];        extremesMax = values[2*interval,2];
        probQueuesMin = values[2*interval-1,3];        probQueuesMax = values[2*interval,3];
        
        # predict the value of the polynomial for each observation
        polyVal = vapply(1:length(x), function(i) {
          p = polynom[interval[i],,class]
          p = p[!is.na(p)]
          return(predict(as.polynomial(p),x[i]))   },  1);
        
        # probability in each part
        probability1 = (x>=extremesMin & x<rootsMin) * probQueuesMin;
        probability2 = (x>=rootsMin & x<=rootsMax) * polyVal;
        probability3 = (x>rootsMax & x<=extremesMax) * probQueuesMax;
        
        # the sum of the three probabilities
        density = probability1 + probability2 + probability3;
      }
      
      
      else     # discrete or discretized
      {
        # the values of the variable depends on its type
        if (varType[posX]==0)   {   dataVal = prob[[3*posX]];   }
        else          # discretized
        {
          dataVal = 1:(length(prob[[3*posX]])-1);    # the values are 1, 2, 3, 4...
        }
        
        # the used part of the model (with the rows for the class)
        dataMatrix = prob[[3*posX+1]];
        
        #cat('\n dataMatrix',dim(dataMatrix));
        #cat('\n n observaciones = ',length(x));
        #cat('\n interval =',interval);
        #cat('\n class = ',class);
        #cat('\n matriz grande =  = ',dim(prob[[3*posX+1]]));
        
        
        # obtain the position in the matrix of probabilities
        density = vapply(1:length(x), function(i) {
          if (x[i] %in% dataVal)   {   density = dataMatrix[dataVal==x[i],interval[i],class];  }
          else   {   density = 0;   }
          return(density)   },  1);
      }
    }
    
    return(density);
  }
  
  
  
  
  
  # matrix for the conditional probabilities in each variable in each observation for a
  # concrete factor of the class
  probCond = matrix(NA, nrow=nObservations, ncol=nVar+1);
  
  # matrix for the conditional probabilities in each observation in each factor of the class
  probClassValues = matrix(NA, nrow=nObservations, ncol=nClass)
  
  # for each factor of the class
  for (i in 1:nClass)
  {
    # for each prediction variable
    for (j in 1:nVar)
    {
      # each column of the matrix corresponds to a variable
      if (nParents[j]==0)      # the node variable
      {
        # calculate the probabilities for all the observations
        probCond[,j] = probXC(data[,j],i,j);
      }
      else   # 2 parents: a parent apart from the discrete class
      {
        # obtain the parent of the variable
        posParent = which(links[,j]==1);
        
        # calculate the probabilities for all the observations
        probCond[,j] = probXYC(data[,j],data[,posParent],i,j,posParent);
      }
      
      # probability of the class
      probCond[,nVar+1] = probC(classVal[i]);
      
      # multiply the probabilities of all the prediction variables
      probClassValues[,i] = apply(probCond, 1, prod);
      
    }
  }
  
  # select the factor of the class with maximum probability for each observation
  predictions = classVal[apply(probClassValues, 1, which.max)];
  
  
  
  #cat('\n vector de probabilidades = ');  print(probClassValues);
  #cat('\n prediction = ',predictions);
  return(predictions);
}
############################################################











PredictDataRegressionEF = function(prob,data,verbose=FALSE,plots=FALSE)
{
  # it predicts the continuous class of a set of observations data using the model
  # INPUT: - prob2 = the TAN model
  #        - data = set of observations without the class
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - prediction = the selected value of the class for each observation
  
  
  # change the structure of the model
  #prob = StructureModelIverse(prob2);
  
  # type and number of variables
  varType = prob[[1]];              nVar = length(varType) - 1;
  
  # matrix with the arcs of the tree
  links = prob[[2]][-(nVar+1),-(nVar+1)];
  
  # sum the columns of links to get the number of parents of each variable
  nParents = apply(links,2,sum);
  
  # extract the information of the class (values and polynomial)
  classVal = prob[[length(prob)-1]];
  classPolynom = prob[[length(prob)]];
  classPolynom = classPolynom[!is.na(classPolynom)];
  
  # extract the roots and probabilities of the queues from the matrix values
  rootsClass = classVal[,1];        extremesClass = classVal[,2];
  probQueuesClass = classVal[,3];
  
  # calculate all the possible cut points in the discretization of the class
  cutPointsClass = extremesClass[1] + 0:8 * diff(extremesClass) / 8;
  nClass = length(cutPointsClass) - 1;
  
  # number of observations
  nObservations = nrow(data);
  
  
  #cat('\n classVal = ');   print(classVal);
  #cat('\n nParents = ',nParents);
  #cat('\n nObservations = ',nObservations);
  
  
  
  
  # function with the probability of the class for one observation
  ####### NO SE USA
  #probC = function(x)
  #{
  #  # predict the value of the polynomial for each observation
  #  polyVal = predict(as.polynomial(classPolynom), x);
  #  
  #  # probability in each part
  #  probability1 = (x>=extremesClass[1] & x<rootsClass[1]) * probQueuesClass[1];
  #  probability2 = (x>=rootsClass[1] & x<=rootsClass[2]) * polyVal;
  #  probability3 = (x>rootsClass[2] & x<=extremesClass[2]) * probQueuesClass[2];
  #  
  #  # the sum of the three probabilities
  #  density = probability1 + probability2 + probability3;
  #  
  #  return(density)
  #}
  
  
  
  # function with the integral of the class
  integralpC = function(x,k)
  {
    # INPUT: - x = observations of the class
    #        - k = 1 for p(C) and other for C*p(C)
    
    
    if (k==1)     # p(C)
    {
      # integral of each part
      i1 = integral(as.polynomial(probQueuesClass[1]));
      i2 = integral(as.polynomial(classPolynom));
      i3 = integral(as.polynomial(probQueuesClass[2]));
    }
    else     # C * p(C)
    {
      # integral of each part
      i1 = integral(c(0,1) * as.polynomial(probQueuesClass[1]));
      i2 = integral(c(0,1) * as.polynomial(classPolynom));
      i3 = integral(c(0,1) * as.polynomial(probQueuesClass[2]));
    }
    
    # integral of each part adjusted to the starting point
    # the first integral starts in 0
    integral1 = i1 - predict(i1,extremesClass[1]);
    
    # the second integral starts where the first one ends
    integral2 = i2 - predict(i2,rootsClass[1]) + predict(integral1,rootsClass[1]);
    
    # the third integral starts where the second one ends
    integral3 = i3 - predict(i3,rootsClass[2]) + predict(integral2,rootsClass[2]);
    
    # probability in each part
    probability1 = (x>=extremesClass[1] & x<rootsClass[1]) * predict(integral1,x);
    probability2 = (x>=rootsClass[1] & x<=rootsClass[2]) * predict(integral2,x);
    probability3 = (x>rootsClass[2] & x<=extremesClass[2]) * predict(integral3,x);
    
    integralClass = probability1 + probability2 + probability3;
    
    return(integralClass);
  }
  
  
  
  
  
  
  # function with the probability of the variable with one parent for one observation
  probXC = function(x,class,posX)
  {
    # INPUT: - x = observations of the concrete variable
    #        - class = one of the possible values of the class
    #        - posX = position of the variable X in the dataframe
    
    
    # cut points of the discretization of the class
    cutPClass = prob[[3*posX+1]];
    
    # select the used interval of the class (class is the mid point of the interval)
    intervalClass = which(class < cutPClass)[1] - 1;
    
    #cat('\n puntos de corte de la clase = ',cutPClass);
    #cat('\n class = ',class);
    #cat('\n intervalo de la clase = ',intervalClass);
    
    # the probabilities depends on the type of variable
    if (varType[posX]==1)     # continuous
    {
      # extract the roots, extremes and probabilities of the queues from the matrix values
      values = prob[[3*posX]][,(3*intervalClass-2):(3*intervalClass)];
      roots = values[,1];        extremes = values[,2];        probQueues = values[,3];
      
      # extract the coefficients that are not NA (degree lower than max)
      p = prob[[3*posX+2]][intervalClass,];
      poly = p[!is.na(p)];
      
      # predict the polynomial for the exact value
      polyVal = predict(as.polynomial(poly),x);
      
      # probability in each part
      probability1 = (x>=extremes[1] & x<roots[1]) * probQueues[1];
      probability2 = (x>=roots[1] & x<=roots[2]) * polyVal;
      probability3 = (x>roots[2] & x<=extremes[2])*probQueues[2];
      
      # the sum of the three probabilities
      density = probability1 + probability2 + probability3;
    }
    else     # discrete or discretized
    {
      # the values of the variable depends on its type
      if (varType[posX]==0)   {   dataVal = prob[[3*posX]];   }
      else          # discretized
      {
        dataVal = 1:(length(prob[[3*posX]])-1);    # the values are 1, 2, 3, 4...
      }
      
      # the used part of the model (with the columns for the class)
      dataMatrix = prob[[3*posX+2]][,intervalClass];
      
      # obtain the position in the matrix of probabilities
      density = vapply(x, function(i) {
        if (i %in% dataVal)   {   density = dataMatrix[dataVal == i];   }
        else   {   density = 0;   }
        return(density)   },  1);
    }
    
    return(density);
  }
  
  
  
  
  # function with the probability of the variable with two parents for one observation
  probXYC = function(x,parent,class,posX,posParent)
  {
    # INPUT: - x = observations of the concrete variable
    #        - parent = observations of the parent
    #        - class = one of the possible values of the class
    #        - posX, posParent = position of the variable X and its parent in the dataframe
    
    
    # the model depends on the type of parent
    if (varType[posParent]==1)   # 2 parents (the class and another continuous)
    {
      # extract the extremes of the continuous parent
      extremesParent = prob[[3*posParent]][1:2,2];
      
      # calculate all the possible cut points in the discretization of the parent
      cutPointsParent = extremesParent[1] + 0:8 * diff(extremesParent) / 8;
      
      # select the used interval of the class (class is the mid point of the interval)
      intervalClass = which(class < cutPointsClass)[1] - 1;
      
      # select the used interval of the parent
      intervalParent = vapply(parent, function(i) {
        which(i<=cutPointsParent[-1])[1]     },1);
      
      # select the position of the discretization of both parents
      indexDistr = prob[[3*posX+1]][intervalParent, intervalClass];
      
      
      
      
      #cat('\n cutPointsParent =',cutPointsParent);
      #cat('\n cutPointsClass =',cutPointsClass);
      #cat('\n clase = ',class,'y padre = ',parent);
      #cat('\n padre = ',parent);
      #cat('\n intervalo del padre en el que debería estar = \n',intervalParent);
      
      
      #cat('\n intervalClass =',intervalClass);
      #cat('\n indexDistr =',indexDistr);
      #cat('\n x =',x);
      
      
      # the model depends on the type of variable
      if (varType[posX]==1)  # when the variable is continuous
      {
        # the used part of the model (with the rows for the class and the other parent)
        values = prob[[3*posX]];
        polynom = prob[[3*posX+2]];
        
        # extract the roots and probabilities of the queues from the matrix values
        rootsMin = values[1,3*indexDistr-2];        rootsMax = values[2,3*indexDistr-2];
        extremesMin = values[1,3*indexDistr-1];     extremesMax = values[2,3*indexDistr-1];
        probQueuesMin = values[1,3*indexDistr];     probQueuesMax = values[2,3*indexDistr];
        
        # predict the value of the polynomial for each observation
        polyVal = vapply(1:length(x), function(i) {
          p = polynom[indexDistr[i],];
          p = p[!is.na(p)];
          return(predict(as.polynomial(p),x[i]))   },  1);
        
        
        #cat('\n raíces izq =',rootsMin);
        #cat('\n raíces dcha =',rootsMax);
        #cat('\n polyval =',polyVal);
        
        
        
        # probability in each part
        probability1 = (x>=extremesMin & x<rootsMin) * probQueuesMin;
        probability2 = (x>=rootsMin & x<=rootsMax) * polyVal;
        probability3 = (x>rootsMax & x<=extremesMax)*probQueuesMax;
        
        # the sum of the three probabilities
        density = probability1 + probability2 + probability3;
        
      }
      
      else     # discrete or discretized
      {
        # the values of the variable depends on its type
        if (varType[posX]==0)   {   dataVal = prob[[3*posX]];   }
        else          # discretized
        {
          dataVal = 1:(length(prob[[3*posX]])-1);    # the values are 1, 2, 3, 4...
        }
        
        # the used part of the model
        dataMatrix = prob[[3*posX+2]];
        
        # obtain the position in the matrix of probabilities
        density = vapply(1:length(x), function(i) {
          if (x[i] %in% dataVal)   {   density = dataMatrix[indexDistr[i],dataVal==x[i]];  }
          else   {   density = 0;   }
          return(density)   },  1);
        
        #cat('\n density = ',density);
      }
    }
    
    
    
    else    # 2 parents (the continuous class and another discrete or discretized)
    {
      # extract the position of the observation of the discrete parent
      if (varType[posParent]==0)   {   parentVal = prob[[3*posParent]];   }
      else          # discretized
      {
        parentVal = 1:(length(prob[[3*posParent]])-1);    # the values are 1, 2, 3, 4...
      }
      
      # position of the exact value of the parent
      intervalParent = vapply (parent, function(i) {   which(i==parentVal)    }, 1);
      
      #cat('\n variable ',posX,'con padre ',posParent);
      #cat('\n padre = ',parent,'que está en la posición = ',intervalParent);
      #cat('\n clase = ',class);
      #cat('\n prob[[3*posX+1]] = \n ');   print(prob[[3*posX+1]]);
      
      
      
      
      # interval where the class belongs depending on the interval of the parent
      intervalClass = vapply(intervalParent, function(i) {
        #cat('\n \n i = ',i);
        if (is.vector(prob[[3*posX+1]]))   {   cutPoints = prob[[3*posX+1]];   }
        else   {   cutPoints = prob[[3*posX+1]][i,];   }
        cutPoints = cutPoints[!is.na(cutPoints)];
        return(which(class<=cutPoints[-1])[1])  },  1);
      
      
      #cat('\n intervalClasss = ',intervalClass);
      
      
      
      
      # the model depends on the type of variable
      if (varType[posX]==1)  # when the variable is continuous
      {
        # the used part of the model (with the rows for the class)
        values = prob[[3*posX]];
        polynom = prob[[3*posX+2]];
        
        # the roots of values
        rootsMin = vapply(1:length(x), function(i) {
          return(values[2*intervalParent[i] - 1, 3*intervalClass[i] - 2])  },  1);
        
        rootsMax = vapply(1:length(x), function(i) {
          return(values[2*intervalParent[i], 3*intervalClass[i] - 2])  },  1);
        
        # the extremes of values
        extremesMin = vapply(1:length(x), function(i) {
          return(values[2*intervalParent[i] - 1, 3*intervalClass[i] - 1])  },  1);
        
        extremesMax = vapply(1:length(x), function(i) {
          return(values[2*intervalParent[i], 3*intervalClass[i] - 1])  },  1);
        
        # the probability of the queues of values
        probQueuesMin = vapply(1:length(x), function(i) {
          return(values[2*intervalParent[i] - 1, 3*intervalClass[i]])  },  1);
        
        probQueuesMax = vapply(1:length(x), function(i) {
          return(values[2*intervalParent[i], 3*intervalClass[i]])  },  1);
        
        # the prediction of the used polynom
        polyVal = vapply(1:length(x), function(i) {
          p = polynom[intervalClass[i],,intervalParent[i]];
          p = p[!is.na(p)];
          return(predict(as.polynomial(p),x[i]))   },  1);
        
        # probability in each part
        probability1 = (x>=extremesMin & x<rootsMin) * probQueuesMin;
        probability2 = (x>=rootsMin & x<=rootsMax) * polyVal;
        probability3 = (x>rootsMax & x<=extremesMax) * probQueuesMax;
        
        # the sum of the three probabilities
        density = probability1 + probability2 + probability3;
        
      }
      
      
      else     # discrete or discretized
      {
        # the values of the variable depends on its type
        if (varType[posX]==0)   {   dataVal = prob[[3*posX]];   }
        else          # discretized
        {
          dataVal = 1:(length(prob[[3*posX]])-1);    # the values are 1, 2, 3, 4...
        }
        
        # the used part of the model
        dataMatrix = prob[[3*posX+2]];
        
        # obtain the position in the matrix of probabilities
        density = vapply(1:length(x), function(i) {
          if (x[i] %in% dataVal)   {   
            density = dataMatrix[dataVal==x[i],intervalClass[i],intervalParent[i]];   }
          else   {   density = 0;   }
          return(density)   },  1);
      }
    }
    
    
    
    return(density);
  }
    
  
  
  
  # matrix for the conditional probabilities in each variable in each observation for a
  # concrete part of the class
  probCond = matrix(NA, nrow=nObservations, ncol=nVar);
  
  # matrix for the conditional probabilities in each observation in each part of the class
  probClassValues = matrix(NA, nrow=nObservations, ncol=nClass);
  
  # matrix for the integral in each observation in each part of the class
  classIntegral = matrix(NA, nrow=nObservations, ncol=nClass);        # p(C)
  classExpectation = matrix(NA, nrow=nObservations, ncol=nClass);        # C * p(C)
  
  
  # for each part of the class
  for (i in 1:nClass)
  {
    # mid point of this part of the class
    c = (cutPointsClass[i] + cutPointsClass[i+1]) / 2;
    
    # for each prediction variable
    for (j in 1:nVar)
    {
      
      #cat('\n \n    Variable ',j)
      
      
      
      # each column of the matrix corresponds to a variable
      if (nParents[j]==0)      # the node variable
      {
        # calculate the probabilities for all the observations
        probCond[,j] = probXC(data[,j],c,j);
      }
      else   # 2 parents: a parent apart from the discrete class
      {
        # obtain the parent of the variable
        posParent = which(links[,j]==1);
        
        # calculate the probabilities for all the observations
        probCond[,j] = probXYC(data[,j],data[,posParent],c,j,posParent);
      }
      
      
      
      
      #if (j==1)
      #{
      #  cat('\n \n \n \n   Variable ',j)
      #  #cat('\n nParents = ',nParents[1]);
      #  cat('\n datos en la variable 1 = ',data[,j]);
      #  #cat('\n datos en la variable padre = ',data[,posParent]);
      #  cat('\n probcond = ',probCond[,j]);
      #}
      
      # probability of the class
      ### REVISAR
      #probCond[,nVar+1] = probC(c);
    }
    
    
    # multiply the probabilities of all the prediction variables
    probClassValues[,i] = apply(probCond, 1, prod);
    
    # evaluate the integral in the extremes of the current interval
    integralpCSub = diff(c(integralpC(cutPointsClass[i],1),
                           integralpC(cutPointsClass[i+1],1)));
    integralCpCSub = diff(c(integralpC(cutPointsClass[i],2),
                            integralpC(cutPointsClass[i+1],2)));
    ### CAMBIAR PARA AÑADIR LAS COLAS
    #integralpCSub = diff(predict(integralpC, c(cutPointsClass[i], cutPointsClass[i+1])));
    #integralCpCSub = diff(predict(integralCpC, c(cutPointsClass[i], cutPointsClass[i+1])));
    
    #cat('\n intervalo',i);
    #cat('\n probCond = ');   print(probCond);
    #cat('\n cutPointsClass = ',cutPointsClass);
    #cat('\n integral1 = ',integralpCSub);
    #cat('\n integral2 = ',integralCpCSub);
    #cat('\n probClassValues = ',probClassValues[,i]);
    
    
    # multiply the integral by the density of all the prediction variables
    classIntegral[,i] = probClassValues[,i] * integralpCSub;
    classExpectation[,i] = probClassValues[,i] * integralCpCSub;
  }
  
  
  # add the probability of each part of the class for each observation
  sumClassIntegral = apply(classIntegral,1,sum);
  sumClassExpectation = apply(classExpectation,1,sum);
  
  # divide the expectation by the probability
  predictions = sumClassExpectation / sumClassIntegral;
  
  
  #cat('\n prediction = ',predictions);
  return(predictions);
}
############################################################





  







