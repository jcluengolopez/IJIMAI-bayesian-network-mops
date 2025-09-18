

SimulatetMoP = function(poly,values,nData,verbose=FALSE)
{
  # it simulates one continuous variable estimated using a tMoP
  # INPUT: - poly = vector with the coefficients of the normalized polynomial
  #        - values = matrix with the min and max root and extremes and prob of the tails
  #        - nData = number of data to simulate
  # OUTPUT: - dataSim = data frame with the new simulation
  
  
  # simulate the probability
  probSim = runif(nData);
  
  # initialize the vector for the new data
  dataSim = rep(0,nData);
  
  # each value of the simulation
  for (i in 1:nData)
  {
    # inverse of a polynomial fit with the queues
    dataSim[i] = InversePolynomialFit(as.polynomial(poly),values,probSim[i]);
  }
  
  return(dataSim);
}
############################################################







Simulate = function(prob,nData,verbose=FALSE)
{
  # it creates a dataset of simulated observation using a model of Naive Bayes or TAN
  # INPUT: - prob = the TAN or Naive Bayes model
  #        - nData = number of simulated observations
  #        - verbose = show the 'cats' on screen
  # OUTPUT: - dataSim = data frame of the new observations
  
  
  # change the structure of the model
  prob2 = StructureModelIverse(prob);
  
  # number of variables, its names, its types and links
  varNames = prob$nodes$names;      links = prob2[[2]];
  varType = prob2[[1]];             nVar = length(varNames);
  
  cat('\n varType = ',varType);
  
  
  # sum the columns of links to get the number of parents of each variable
  nParents = apply(links,2,sum);
  
  # dataframe for the simulated data
  dataSim = data.frame(matrix(ncol=nVar,nrow=nData));
  colnames(dataSim) = varNames;
  
  # simulate the class
  if (varType[nVar]==0)        # classification
  {
    # values and matrix of probabilities of the class
    l = length(prob2);           classVal = prob2[[l-1]];          classProb = prob2[[l]];
    
    # simulate using the probabilities
    dataSim[,nVar] = sample(classVal,nData,replace=T,prob=classProb);
  }
  else        # regression
  {
    # values and polynomials of the class
    l = length(prob2);           values = prob2[[l-1]];          polynom = prob2[[l]];
    
    # used polynomial with the used coefficients
    nCoef = length(polynom[!is.na(polynom)]);
    p = as.polynomial(polynom[1:nCoef]);
    
    # simulate the probability
    probSim = runif(nData);
    
    # inverse of a polynomial fit with the queues of each simulated probability
    for (i in 1:nData)
    {
      dataSim[i,nVar] = InversePolynomialFit(p,values,probSim[i]);
    }
  }
  
  # vector with the simulated variables (the class is already done)
  simulated = c(rep(FALSE,nVar-1),TRUE)
  
  # repeat the process until all the variables are simulated
  while (FALSE %in% simulated)
  {
    # check which variable can be simulated
    validVar = apply(links,2,function(x) !(FALSE %in% simulated[which(x==1)]))
    
    # valid variables that have not been simulated yet (use the first one)
    simVar = which(simulated==FALSE & validVar==TRUE)[1];
    
    # parents of the selected variable
    parents = which(links[,simVar]==1);
    
    #cat('\n simVar = ',simVar);
    #cat('\n parents = ',parents);
    
    
    
    # simulate the selected variable and save in the existing dataframe
    if (varType[nVar]==0)        # classification
    {
      dataSim = SimulateVarClassification(prob2,simVar,parents,dataSim,verbose);
    }
    else        # regression
    {
      dataSim = SimulateVarRegression(prob2,simVar,parents,dataSim,verbose);
    }
    
    simulated[simVar] = TRUE;
  }
  
  return(dataSim)
}
############################################################








SimulateVarClassification = function(prob,simVar,parents,dataSim,verbose)
{
  # it simulates one variable in a classification model using the simulation of its 
  # parents and the model
  # INPUT: - prob = the TAN or Naive Bayes internal model
  #        - simVar = position of the variable to simulate
  #        - parents = position of the parents of the variable to simulate
  #        - dataSim = dataframe with the variables that are already simulated
  #        - verbose = show the 'cats' on screen
  # OUTPUT: - dataSim = data frame uploaded with the new simulation
  
  
  # number of variables and its types
  varType = prob[[1]];          nVar = length(varType);
  
  # number of observations
  nData = dim(dataSim)[1];
  
  # values of the class
  classVal = prob[[length(prob)-1]];
  
  # the model depends on the parents
  if (length(parents)==1)        # 1 parent (it mus be the node variable)
  {
    # the model depends on the type of variable
    if (varType[simVar]==1)        # when the variable is continuous
    {
      cat('\n Entra continua sin padres con la clase discreta. Variable ',simVar);
      
      
      # conditioned polynomials and the matrix of values
      polynom = prob[[3*simVar+1]];            values = prob[[3*simVar]];
      
      # simulate the probability
      probSim = runif(nData);
      
      # each value of the simulation
      for (i in 1:nData)
      {
        # position that corresponds to the simulated value in the class
        posClass = which(classVal == dataSim[i,nVar]);
        
        # used part of the values matrix
        v = values[,(3*posClass-2):(3*posClass)];
        
        # used polynomial with the used coefficients
        aux = polynom[posClass,]
        nCoef = length(aux[!is.na(aux)]);
        p = as.polynomial(aux[1:nCoef]);
        
        # inverse of a polynomial fit with the queues
        dataSim[i,simVar] = InversePolynomialFit(p,v,probSim[i]);
      }
    }
    
    else         # when the variable is discrete or discretized
    {
      cat('\n Entra discreta sin padres con la clase discreta. Variable ',simVar);
      
      
      # the values of the variable depends on the type
      if (varType[simVar]==0)   {   dataVal = prob[[3*simVar]];   }
      else          # discretized
      {
        # the values are 1 ... number of split points - 1
        dataVal = 1:(length(prob[[3*simVar]])-1);
      }
      
      #cat('\n dataVal = ',dataVal);
      #cat('\n classaVal = ',classVal);
      #print(dataSim);
      
      # extract its conditional probabilities
      dataMatrix = prob[[3*simVar+1]];
      
      # each value of the simulation
      for (i in 1:dim(dataSim)[1])
      {
        # column that corresponds to the simulated value in the class
        posClass = (classVal == dataSim[i,nVar]);
        
        # simulate using the probabilities of that value
        dataSim[i,simVar] = sample(dataVal,1,replace=T,prob=dataMatrix[,posClass]);
      }
    }
  }
  
  
  else        # 2 parents: a parent apart from the discrete class
  {
    # the model depends on the type of parent
    if (varType[parents[1]]==1)   # 2 parents (the discrete class and another continuous)
    {
      # the model depends on the type of variable
      if (varType[simVar]==1)        # when the variable is continuous
      {
        cat('\n Entra continua con un padre continuo y clase discreta. Variable ',simVar);
        
        
        # conditioned polynomials and the matrix of values
        polynom = prob[[3*simVar+2]];            values = prob[[3*simVar]];
        
        # simulate the probability
        probSim = runif(nData);
        
        # each value of the simulation
        for (i in 1:nData)
        {
          # position that corresponds to the simulated value in the class
          posClass = which(classVal == dataSim[i,nVar]);
          
          # cut points of the discretization of the continuous parent
          aux = prob[[3*simVar+1]][posClass,];
          nPoints = length(aux[!is.na(aux)]);
          splitPoints = aux[1:nPoints];
          
          # position of the discretization where the observation belongs
          aux = max(which(dataSim[i,parents[1]]>=splitPoints));
          interval = min(aux,length(splitPoints)-1);
          
          # used part of the values matrix
          v = values[(2*posClass-1):(2*posClass),(3*interval-2):(3*interval)];
          
          # used polynomial with the used coefficients
          aux = polynom[interval,,posClass]
          nCoef = length(aux[!is.na(aux)]);
          p = as.polynomial(aux[1:nCoef]);
          
          # inverse of a polynomial fit with the queues
          dataSim[i,simVar] = InversePolynomialFit(p,v,probSim[i]);
        }
      }
      
      
      else         # when the variable is discrete or discretized
      {
        cat('\n Entra discreta con un padre continuo y clase discreta. Variable ',simVar);
        
        
        # the values of the variable depends on the type
        if (varType[simVar]==0)   {   dataVal = prob[[3*simVar]];   }
        else          # discretized
        {
          # the values are 1 ... number of split points - 1
          dataVal = 1:(length(prob[[3*simVar]])-1);
        }
        
        # extract its conditional probabilities and the cutPoints
        dataMatrix = prob[[3*simVar+2]];
        
        # each value of the simulation
        for (i in 1:dim(dataSim)[1])
        {
          # position that corresponds to the simulated value in the class
          posClass = which(classVal == dataSim[i,nVar]);
          
          # cut points of the discretization of the continuous parent
          aux = prob[[3*simVar+1]][posClass,];
          nPoints = length(aux[!is.na(aux)]);
          splitPoints = aux[1:nPoints];
          
          # position of the discretization where the observation belongs
          aux = max(which(dataSim[i,parents[1]]>=splitPoints));
          interval = min(aux,length(splitPoints)-1);
          
          # used probabilities
          p = dataMatrix[,interval,posClass];
          
          # simulate using the probabilities of that value
          dataSim[i,simVar] = sample(dataVal,1,replace=T,prob=p);
        }
      }
    }
    
    
    else    # 2 parents (the discrete class and another discrete or discretized)
    {
      # values of the other parent
      if (varType[parents[1]]==0)   {   parentVal = prob[[3*parents[1]]];   }
      else          # discretized
      {
        # the values are 1 ... number of split points - 1
        parentVal = 1:(length(prob[3*parents[1]])-1);
      }
      
      # the model depends on the type of variable
      if (varType[simVar]==1)        # when the variable is continuous
      {
        cat('\n Entra continua con un padre discreto y clase discreta. Variable ',simVar);
        
        
        # conditioned polynomials and the matrix of values
        polynom = prob[[3*simVar+1]];            values = prob[[3*simVar]];
        
        # simulate the probability
        probSim = runif(nData);
        
        # each value of the simulation
        for (i in 1:nData)
        {
          # position that corresponds to the simulated value in the class
          posClass = which(classVal == dataSim[i,nVar]);
          
          # position that corresponds to the simulated value in the other parent
          posParent = which(parentVal == dataSim[i,parents[1]]);
          
          # used part of the values matrix
          v = values[(2*posParent-1):(2*posParent),(3*posClass-2):(3*posClass)];
          
          # used polynomial with the used coefficients
          aux = polynom[posParent,,posClass]
          nCoef = length(aux[!is.na(aux)]);
          p = as.polynomial(aux[1:nCoef]);
          
          # inverse of a polynomial fit with the queues
          dataSim[i,simVar] = InversePolynomialFit(p,v,probSim[i]);
        }
      }
      
      else         # when the variable is discrete or discretized
      {
        cat('\n Entra discreta con un padre discreto y clase discreta. Variable ',simVar);
        
        
        # the values of the variable depends on the type
        if (varType[simVar]==0)   {   dataVal = prob[[3*simVar]];   }
        else          # discretized
        {
          # the values are 1 ... number of split points - 1
          dataVal = 1:(length(prob[[3*simVar]])-1);
        }
        
        # extract its conditional probabilities and the cutPoints
        dataMatrix = prob[[3*simVar+2]];
        
        # each value of the simulation
        for (i in 1:dim(dataSim)[1])
        {
          # position that corresponds to the simulated value in the class
          posClass = which(classVal == dataSim[i,nVar]);
          
          # position that corresponds to the simulated value in the other parent
          posParent = which(parentVal == dataSim[i,parents[1]]);
          
          # used probabilities
          p = dataMatrix[,posParent,posClass];
          
          # simulate using the probabilities of that value
          dataSim[i,simVar] = sample(dataVal,1,replace=T,prob=p);
        }
      }
    }
  }
  
  return(dataSim)
}
############################################################









SimulateVarRegression = function(prob,simVar,parents,dataSim,verbose)
{
  # it simulates one variable in a regression model using the simulation of its 
  # parents and the model
  # INPUT: - prob = the TAN or Naive Bayes internal model
  #        - simVar = position of the variable to simulate
  #        - parents = position of the parents of the variable to simulate
  #        - dataSim = dataframe with the variables that are already simulated
  #        - verbose = show the 'cats' on screen
  # OUTPUT: - dataSim = data frame uploaded with the new simulation
  
  
  # number of variables and its types
  varType = prob[[1]];          nVar = length(varType);
  
  # number of observations
  nData = dim(dataSim)[1];
  
  # the model depends on the parents
  if (length(parents)==1)        # 1 parent (it mus be the node variable)
  {
    # the model depends on the type of variable
    if (varType[simVar]==1)        # when the variable is continuous
    {
      cat('\n Entra continua sin padres con la clase continua Variable ',simVar);
      
      
      # conditioned polynomials and the matrix of values
      polynom = prob[[3*simVar+2]];            values = prob[[3*simVar]];
      
      # simulate the probability
      probSim = runif(nData);
      
      # each value of the simulation
      for (i in 1:nData)
      {
        # cut points of the discretization of the continuous class
        aux = prob[[3*simVar+1]];
        nPoints = length(aux[!is.na(aux)]);
        splitPoints = aux[1:nPoints];
        
        # position of the discretization where the observation belongs
        aux = max(which(dataSim[i,parents]>=splitPoints));
        interval = min(aux,length(splitPoints)-1);
        
        # used part of the values matrix
        v = values[,(3*interval-2):(3*interval)];
        
        # used polynomial with the used coefficients
        aux = polynom[interval,];
        nCoef = length(aux[!is.na(aux)]);
        p = as.polynomial(aux[1:nCoef]);
        
        # inverse of a polynomial fit with the queues
        dataSim[i,simVar] = InversePolynomialFit(p,v,probSim[i]);
      }
    }
    
    
    else         # when the variable is discrete or discretized
    {
      cat('\n Entra discreta sin padres con la clase continua Variable ',simVar);
      
      
      # the values of the variable depends on the type
      if (varType[simVar]==0)   {   dataVal = prob[[3*simVar]];   }
      else          # discretized
      {
        # the values are 1 ... number of split points - 1
        dataVal = 1:(length(prob[[3*simVar]])-1);
      }
      
      # extract its conditional probabilities and the cutPoints
      dataMatrix = prob[[3*simVar+2]];
      
      # each value of the simulation
      for (i in 1:dim(dataSim)[1])
      {
        # cut points of the discretization of the continuous class
        aux = prob[[3*simVar+1]];
        nPoints = length(aux[!is.na(aux)]);
        splitPoints = aux[1:nPoints];
        
        # position of the discretization where the observation belongs
        aux = max(which(dataSim[i,parents[1]]>=splitPoints));
        interval = min(aux,length(splitPoints)-1);
        
        # used probabilities
        p = dataMatrix[,interval];
        
        # simulate using the probabilities of that value
        dataSim[i,simVar] = sample(dataVal,1,replace=T,prob=p);
      }
    }
  }
  
  
  else        # 2 parents: a parent apart from the continuous class
  {
    # the model depends on the type of parent
    if (varType[parents[1]]==1)   # 2 continuous parents (the class and another continuous)
    {
      # min and max value of the continuous class (in all of its divisions)
      classMin = prob[[3*nVar]][1,2];               classMax = prob[[3*nVar]][2,2];
      
      # calculate all the possible cut points in the discretization of the class
      distClass = (classMax - classMin) / 8;
      splitPointsClass = classMin + distClass*(0:8);
      
      # min and max value of the parent (in all of its divisions)
      parentMin = prob[[3*parents[1]]][1,2];         parentMax = prob[[3*parents[1]]][2,2];
      
      # calculate all the possible cut points in the discretization of the parent
      distParent = (parentMax - parentMin) / 8;
      splitPointsParent = parentMin + distParent*(0:8);
      
      # the model depends on the type of variable
      if (varType[simVar]==1)        # when the variable is continuous
      {
        cat('\n Entra continua con un padre continuo y clase continua. Variable ',simVar);
        
        
        # conditioned polynomials and the matrix of values
        polynom = prob[[3*simVar+2]];            values = prob[[3*simVar]];
        
        # simulate the probability
        probSim = runif(nData);
        
        # each value of the simulation
        for (i in 1:nData)
        {
          # position of the discretization of the parent where the observation belongs
          aux = max(which(dataSim[i,nVar]>=splitPointsClass));
          posClass = min(aux,length(splitPointsClass)-1);
          
          # position of the discretization of the parent where the observation belongs
          aux = max(which(dataSim[i,parents[1]]>=splitPointsParent));
          posParent = min(aux,length(splitPointsParent)-1);
          
          # position in the matrix
          n = prob[[3*simVar+1]][posParent,posClass];
          
          # used part of the values matrix
          v = values[,(3*n-2):(3*n)];
          
          # used polynomial with the used coefficients
          aux = polynom[n,];
          nCoef = length(aux[!is.na(aux)]);
          p = as.polynomial(aux[1:nCoef]);
          
          # inverse of a polynomial fit with the queues
          dataSim[i,simVar] = InversePolynomialFit(p,v,probSim[i]);
        }
      }
      
      
      else         # when the variable is discrete or discretized
      {
        cat('\n Entra discreta con un padre continuo y clase continua. Variable ',simVar);
        
        
        # the values of the variable depends on the type
        if (varType[simVar]==0)   {   dataVal = prob[[3*simVar]];   }
        else          # discretized
        {
          # the values are 1 ... number of split points - 1
          dataVal = 1:(length(prob[[3*simVar]])-1);
        }
        
        # extract its conditional probabilities and the cutPoints
        dataMatrix = prob[[3*simVar+2]];
        
        # each value of the simulation
        for (i in 1:nData)
        {
          # position of the discretization of the parent where the observation belongs
          aux = max(which(dataSim[i,nVar]>=splitPointsClass));
          posClass = min(aux,length(splitPointsClass)-1);
          
          # position of the discretization of the parent where the observation belongs
          aux = max(which(dataSim[i,parents[1]]>=splitPointsParent));
          posParent = min(aux,length(splitPointsParent)-1);
          
          # position in the matrix
          n = prob[[3*simVar+1]][posParent,posClass];
          
          # used probabilities
          p = dataMatrix[n,];
          
          # simulate using the probabilities of that value
          dataSim[i,simVar] = sample(dataVal,1,replace=T,prob=p);
        }
      }
    }
    
    
    else    # 2 parents (the continuous class and another discrete or discretized)
    {
      # values of the other parent
      if (varType[parents[1]]==0)   {   parentVal = prob[[3*parents[1]]];   }
      else          # discretized
      {
        # the values are 1 ... number of split points - 1
        parentVal = 1:(length(prob[3*parents[1]])-1);
      }
      
      # the model depends on the type of variable
      if (varType[simVar]==1)        # when the variable is continuous
      {
        cat('\n Entra continua con un padre discreto y clase continua. Variable ',simVar);
        
        
        # conditioned polynomials and the matrix of values
        polynom = prob[[3*simVar+2]];            values = prob[[3*simVar]];
        
        # simulate the probability
        probSim = runif(nData);
        
        # each value of the simulation
        for (i in 1:nData)
        {
          # position that corresponds to the simulated value in the discrete parent
          posParent = which(parentVal == dataSim[i,parents[1]]);
          
          # cut points of the discretization of the continuous class
          aux = prob[[3*simVar+1]][posParent,];
          nPoints = length(aux[!is.na(aux)]);
          splitPoints = aux[1:nPoints];
          
          # position of the discretization where the observation belongs
          aux = max(which(dataSim[i,nVar]>=splitPoints));
          interval = min(aux,length(splitPoints)-1);
          
          # used part of the values matrix
          v = values[(2*posParent-1):(2*posParent),(3*interval-2):(3*interval)];
          
          # used polynomial with the used coefficients
          aux = polynom[interval,,posParent]
          nCoef = length(aux[!is.na(aux)]);
          p = as.polynomial(aux[1:nCoef]);
          
          # inverse of a polynomial fit with the queues
          dataSim[i,simVar] = InversePolynomialFit(p,v,probSim[i]);
        }
      }
      
      
      else         # when the variable is discrete or discretized
      {
        cat('\n Entra discreta con un padre discreto y clase continua. Variable ',simVar);
        
        
        # the values of the variable depends on the type
        if (varType[simVar]==0)   {   dataVal = prob[[3*simVar]];   }
        else          # discretized
        {
          # the values are 1 ... number of split points - 1
          dataVal = 1:(length(prob[[3*simVar]])-1);
        }
        
        # extract its conditional probabilities and the cutPoints
        dataMatrix = prob[[3*simVar+2]];
        
        # each value of the simulation
        for (i in 1:dim(dataSim)[1])
        {
          # position that corresponds to the simulated value in the discrete parent
          posParent = which(parentVal == dataSim[i,parents[1]]);
          
          # cut points of the discretization of the continuous class
          aux = prob[[3*simVar+1]][posParent,];
          nPoints = length(aux[!is.na(aux)]);
          splitPoints = aux[1:nPoints];
          
          # position of the discretization where the observation belongs
          aux = max(which(dataSim[i,nVar]>=splitPoints));
          interval = min(aux,length(splitPoints)-1);
          
          # used probabilities
          p = dataMatrix[,interval,posParent];
          
          # simulate using the probabilities of that value
          dataSim[i,simVar] = sample(dataVal,1,replace=T,prob=p);
        }
      }
    }
  }
  
  return(dataSim)
}
############################################################




  
  



InversePolynomialFit = function(polynom,values,probSim)
{
  # it calculates the inverse of a continuous fit with two tails and a polynomial
  # INPUT: - values = matrix with roots, extremes and prob. of the queues
  #        - polynomial = polynomial of the central part of the fit
  #        - probSim = simulated probability
  # OUTPUT: - dataSim = data frame of the new observations
  
  
  # separate roots, extremes and probabilities of the queues
  roots = values[,1];        extremes = values[,2];        probQueues = values[,3];
  
  # cumulative probability in each queue
  probQueue1 = probQueues[1] * abs(roots[1] - extremes[1]);
  probQueue2 = probQueues[2] * abs(extremes[2] - roots[2]);
  
  # check if the simulation prob corresponds to the polynomial or the queues
  if (probSim < probQueue1)          # first queue
  {
    # probability function of the queue starting in 0 (probQ * x - probQ * extr)
    F1 = integral(polynomial(probQueues[1])) - probQueues[1] * extremes[1];
    
    # function to solve (F1 - probSim)
    F1Solve = F1 - polynomial(probSim);
    
    # real zeroes of that function (only one because it's a linear function)
    zeroes = solve(F1Solve);
    dataSim = Re(zeroes[Im(zeroes) == 0]);
  }
  else if (probSim > (1-probQueue2))          # second queue
  {
    # probability function of the queue finishing in 1 (prob * x - probQ * extr + 1)
    F2 = integral(polynomial(probQueues[2])) - probQueues[2] * extremes[2] + 1;
    
    # function to solve (F2 - probSim)
    F2Solve = F2 - polynomial(probSim);
    
    # real zeroes of that function (only one because it's a linear function)
    zeroes = solve(F2Solve);
    dataSim = Re(zeroes[Im(zeroes) == 0]);
  }
  else          # polynomial
  {
    # probability function of the central part starting in probQ1 (p - p(root) + probQueue1)
    F3 = integral(polynom) - predict(integral(polynom),roots[1]) + probQueue1;
    
    # function to solve (F3 - probSim)
    F3Solve = F3 - polynomial(probSim);
    
    # real zeroes of that function (only one because it is increasing)
    zeroes = solve(F3Solve);
    realRoots = Re(zeroes[Im(zeroes) == 0]);
    dataSim = realRoots[realRoots>roots[1] & realRoots<roots[2]];
  }
  
  return(dataSim)
}
############################################################








predicttMoP = function(tMoP,x)
{
  # it calculates the probability of a x-value in a tMoP
  # INPUT: - tMoP = model of a tMoP
  #        - x = vector of values to predict
  # OUTPUT: - y = predicted values
  
  
  # vector for the predictions
  y = vector(mode='numeric',length=length(x));
  
  # extract information from the tMoP
  roots = tMoP[[1]];        poly = tMoP[[2]];        probTails = tMoP[[3]];
  
  # values of x to predict using the polynomial
  selection1 = x>roots[1] & x<roots[2]; 
  y[selection1] = predict(as.polynomial(poly),x[selection1]);
  
  # values of x in the left tail
  selection2 = x<=roots[1];            y[selection2] = probTails[1];
  
  # values of x in the right tail
  selection3 = x>=roots[2];            y[selection3] = probTails[2];
  
  return(y);
}
############################################################









