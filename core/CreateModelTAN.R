

CreateModelTANClassification = function(data,dataVal,dataMin,dataMax,varType,splitPoints,
                               maxDegree,mutualInfoClass,verbose=FALSE,plots=FALSE)
{
  # it creates the TAN model for classification with a discrete class using data 
  # (the selected variables in each case). It contains the conditional probabilities
  # in discrete variables and the conditional polynomials in continuous variables
  # INPUT: - data = a data-frame with the data of study
  #        - dataVal = list with the different values of the discrete variables
  #        - dataMin = the minimum value in the variable
  #        - dataMax = the maximum value in the variable
  #        - varType = a vector with 0 if the variable is discrete, 1 if it's
  #       continuous and 2 if it has to be discretized. Program calculates it
  #        - splitPoints = the cut points in discretized variables (it can be NULL)
  #        - maxDegree = maximum degree for the polynomial adjust
  #        - mutualInfoClass = matrix with the mutual information of each pair of variables
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - prob: list with 3 elements per variable, the links, the type of variables
  #       and the information of the class
  
  
  # the class is the last position of the data frame
  C = data[,ncol(data)];
  
  # create the list for probabilities and polynomials (3 each variable and the class, 2) 
  prob = vector("list",((ncol(data)-1)*3+4));
  
  # matrix with the probabilities of C (last position)
  classVal = dataVal[[length(dataVal)]];
  prob[[length(prob)-1]] = classVal;
  prob[[length(prob)]] = DiscreteProb(data[,ncol(data)],classVal,verbose)[,2];
  
  # create the maximum spanning tree avoiding loops
  links = CreateTree(mutualInfoClass,verbose);
  
  # put the structure of the model in prob
  prob[[1]] = varType;            prob[[2]] = links;
  
  # from now on the matrix of links avoid the last column and row referring to the class
  links = links[-ncol(data),-ncol(data)];
  
  # sum the columns of links to get the number of parents of each variable
  nParents = apply(links,2,sum);
  
  # for each variable except the class
  for (i in 1:(ncol(data)-1))
  {
    # the model depends on the parents
    if (nParents[i]==0)   # 0 parent (it mus be the node variable)
    {
      # the adjust depends on the type of variable
      if (varType[i]==0)   # when the variable is discrete
      {
        # the function that calculates the conditional probability and the values
        p = Adjust1ParentDD(data[,i],dataVal[[i]],C,classVal,verbose,plots);
        
        # insert the results in the list
        prob[[3*i]] = dataVal[[i]];
        prob[[3*i+1]] = p[[1]];
      }
      else if (varType[i]==1)   # when the variable is continuous
      {
        # division for the plots
        if (plots==TRUE)   {
          par(mfrow=c(length(classVal)+1,2),mar=c(1,1,1,1),oma=c(0.4,0.4,0.4,0.4));   }
        
        # calculate the polynomials conditioned to the class
        p = Adjust1ParentCD(data[,i],dataMin[i],dataMax[i],C,classVal,i,maxDegree,verbose,
                            plots);
        
        # insert the results in the list
        prob[[3*i]] = p[[1]];
        prob[[3*i+1]] = p[[2]];
      }
      else   # when the variable is discretized
      {
        # calculate the conditional probabilities and the values
        p = Adjust1ParentDD(data[,i],dataVal[[i]],C,classVal,verbose,plots);
        
        # insert the cut points in the list prob. When there're no cut points
        prob[[3*i]] = splitPoints[[i]];
        
        # insert the probabilities in the list prob
        prob[[3*i+1]] = p[[1]];
      }
    }
    else   # 2 parents: a parent apart from the class
    {
      # obtain the parent of the variable
      posParent = which(links[,i]==1);
      
      # the model depends on the type of parent
      if (varType[posParent]==1)   # 2 parents (the class and another continuous)
      {
        # the adjust depends on the type of variable
        if (varType[i]==0)   # when the variable is discrete
        {
          # the function that calculates the conditional probability and the values
          p = Adjust2ParentsDCD(data[,i],dataVal[[i]],data[,posParent],dataMin[posParent],
                                dataMax[posParent],C,classVal,verbose,plots);
          
          # insert the results in the list
          prob[[3*i]] = dataVal[[i]];
          prob[[3*i+1]] = p[[2]];    # cut points
          prob[[3*i+2]] = p[[1]];    # dataMatrix
        }
        else if (varType[i]==1)   # when the variable is continuous
        {
          if (verbose==TRUE)  {
            cat('\n variable ',i,'con padre',posParent,' y la clase discreta');
            cat('\n padre ',posParent);          
            cat('\n padreMin = ',dataMin[posParent],' y padreMax = ',dataMax[posParent])
          }
          
          # calculate the polynomials conditioned to the class
          p = Adjust2ParentsCCD(data[,i],dataMin[i],dataMax[i],data[,posParent],
                        dataMin[posParent],dataMax[posParent],C,classVal,i,maxDegree,
                        verbose,plots);
          
          # insert the results in the list
          prob[[3*i]] = p[[1]];       # values
          prob[[3*i+1]] = p[[3]];     # cutPoints
          prob[[3*i+2]] = p[[2]];     # polynomial
        }
        else   # when the variable is discretized
        {
          # the function that calculates the conditional probability and the values
          p = Adjust2ParentsDCD(data[,i],dataVal[[i]],data[,posParent],dataMin[posParent],
                                dataMax[posParent],C,classVal,verbose,plots);
          
          # insert the cut points in the list prob. When there're no cut points
          prob[[3*i]] = splitPoints[[i]];
          
          # insert the probabilities and the cut points in the list prob
          prob[[3*i+1]] = p[[2]];
          prob[[3*i+2]] = p[[1]];
        }
      }
      else   # 2 parents (the class and another discrete)
      {
        # the adjust depends on the type of variable
        if (varType[i]==0)   # when the variable is discrete
        {
          # calculate the conditional probabilities and the values
          p = Adjust2ParentsDDD(data[,i],dataVal[[i]],data[,posParent],dataVal[[posParent]],
                                C,classVal,verbose,plots);
          
          # insert the results in the list
          prob[[3*i]] = dataVal[[i]];
          prob[[3*i+1]] = p[[1]];              # dataMatrix
        }
        else if (varType[i]==1)   # when the variable is continuous
        {
          # calculate the polynomials conditioned to the class
          p = Adjust2ParentsCDD(data[,i],dataMin[i],dataMax[i],data[,posParent],
                                dataVal[[posParent]],C,classVal,i,maxDegree,verbose,plots);
          
          # insert the results in the list
          prob[[3*i]] = p[[1]];       # values
          prob[[3*i+1]] = p[[2]];     # polynomial
        }
        else   # when the variable is discretized
        {
          # calculate the conditional probabilities and the values
          p = Adjust2ParentsDDD(data[,i],dataVal[[i]],data[,posParent],
                                dataVal[[posParent]],C,classVal,verbose,plots);
          
          # insert the cut points in the list prob. When there're no cut points
          prob[[3*i]] = splitPoints[[i]];
          
          # insert the probabilities in the list prob
          prob[[3*i+1]] = p[[1]];
        }
      }
    }
  }
  
  return(prob);
}
############################################################








CreateModelTANRegression = function(data,dataVal,dataMin,dataMax,varType,splitPoints,
                           maxDegree,mutualInfoClass,verbose=FALSE,plots=FALSE)
{
  # it creates the TAN model for regression with a continuous class using data 
  # (the selected variables in each case). It contains the conditional probabilities
  # in discrete variables and the conditional polynomials in continuous variables
  # INPUT: - data = a data-frame with the data of study
  #        - dataVal = list with the different values of the discrete variables
  #        - dataMin = the minimum value in the variable
  #        - dataMax = the maximum value in the variable
  #        - varType = a vector with 0 if the variable is discrete, 1 if it's
  #       continuous and 2 if it has to be discretized. Program calculates it
  #        - splitPoints = the cut points in discretized variables (it can be NULL)
  #        - maxDegree = maximum degree for the polynomial adjust
  #        - mutualInfoClass = matrix with the mutual information of each pair of variables
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - prob: list with 3 elements per variable, the links, the type of variables
  #      and the information of the class
  
  
  # the class is the last position of the data frame
  C = data[,ncol(data)];
  classMin = dataMin[[ncol(data)]];             classMax = dataMax[[ncol(data)]];
  
  # create the list for probabilities and polynomials (3 each variable and the class, 2) 
  prob = vector("list",((ncol(data)-1)*3+4));
  
  # matrix with the probabilities of C (last position). Strings are not factors
  adjustClass = Adjust1ParentCD(C,classMin,classMax,rep(0,length(C)),0,ncol(data),
                                maxDegree,verbose,plots);
  
  prob[[length(prob)-1]] = adjustClass[[1]];
  prob[[length(prob)]] = adjustClass[[2]];
  
  # create the maximum spanning tree avoiding loops
  links = CreateTree(mutualInfoClass,verbose);
  
  # put the structure of the model in prob
  prob[[1]] = varType;
  prob[[2]] = links;
  
  # from now on the matrix of links avoid the last column and row referring to the class
  links = links[-ncol(data),-ncol(data)];
  
  # sum the columns of links to get the number of parents of each variable
  nParents = apply(links,2,sum);
  
  # for each variable except the class
  for (i in 1:(ncol(data)-1))
  {
    # the model depends on the parents
    if (nParents[i]==0)   # 0 parent (it mus be the node variable)
    {
      # the adjust depends on the type of variable
      if (varType[i]==0)   # when the variable is discrete
      {
        if (verbose)   {
          cat('\n VARIABLE ',i,' DISCRETA SIN PADRES');   }
        
        # calculate the conditional probabilities and the values
        p = Adjust1ParentDC(data[,i],dataVal[[i]],C,classMin,classMax,verbose,plots);
        
        # insert the results in the list
        prob[[3*i]] = dataVal[[i]];
        prob[[3*i+1]] = p[[2]];    # cut points
        prob[[3*i+2]] = p[[1]];    # dataMatrix
      }
      else if (varType[i]==1)   # when the variable is continuous
      {
        if (verbose)   {
          cat('\n VARIABLE ',i,' CONTINUA SIN PADRES');   }
        
        # calculate the polynomials conditioned to the class
        p = Adjust1ParentCC(data[,i],dataMin[i],dataMax[i],C,classMin,classMax,i,
                            maxDegree,verbose,plots);
        
        # insert the results in the list
        prob[[3*i]] = p[[1]];       # values
        prob[[3*i+1]] = p[[3]];     # cutPoints
        prob[[3*i+2]] = p[[2]];     # polynomial
      }
      else   # when the variable is discretized
      {
        if (verbose)   {
          cat('\n VARIABLE ',i,' DISCRETIZADA SIN PADRES');   }
        
        # calculate the conditional probabilities and the values
        p = Adjust1ParentDC(data[,i],dataVal[[i]],C,classMin,classMax,verbose,plots);
        
        # insert the cut points in the list prob. When there're no cut points
        prob[[3*i]] = splitPoints[[i]];
        
        # insert the results in the list
        prob[[3*i+1]] = p[[2]];    # cut points
        prob[[3*i+2]] = p[[1]];    # dataMatrix
      }
    }
    else   # 2 parents: a parent apart from the class
    {
      # obtain the parent of the variable
      posParent = which(links[,i]==1);
      
      # the model depends on the type of parent
      if (varType[posParent]==1)   # 2 parents (the class and another continuous)
      {
        # the adjust depends on the type of variable
        if (varType[i]==0)   # when the variable is discrete
        {
          if  (verbose)   {
            cat('\n VARIABLE ',i,' DISCRETA CON 2 PADRES CONTINUOS');   }
          
          # calculate the conditional probabilities and the values
          p = Adjust2ParentsDCC(data[,i],dataVal[[i]],data[,posParent],dataMin[posParent],
                      dataMax[posParent],C,classMin,classMax,i,maxDegree,verbose,plots);
          
          # insert the results in the list
          prob[[3*i]] = dataVal[[i]];
          prob[[3*i+1]] = p[[2]];    # cut points
          prob[[3*i+2]] = p[[1]];    # dataMatrix
        }
        else if (varType[i]==1)   # when the variable is continuous
        {
          if (verbose)   {
            cat('\n VARIABLE ',i,' CONTINUA CON 2 PADRES CONTINUOS');   }
          
          # calculate the polynomials conditioned to the class 
          p = Adjust2ParentsCCC(data[,i],dataMin[i],dataMax[i],data[,posParent],
                                dataMin[posParent],dataMax[posParent],C,classMin,classMax,
                                i,maxDegree,verbose,plots);
          
          # insert the results in the list
          prob[[3*i]] = p[[1]];       # values
          prob[[3*i+1]] = p[[3]];     # cutPoints
          prob[[3*i+2]] = p[[2]];     # polynomial
        }
        else   # when the variable is discretized
        {
          if (verbose)   {
            cat('\n VARIABLE ',i,' DISCRETIZADA CON 2 PADRES CONTINUOS');   }
          
          # calculate the conditional probabilities and the values
          p = Adjust2ParentsDCC(data[,i],dataVal[[i]],data[,posParent],
                                dataMin[posParent],dataMax[posParent],C,classMin,classMax,
                                i,maxDegree,verbose,plots);
          
          # insert the cut points in the list prob. When there're no cut points
          prob[[3*i]] = splitPoints[[i]];
          
          # insert the results in the list
          prob[[3*i+1]] = p[[2]];    # cut points
          prob[[3*i+2]] = p[[1]];    # dataMatrix
        }
      }
      else   # 2 parents (the class and another discrete)
      {
        # the adjust depends on the type of variable
        if (varType[i]==0)   # when the variable is discrete
        {
          if (verbose)   {
            cat('\n VARIABLE ',i,' DISCRETA CON 1 PADRE DICRETO Y 1 PADRE CONTINUO');   }
          
          # calculate the conditional probabilities and the values
          p = Adjust2ParentsDCD(data[,i],dataVal[[i]],C,classMin,classMax,
                                data[,posParent],dataVal[[posParent]],verbose,plots);
          
          # insert the results in the list
          prob[[3*i]] = dataVal[[i]];
          prob[[3*i+1]] = p[[2]];    # cut points
          prob[[3*i+2]] = p[[1]];    # dataMatrix
        }
        else if (varType[i]==1)   # when the variable is continuous
        {
          if (verbose)   {
            cat('\n VARIABLE ',i,' CONTINUA CON 1 PADRE DICRETO Y 1 PADRE CONTINUO');   }
          
          # calculate the polynomials conditioned to the class 
          p = Adjust2ParentsCCD(data[,i],dataMin[i],dataMax[i],C,classMin,classMax,
                        data[,posParent],dataVal[[posParent]],i,maxDegree,verbose,plots);
          
          # insert the results in the list
          prob[[3*i]] = p[[1]];       # values
          prob[[3*i+1]] = p[[3]];     # cutPoints
          prob[[3*i+2]] = p[[2]];     # polynomial
        }
        else   # when the variable is discretized
        {
          if (verbose)   {
            cat('\n VARIABLE ',i,' DISCRETIZADA CON 1 PADRE DICRETO Y 1 PADRE CONTINUO'); }
          
          # calculate the conditional probabilities and the values
          p = Adjust2ParentsDCD(data[,i],dataVal[[i]],C,classMin,classMax,
                                data[,posParent],dataVal[[posParent]],verbose,plots);
          
          # insert the cut points in the list prob. When there're no cut points
          prob[[3*i]] = splitPoints[[i]];
          
          # insert the results in the list
          prob[[3*i+1]] = p[[2]];    # cut points
          prob[[3*i+2]] = p[[1]];    # dataMatrix
        }
      }
    }
  }
  
  return(prob);
}
############################################################










