
# CARGAR DEPENDENCIAS AL PRINCIPIO 


LearnParameters = function(data,dataVal,dataMin,dataMax,varType,splitPoints,maxDegree,links,
                           fX,fXY,fXYZ,BICpartX,BICpartXY,BICpartXYZ,verbose=FALSE)
{
  # it obtains all the distribution of a given structure using some distributions that have
  # already been estimated in previous steps. It also computes the BIC using other partial BIC
  # INPUT: - data = a data-frame with the data of study
  #        - dataVal = list with the different values of the discrete variables
  #        - dataMin = the minimum value in the variable
  #        - dataMax = the maximum value in the variable
  #        - varType = a vector with 0 if the variable is discrete, 1 if it's
  #       continuous and 2 if it has to be discretized. Program calculates it
  #        - splitPoints = the cut points in discretized variables (it can be NULL)
  #        - maxDegree = maximum degree for the polynomial fit
  #        - links = matrix of the links to build the model
  #        - fX,fXY,fXYZ = previously obtained distributions
  #        - BICpartX,BICpartXY,BICpartXYZ = previously calculated partial BIC
  #        - verbose (optional) = show the 'cats' on screen
  # OUTPUT:  - model = a list with the structure of the model, values and probabilities
  #          - BIC = BIC of the new model
  #          - fX,fXY,fXYZ = uploaded distributions
  #          - BICpartX,BICpartXY,BICpartXYZ = uploaded partial BIC
  
  
  # sum of each row to know if the variable must be included in the model
  sumRow = apply(links,1,sum);
  sumCol = apply(links,2,sum);
  
  # set of variables to use in the model
  usedVar = c(which(sumRow>0),which(sumCol>0));
  usedVar = sort(unique(usedVar));
  
  # set of variables not used in the model
  aux = 1:length(dataVal);
  notUsedVar = aux[-usedVar];
  
  #cat('\n variables utilizadas = ',usedVar);
  
  # initialize the bic
  bic = 0;
  
  # number of variables of the model
  nVar = length(usedVar);
  
  # create the list for probabilities and polynomials (3 each variable, links and varType
  model = vector("list",ncol(data)*3+2);
  
  # put the structure of the model in prob
  model[[1]] = varType;            model[[2]] = links;
  
  # for each variable of the model
  for (i in usedVar)
  {
    # get the number of parents
    nParents = sumCol[i];
    
    #cat('\n número de padres =',nParents,' en variable ',i);
    
    
    if (nParents==0)      # marginal distribution are already estimated
    {
      model[[3*i]] = fX[[i]][[1]];
      model[[3*i+1]] = fX[[i]][[2]];
      
      # upload the BIC
      bic = bic + BICpartX[i];
      
      # # check if it was previously estimated and the BIC was previously calculated
      # if (length(fX[[i]])>0)   # done
      # {
      #   #cat('\n Sin padres. La variable ',i,' ya está estimada');
      #   
      #   
      #   model[[3*i]] = fX[[i]][[1]];
      #   model[[3*i+1]] = fX[[i]][[2]];
      # }
      # else       # not done
      # {
      #   # the fit depends on the type of variable
      #   if (varType[i]==0)   # when the variable is discrete
      #   {
      #     # the function that calculates the conditional probability and the values
      #     p = Adjust1ParentDD(data[,i],dataVal[[i]],rep(0,nrow(data)),0,verbose);
      #     
      #     # insert the results in the list
      #     model[[3*i]] = dataVal[[i]];
      #     model[[3*i+1]] = p[[1]];
      #     
      #     # the BIC is previously calculated
      #   }
      #   else if (varType[i]==1)   # when the variable is continuous
      #   {
      #     # calculate the coefficients of the polynomial and the tails
      #     p = Adjust1ParentCD(data[,i],dataMin[i],dataMax[i],rep(0,nrow(data)),0,i,
      #                         maxDegree,verbose);
      #     
      #     # insert the results in the list
      #     model[[3*i]] = p[[1]];
      #     model[[3*i+1]] = p[[2]];
      #     
      #     # the BIC is previously calculated
      #   }
      #   else   # when the variable is discretized
      #   {
      #     # the function that calculates the conditional probability and the values
      #     p = Adjust1ParentDD(data[,i],dataVal[[i]],rep(0,nrow(data)),0,verbose);
      #     
      #     # insert the cut points and the probabilities in the list prob
      #     model[[3*i]] = splitPoints[[i]];
      #     model[[3*i+1]] = p[[1]];
      #     
      #     # the BIC is previously calculated
      #   }
      #   
      #   # upload the database of distributions
      #   fX[[i]] = vector(mode='list',length=2);
      #   fX[[i]][[1]] = model[[3*i]];
      #   fX[[i]][[2]] = model[[3*i+1]];
      # }
      # 
      # # upload the BIC
      # bic = bic + BICpartX[i];
    }
    
    
    
    
    else if (nParents==1)      # 1 parent distribution
    {
      # obtain the position of the parent
      posParent = which(links[,i]==1);
      
      # obtain the position in the list of distributions
      posfXY = i + (posParent - 1) * ncol(data);
      
      # check if it was previously estimated and the BIC was previously calculated
      if (length(fXY[[posfXY]])>0)   # done
      {
        #cat('\n Con un padre. La variable ',i,' ya está estimada');
        
        # the size of the model depends on the type of parent
        if (varType[posParent]==1)   # 1 continuous parent
        {
          model[[3*i]] = fXY[[posfXY]][[1]];
          model[[3*i+1]] = fXY[[posfXY]][[2]];
          model[[3*i+2]] = fXY[[posfXY]][[3]];
        }
        else   # 1 discrete or discretized parent
        {
          model[[3*i]] = fXY[[posfXY]][[1]];
          model[[3*i+1]] = fXY[[posfXY]][[2]];
        }
      }
      else       # not done
      {
        # the model depends on the type of parent
        if (varType[posParent]==1)   # 1 continuous parent
        {
          # the fit depends on the type of variable
          if (varType[i]==0)   # when the variable is discrete
          {
            # calculate the conditional probabilities and the values
            p = Adjust1ParentDC(data[,i],dataVal[[i]],data[,posParent],dataMin[posParent],
                                dataMax[posParent],verbose);
            
            # insert the results in the list
            model[[3*i]] = dataVal[[i]];
            model[[3*i+1]] = p[[2]];    # cut points
            model[[3*i+2]] = p[[1]];    # dataMatrix
            
            # upload the BIC
            #BICpartXY[posfXY] = p[[3]];
            BICpartXY[posfXY] = BIC1ParentC(p[[1]],data[,i],dataVal[[i]],data[,posParent],
                                            p[[2]],varType[i]);
          }
          else if (varType[i]==1)   # when the variable is continuous
          {
            # calculate the polynomials conditioned to the class
            p = Adjust1ParentCC(data[,i],dataMin[i],dataMax[i],data[,posParent],
                                dataMin[posParent],dataMax[posParent],i,maxDegree,verbose);
            
            # insert the results in the list
            model[[3*i]] = p[[1]];       # values
            model[[3*i+1]] = p[[3]];     # cutPoints
            model[[3*i+2]] = p[[2]];     # polynomial
            
            # upload the BIC
            #BICpartXY[posfXY] = p[[4]];
            BICpartXY[posfXY] = BIC1ParentC(p[[2]],data[,i],p[[1]],data[,posParent],
                                            p[[3]],varType[i]);
          }
          else   # when the variable is discretized
          {
            # calculate the conditional probabilities and the values
            p = Adjust1ParentDC(data[,i],dataVal[[i]],data[,posParent],dataMin[posParent],
                                dataMax[posParent],verbose);
            
            # insert the cut points in the list prob. When there're no cut points
            model[[3*i]] = splitPoints[[i]];
            
            # insert the results in the list
            model[[3*i+1]] = p[[2]];    # cut points
            model[[3*i+2]] = p[[1]];    # dataMatrix
            
            # upload the BIC
            #BICpartXY[posfXY] = p[[3]];
            BICpartXY[posfXY] = BIC1ParentC(p[[1]],data[,i],dataVal[[i]],data[,posParent],
                                            p[[2]],varType[i]);
          }
          
          # upload the database of distributions
          fXY[[posfXY]] = vector(mode='list',length=3);
          fXY[[posfXY]][[1]] = model[[3*i]];
          fXY[[posfXY]][[2]] = model[[3*i+1]];
          fXY[[posfXY]][[3]] = model[[3*i+2]];
          
          # the continuous parent is discretized with values 1,2,3...
          parentVal = 1:(length(model[[3*i+1]]) - 1);       
          
          # upload the BIC
          #BICpartXY[posfXY] = BIC1Parent(model[[3*i+2]],data[,i],model[[3*i]],
          #                               data[,posParent],parentVal,varType[i],verbose);
        }
        else   # 1 discrete or discretized parent
        {
          # the fit depends on the type of variable
          if (varType[i]==0)   # when the variable is discrete
          {
            # the function that calculates the conditional probability and the values
            p = Adjust1ParentDD(data[,i],dataVal[[i]],data[,posParent],dataVal[[posParent]],
                                verbose);
            
            # insert the results in the list
            model[[3*i]] = dataVal[[i]];
            model[[3*i+1]] = p[[1]];
            
            # upload the BIC
            #BICpartXY[posfXY] = p[[2]];
            BICpartXY[posfXY] = BIC1ParentD(p[[1]],data[,i],dataVal[[i]],data[,posParent],
                                            dataVal[[posParent]],varType[i]);
          }
          else if (varType[i]==1)   # when the variable is continuous
          {
            # calculate the polynomials conditioned to the class
            p = Adjust1ParentCD(data[,i],dataMin[i],dataMax[i],data[,posParent],
                                dataVal[[posParent]],i,maxDegree,verbose);
            
            # insert the results in the list
            model[[3*i]] = p[[1]];
            model[[3*i+1]] = p[[2]];
            
            # upload the BIC
            #BICpartXY[posfXY] = p[[4]];
            BICpartXY[posfXY] = BIC1ParentD(p[[2]],data[,i],p[[1]],data[,posParent],
                                            dataVal[[posParent]],varType[i]);
          }
          else   # when the variable is discretized
          {
            # calculate the conditional probabilities and the values
            p = Adjust1ParentDD(data[,i],dataVal[[i]],data[,posParent],
                                dataVal[[posParent]],verbose);
            
            # insert the cut points in the list prob. When there're no cut points
            model[[3*i]] = splitPoints[[i]];
            
            # insert the probabilities in the list prob
            model[[3*i+1]] = p[[1]];
            
            # upload the BIC
            #BICpartXY[posfXY] = p[[2]];
            BICpartXY[posfXY] = BIC1ParentD(p[[1]],data[,i],dataVal[[i]],data[,posParent],
                                            dataVal[[posParent]],varType[i]);
          }
          
          # upload the database of distributions
          fXY[[posfXY]] = vector(mode='list',length=2);
          fXY[[posfXY]][[1]] = model[[3*i]];
          fXY[[posfXY]][[2]] = model[[3*i+1]];
          
          # upload the BIC
          #cat('\n     Calcula el BIC con 1 padre discreto');
          
          
          # upload the BIC
          #BICpartXY[posfXY] = BIC1Parent(model[[3*i+1]],data[,i],model[[3*i]],
          #                               data[,posParent],model[[3*posParent]],varType[i]);
        }
      }
      
      # upload the BIC
      bic = bic + BICpartXY[posfXY];
    }
    
    
    
    
    else      # 2 parents distribution
    {
      # obtain the position of both parents
      posParent1 = which(links[,i]==1)[1];
      posParent2 = which(links[,i]==1)[2];
      
      # obtain the position in the list of distributions
      posfXYZ = i + (posParent1 - 1) * ncol(data) + (posParent2 - 1) * ncol(data)^2;
      
      # check if it was previously estimated and the BIC was previously calculated
      if (length(fXYZ[[posfXYZ]])>0)   # done
      {
        #cat('\n Con dos padres. La variable ',i,' ya está estimada');
        
        # the size of the model depends on the type of parent
        if (varType[posParent1]!=1 & varType[posParent2]!=1)   # 2 non-continuous parents
        {
          model[[3*i]] = fXYZ[[posfXYZ]][[1]];
          model[[3*i+1]] = fXYZ[[posfXYZ]][[2]];
        }
        else   # 2 parents including one continuous
        {
          model[[3*i]] = fXYZ[[posfXYZ]][[1]];
          model[[3*i+1]] = fXYZ[[posfXYZ]][[2]];
          model[[3*i+2]] = fXYZ[[posfXYZ]][[3]];
        }
      }
      
      
      else       # not done
      {
        # the model depends on the type of parents
        if (varType[posParent1]!=1 & varType[posParent2]!=1)   # 2 non continuous parents
        {
          # the fit depends on the type of variable
          if (varType[i]==0)   # when the variable is discrete
          {
            # calculate the conditional probabilities and the values
            p = Adjust2ParentsDDD(data[,i],dataVal[[i]],data[,posParent1],
                        dataVal[[posParent1]],data[,posParent2],dataVal[[posParent2]]);
            
            # insert the results in the list
            model[[3*i]] = dataVal[[i]];
            model[[3*i+1]] = p[[1]];              # dataMatrix
            
            # upload the BIC
            #BICpartXYZ[posfXYZ] = p[[2]];
            BICpartXYZ[posfXYZ] = BIC2ParentsDD(p[[1]],data[,i],dataVal[[i]],
                                    data[,posParent1],dataVal[[posParent1]],
                                    data[,posParent2],dataVal[[posParent2]],varType[i]);
          }
          else if (varType[i]==1)   # when the variable is continuous
          {
            # calculate the polynomials conditioned to the class
            p = Adjust2ParentsCDD(data[,i],dataMin[i],dataMax[i],data[,posParent1],
                                  dataVal[[posParent1]],data[,posParent1],
                                  dataVal[[posParent2]],i,maxDegree,verbose);
            
            # insert the results in the list
            model[[3*i]] = p[[1]];       # values
            model[[3*i+1]] = p[[2]];     # polynomial
            
            # upload the BIC
            #BICpartXYZ[posfXYZ] = p[[4]];
            BICpartXYZ[posfXYZ] = BIC2ParentsDD(p[[2]],data[,i],p[[1]],
                                     data[,posParent1],dataVal[[posParent1]],
                                     data[,posParent2],dataVal[[posParent2]],varType[i]);
          }
          else   # when the variable is discretized
          {
            # calculate the conditional probabilities and the values
            p = Adjust2ParentsDDD(data[,i],dataVal[[i]],data[,posParent1],
                      dataVal[[posParent1]],data[,posParent2],dataVal[[posParent2]]);
            
            # insert the cut points in the list prob
            model[[3*i]] = splitPoints[[i]];
            
            # insert the probabilities in the list prob
            model[[3*i+1]] = p[[1]];
            
            # upload the BIC
            #BICpartXYZ[posfXYZ] = p[[2]];
            BICpartXYZ[posfXYZ] = BIC2ParentsDD(p[[1]],data[,i],dataVal[[i]],
                                      data[,posParent1],dataVal[[posParent1]],
                                      data[,posParent2],dataVal[[posParent2]],varType[i]);
          }
          
          # upload the database of distributions
          fXYZ[[posfXYZ]] = vector(mode='list',length=2);
          fXYZ[[posfXYZ]][[1]] = model[[3*i]];
          fXYZ[[posfXYZ]][[2]] = model[[3*i+1]];
        }
        
        else if (varType[posParent1]==1 & varType[posParent2]==1)   # 2 continuous parents
        {
          # the fit depends on the type of variable
          if (varType[i]==0)   # when the variable is discrete
          {
            # calculate the conditional probabilities and the values
            p = Adjust2ParentsDCC(data[,i],dataVal[[i]],data[,posParent1],dataMin[posParent1],
                                  dataMax[posParent1],data[,posParent2],dataMin[posParent2],
                                  dataMax[posParent2],i,maxDegree,verbose);
            
            # insert the results in the list
            model[[3*i]] = dataVal[[i]];
            model[[3*i+1]] = p[[2]];    # cut points
            model[[3*i+2]] = p[[1]];    # dataMatrix
            
            # upload the BIC
            #BICpartXYZ[posfXYZ] = p[[3]];
            BICpartXYZ[posfXYZ] = BIC2ParentsCC(p[[1]],data[,i],dataVal[[i]],
                                    data[,posParent1],data[,posParent2],p[[2]],varType[i]);
          }
          else if (varType[i]==1)   # when the variable is continuous
          {
            # calculate the polynomials conditioned to both parents
            p = Adjust2ParentsCCC(data[,i],dataMin[i],dataMax[i],data[,posParent1],
                                  dataMin[posParent1],dataMax[posParent1],data[,posParent2],
                                  dataMin[posParent2],dataMax[posParent2],i,maxDegree);
            
            # insert the results in the list
            model[[3*i]] = p[[1]];       # values
            model[[3*i+1]] = p[[3]];     # cutPoints
            model[[3*i+2]] = p[[2]];     # polynomial
            
            # upload the BIC
            #BICpartXYZ[posfXYZ] = p[[4]];
            BICpartXYZ[posfXYZ] = BIC2ParentsCC(p[[2]],data[,i],p[[1]],data[,posParent1],
                                        data[,posParent2],p[[3]],varType[i]);
          }
          else   # when the variable is discretized
          {
            # calculate the conditional probabilities and the values
            p = Adjust2ParentsDCC(data[,i],dataVal[[i]],data[,posParent1],dataMin[posParent1],
                                  dataMax[posParent1],data[,posParent2],dataMin[posParent2],
                                  dataMax[posParent2],i,maxDegree,verbose);
            
            # insert the cut points in the list prob. When there're no cut points
            model[[3*i]] = splitPoints[[i]];
            
            # insert the results in the list
            model[[3*i+1]] = p[[2]];    # cut points
            model[[3*i+2]] = p[[1]];    # dataMatrix
            
            # upload the BIC
            #BICpartXYZ[posfXYZ] = p[[3]];
            BICpartXYZ[posfXYZ] = BIC2ParentsCC(p[[1]],data[,i],dataVal[[i]],
                                    data[,posParent1],data[,posParent2],p[[2]],varType[i]);
          }
          
          # upload the database of distributions
          fXYZ[[posfXYZ]] = vector(mode='list',length=3);
          fXYZ[[posfXYZ]][[1]] = model[[3*i]];
          fXYZ[[posfXYZ]][[2]] = model[[3*i+1]];
          fXYZ[[posfXYZ]][[3]] = model[[3*i+2]];
          
          #cat('\n     Calcula el BIC con 2 padres continuos');
        }
        
        else if (varType[posParent1]==1 & varType[posParent2]!=1)   # 2 parents (cont-disc)
        {
          # the fit depends on the type of variable
          if (varType[i]==0)   # when the variable is discrete
          {
            # the function that calculates the conditional probability and the values
            p = Adjust2ParentsDCD(data[,i],dataVal[[i]],data[,posParent1],dataMin[posParent1],
                        dataMax[posParent1],data[,posParent2],dataVal[[posParent2]]);
            
            # insert the results in the list
            model[[3*i]] = dataVal[[i]];
            model[[3*i+1]] = p[[2]];    # cut points
            model[[3*i+2]] = p[[1]];    # dataMatrix
            
            # upload the BIC
            #BICpartXYZ[posfXYZ] = p[[3]];
            BICpartXYZ[posfXYZ] = BIC2ParentsCD(p[[1]],data[,i],dataVal[[i]],
                                    data[,posParent1],p[[2]],data[,posParent2],
                                    dataVal[[posParent2]],varType[i]);
          }
          else if (varType[i]==1)   # when the variable is continuous
          {
            # calculate the polynomials conditioned to both variables
            p = Adjust2ParentsCCD(data[,i],dataMin[i],dataMax[i],data[,posParent1],
                                  dataMin[posParent1],dataMax[posParent1],data[,posParent2],
                                  dataVal[[posParent2]],i,maxDegree,verbose);
            
            # insert the results in the list
            model[[3*i]] = p[[1]];       # values
            model[[3*i+1]] = p[[3]];     # cutPoints
            model[[3*i+2]] = p[[2]];     # polynomial
            
            # upload the BIC
            #BICpartXYZ[posfXYZ] = p[[4]];
            BICpartXYZ[posfXYZ] = BIC2ParentsCD(p[[2]],data[,i],p[[1]],
                                      data[,posParent1],p[[3]],data[,posParent2],
                                      dataVal[[posParent2]],varType[i]);
          }
          else   # when the variable is discretized
          {
            # the function that calculates the conditional probability and the values
            p = Adjust2ParentsDCD(data[,i],dataVal[[i]],data[,posParent1],dataMin[posParent1],
                        dataMax[posParent1],data[,posParent2],dataVal[[posParent2]]);
            
            # insert the cut points in the list prob. When there're no cut points
            model[[3*i]] = splitPoints[[i]];
            
            # insert the probabilities and the cut points in the list prob
            model[[3*i+1]] = p[[2]];
            model[[3*i+2]] = p[[1]];
            
            # upload the BIC
            #BICpartXYZ[posfXYZ] = p[[3]];
            BICpartXYZ[posfXYZ] = BIC2ParentsCD(p[[1]],data[,i],dataVal[[i]],
                                      data[,posParent1],p[[2]],data[,posParent2],
                                      dataVal[[posParent2]],varType[i]);
          }
          
          # upload the database of distributions
          fXYZ[[posfXYZ]] = vector(mode='list',length=3);
          fXYZ[[posfXYZ]][[1]] = model[[3*i]];
          fXYZ[[posfXYZ]][[2]] = model[[3*i+1]];
          fXYZ[[posfXYZ]][[3]] = model[[3*i+2]];
        }
        
        else   # 2 parents (disc-cont)
        {
          # the fit depends on the type of variable
          if (varType[i]==0)   # when the variable is discrete
          {
            # the function that calculates the conditional probability and the values
            p = Adjust2ParentsDCD(data[,i],dataVal[[i]],data[,posParent2],dataMin[posParent2],
                        dataMax[posParent2],data[,posParent1],dataVal[[posParent1]]);
            
            # insert the results in the list
            model[[3*i]] = dataVal[[i]];
            model[[3*i+1]] = p[[2]];    # cut points
            model[[3*i+2]] = p[[1]];    # dataMatrix
            
            # upload the BIC
            #BICpartXYZ[posfXYZ] = p[[3]];
            BICpartXYZ[posfXYZ] = BIC2ParentsCD(p[[1]],data[,i],dataVal[[i]],
                                      data[,posParent2],p[[2]],data[,posParent1],
                                      dataVal[[posParent1]],varType[i]);
          }
          else if (varType[i]==1)   # when the variable is continuous
          {
            # calculate the polynomials conditioned to both variables
            p = Adjust2ParentsCCD(data[,i],dataMin[i],dataMax[i],data[,posParent2],
                                  dataMin[posParent2],dataMax[posParent2],data[,posParent1],
                                  dataVal[[posParent1]],i,maxDegree,verbose);
            
            # insert the results in the list
            model[[3*i]] = p[[1]];       # values
            model[[3*i+1]] = p[[3]];     # cutPoints
            model[[3*i+2]] = p[[2]];     # polynomial
            
            # upload the BIC
            #BICpartXYZ[posfXYZ] = p[[4]];
            BICpartXYZ[posfXYZ] = BIC2ParentsCD(p[[2]],data[,i],p[[1]],
                                      data[,posParent2],p[[3]],data[,posParent1],
                                      dataVal[[posParent1]],varType[i]);
          }
          else   # when the variable is discretized
          {
            # the function that calculates the conditional probability and the values
            p = Adjust2ParentsDCD(data[,i],dataVal[[i]],data[,posParent2],dataMin[posParent2],
                        dataMax[posParent2],data[,posParent1],dataVal[[posParent1]]);
            
            # insert the cut points in the list prob. When there're no cut points
            model[[3*i]] = splitPoints[[i]];
            
            # insert the probabilities and the cut points in the list prob
            model[[3*i+1]] = p[[2]];
            model[[3*i+2]] = p[[1]];
            
            # upload the BIC
            #BICpartXYZ[posfXYZ] = p[[3]];
            BICpartXYZ[posfXYZ] = BIC2ParentsCD(p[[1]],data[,i],dataVal[[i]],
                                      data[,posParent2],p[[2]],data[,posParent1],
                                      dataVal[[posParent1]],varType[i]);
          }
          
          # upload the database of distributions
          fXYZ[[posfXYZ]] = vector(mode='list',length=3);
          fXYZ[[posfXYZ]][[1]] = model[[3*i]];
          fXYZ[[posfXYZ]][[2]] = model[[3*i+1]];
          fXYZ[[posfXYZ]][[3]] = model[[3*i+2]];
        }
      }
      
      # upload the BIC
      bic = bic + BICpartXYZ[posfXYZ];
    }
  }
  
  # update the BIC with the non used variables
  bic = bic + sum(BICpartX[notUsedVar]);
  
  cat('\n variables usadas = ',usedVar);
  cat('\n variables no usadas = ',notUsedVar);
  cat('\n BIC = ',bic);
  
  return(list(model,fX,fXY,fXYZ,BICpartX,BICpartXY,BICpartXYZ,bic))
}
############################################################









FitMarginalDistributiuons = function(data,dataVal,dataMin,dataMax,varType,splitPoints,
                                     maxDegree,verbose=FALSE)
{
  # it estimates the marginal distribution and BIC of each variable of the model
  # INPUT: - data = a data-frame with the data of study
  #        - dataVal = list with the different values of the discrete variables
  #        - dataMin = the minimum value in the variable
  #        - dataMax = the maximum value in the variable
  #        - varType = a vector with 0 if the variable is discrete, 1 if it's
  #       continuous and 2 if it has to be discretized. Program calculates it
  #        - splitPoints = the cut points in discretized variables (it can be NULL)
  #        - maxDegree = maximum degree for the polynomial fit
  #        - verbose (optional) = show the 'cats' on screen
  # OUTPUT:  - fX = marginal distribution
  #          - bic = BIC of each distribution
  
  
  # size of the dataset
  N = nrow(data) * ncol(data);
  
  # create the list for the distributions and a vector for the BIC
  fX = vector(mode='list',length=length(varType));
  bic = rep(NA,length=length(varType));
  
  for (i in 1:length(varType))
  {
    # two positions for each distribution
    fX[[i]] = vector(mode='list',length=2);
    
    # the fit depends on the type of variable
    if (varType[i]==0)   # when the variable is discrete
    {
      # the function that calculates the conditional probability and the values
      p = Adjust1ParentDD(data[,i],dataVal[[i]],rep(0,nrow(data)),varType[i],verbose);
      
      # upload the database of distributions
      fX[[i]][[1]] = dataVal[[i]];
      fX[[i]][[2]] = p[[1]];
      
      # upload the BIC of this distribution
      #bic[i] = p[[2]];
      bic[i] = BIC0Parent(p[[1]],data[,i],dataVal[[i]],varType[i]);
    }
    else if (varType[i]==1)   # when the variable is continuous
    {
      # calculate the coefficients of the polynomial and the tails
      p = Adjust1ParentCD(data[,i],dataMin[i],dataMax[i],rep(0,nrow(data)),0,i,
                          maxDegree,verbose);
      
      # upload the database of distributions
      fX[[i]][[1]] = p[[1]];
      fX[[i]][[2]] = p[[2]];
      
      # upload the BIC of this distribution
      #bic[i] = p[[4]];
      bic[i] = BIC0Parent(p[[2]],data[,i],p[[1]],varType[i]);
    }
    else   # when the variable is discretized
    {
      # the function that calculates the conditional probability and the values
      p = Adjust1ParentDD(data[,i],dataVal[[i]],rep(0,nrow(data)),0,verbose);
      
      # upload the database of distributions
      fX[[i]][[1]] = splitPoints[[i]];
      fX[[i]][[2]] = p[[1]];
      
      # upload the BIC of this distribution
      #bic[i] = p[[2]];
      bic[i] = BIC0Parent(p[[1]],data[,i],dataVal[[i]],0);
    }
  }
  
  return(list(fX,bic))
}
############################################################

















