
# CARGAR DEPENDENCIAS AL PRINCIPIO
library(classInt);

NBTAN = function(data,varType,method='TAN',discretize=0,discretizeType=0,selectVar=TRUE,
                 train=0.8,nValidations=5,maxDegree=10,verbose=FALSE,plots=FALSE)
{
  # It makes the Naive Bayes or TAN model starting from a database composed by a set 
  # of discrete and continuous variables. Depending on the class, it makes a 
  # classification or regression model with or without selection of variables and cross
  # validation
  # INPUT: - data = the database in a .RData (class is the last one)
  #        - varType (optional) = a vector with 0 if the variable is discrete, 1
  #       if it's continuous and 2 if it has to be discretized. Program calculates it
  #        - method (optional) = a string with the method ('TAN' or 'NaiveBayes')
  #        - discretize (optional) = it decides when use a discretization. 
  #       0=never, 1=between 20 and 5% of the values, 2=more than 20 values
  #        - discretizeType (optional) = it decides the type of discretization. 0=dynamic,
  #       1=equal width
  #        - selectVar (optional) = it decides when use the selection of variables
  #        - train (optional) = percentage of dataTrain when cross validation is not used
  #        - nValidations (optional) = number of cross validations
  #        - maxDegree (optional) = maximum degree for the polynomial fit
  #        - verbose (optional) = show the 'cats' on screen
  #        - plots (optional) = draw the graphs
  # OUTPUT:  - prob = a list with the structure of the model, values and probabilities 
  #          - test = the percentage of well classified test data on average
  #          - errors = confusion matrix for classification models
  
  
  # number of variables
  nVar = ncol(data);
  
  # if there's not varType, calculate it
  if (missing(varType))   {   varType = DiscOrCont(data,discretize);   }
  
  # define the kind of model
  if (varType[nVar]==0)   {   model = 'classification';   }
  else   {   model = 'regression';   }
  
  # calculate the min and max of all the variables (NAs introduced)
  dataMin = as.numeric(sapply(data,min)) - 0.000000001;
  dataMax = as.numeric(sapply(data,max)) + 0.000000001;
  
  # shuffle the data set
  set.seed(20)
  data = data[sample(1:nrow(data)),];
  
  # keep the class before the discretization
  class = data[,nVar];
  
  # check if there are variables to discretize
  if (length(which(varType==2))>0)
  {
    # do the discretization of the data frame 
    discretized = Discretization(data,varType,discretizeType,dataMin,dataMax,model);
    
    # extract the information of the discretization
    data = discretized[[1]];               splitPoints = discretized[[2]];
  }
  
  
  cat('\n varType = ',varType);
  
  # dataVal of discrete variables (it must be a list)
  dataVal = vector(mode='list',length=nVar);
  
  for (i in 1:nVar)
  {
    if (varType[i]==0)       # discrete
    {
      dataVal[[i]] = as.data.frame(table(data[,i]),stringsAsFactors=FALSE)[,1];
    }
    else       # discretized
    {
      dataVal[[i]] = 1:max(data[,i]);
    }
  }
  
  # check if the selection of variables has to be done
  if (selectVar==FALSE)       # no selection of variables
  {
    if (method=='TAN')     # TAN model
    {
      # calculate the mutual information conditioned to the class (only for TAN)
      MIC = MutualInfoClass(data,varType,verbose);
      
      # TAN model with all the variables
      TAN = TAN(data,class,dataVal,dataMin,dataMax,varType,splitPoints,maxDegree,nValidations,
                train,MIC,verbose,plots);
      
      # model and the well classified - error
      prob = TAN[[1]];          test = TAN[[2]];          errors = TAN[[3]];
    }
    else     # Naive Bayes model
    {
      # Naive Bayes model with all the variables
      NB = NaiveBayes(data,class,dataVal,dataMin,dataMax,varType,splitPoints,maxDegree,
                      nValidations,train,verbose,plots);
      
      # model and the well classified - error
      prob = NB[[1]];          test = NB[[2]];          errors = NB[[3]];
    }
    
    # reorganize the structure of the model to return
    prob2 = StructureModel(prob,names(data),varType,discretize,discretizeType,method);
  }
  else     # include the selection of variables
  {
    # calculate the mutual information of each variable and the class
    MI = MutualInfo(data,varType);
    
    # sort the variables
    rankingMI = sort(MI,decreasing=TRUE,index.return=TRUE);
    
    cat('\n rankingMI = ',rankingMI$ix);
    
    if (verbose)   {   cat('\n ranking de variables = ',rankingInfo$ix);   }
    
    # include only the class and another variable to start the process
    select = c(rankingMI$ix[1],nVar);
    
    # the process depends on the type of model
    if (method=='TAN')     # TAN model
    {
      # calculate the mutual information conditioned to the class (only for TAN)
      MIC = MutualInfoClass(data,varType,verbose);
      
      # make the ranking when there are more than two variables
      if (length(MI)<3)   {   rankingMIC = rankingMI$ix[2];   }
      else
      {
        # sort the variables according to its MIC with the variables included in the model
        rankingMIC = SortVariablesMIC(MIC,select[-length(select)],verbose);
      }
      
      # initialize the % of well-classified or error to 0 (always in the loop)
      test = -1e+30;
      
      # initialize the position of the variable in the ranking of MIC
      n = 1;
      
      # add new variables to the selection of variables process
      while (length(select)<nVar & n<=length(rankingMIC))
      {
        # include the new variable
        newSelect = sort(c(rankingMIC[n],select));
        newSelectMIC = newSelect[-length(newSelect)];
        
        # TAN model with the selected variables
        TAN = TAN(data[,newSelect],class,dataVal[newSelect],dataMin[newSelect],
                  dataMax[newSelect],varType[newSelect],splitPoints[newSelect],maxDegree,
                  nValidations,train,MIC[newSelectMIC,newSelectMIC],verbose,plots);
        
        # model, well classified - error and table of errors
        newProb = TAN[[1]];          newTest = TAN[[2]];          newErrors = TAN[[3]];
        
        # check if the results are better and change the parameters
        if (newTest > test)
        {
          prob = newProb;              test = newTest;
          select = newSelect;          errors = newErrors;
          
          # sort the variables according to its MIC with the variables included in the model
          rankingMIC = SortVariablesMIC(MIC,select[-length(select)],verbose);
          
          # start with the first variable
          n = 1;
        }
        else   {   n = n + 1;   }     # try using a new variable
      }
    }
    
    else      # Naive Bayes model
    {
      # initialize the % of well-classified or error to 0 (always in the loop)
      test = -1e+30;
      
      # delete the first variable in the ranking of MI
      rankingMI$ix = rankingMI$ix[-1];
      
      # initialize the position of the variable in the ranking of MI
      n = 1;
      
      # add new variables to the selection of variables process
      while (length(select)<nVar & n<=length(rankingMI$ix))
      {
        # include the new variable
        newSelect = sort(c(rankingMI$ix[n],select));
        
        # Naive Bayes model with the selected variables
        NB = NaiveBayes(data[,newSelect],class,dataVal[newSelect],dataMin[newSelect],
                        dataMax[newSelect],varType[newSelect],splitPoints[newSelect],
                        maxDegree,nValidations,train,verbose,plots);
        
        # model and well classified - error
        newProb = NB[[1]];          newTest = NB[[2]];          newErrors = NB[[3]];
        
        # check if the results are better and change the parameters
        if (newTest>test)
        {
          prob = newProb;              test = newTest;
          select = newSelect;          errors = newErrors;
          
          # delete the used variable from the ranking of MI
          rankingMI$ix = rankingMI$ix[-n];
          
          # the first variable to include is the size of select (it includes the class)
          n = 1;
        }
        else   {   n = n + 1;   }
      }
    }
    
    cat('\n select = ',select);
    
    # reorganize the structure of the model to return
    prob2 = StructureModel(prob,names(data)[select],varType[select],discretize,
                           discretizeType,method);
  }
  
  # make the opposite of the error in regression models
  if (varType[nVar]>0)   {   test = -test;   }
  
  return(list(prob2,test,errors));
}
############################################################










