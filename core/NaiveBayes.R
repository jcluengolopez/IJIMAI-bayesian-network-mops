
# CARGAR DEPENDENCIAS AL PRINCIPIO (discretización)
library(discretization);


NaiveBayes = function(data,class,dataVal,dataMin,dataMax,varType,splitPoints,maxDegree,
                      nValidations,train,verbose=FALSE,plots=FALSE)
{
  # It makes the Naive Bayes model with or without cross validation. Depending on the class, 
  # it makes a classification or regression model
  # INPUT: - data = the database in a .RData (class is the last one)
  #        - class = the class without discretization in case of regression
  #        - dataVal = list with the different values of the discrete variables
  #        - dataMin = the minimum value in the variable
  #        - dataMax = the maximum value in the variable
  #        - varType = a vector with 0 if the variable is discrete, 1 if it's
  #       continuous and 2 if it has to be discretized. Program calculates it
  #        - splitPoints = the cut points in discretized variables (it can be NULL)
  #        - maxDegree = maximum degree for the polynomial fit
  #        - nValidations (optional) = number of cross validations
  #        - train (optional) = percentage of dataTrain when cross validation is not used
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - prob = a list with the structure of the model, values and probabilities 
  #         - testM = the percentage of well classified or the MSE (negative)
  
  
  # size of data
  dataDim = dim(data);
  
  # separate when the cross validation is made or not
  if (nValidations==0)          # no cross validation and no test
  {
    if (varType[dataDim[2]]==0)       # classification
    {
      # create the model
      prob = CreateModelNBClassification(data,dataVal,dataMin,dataMax,varType,splitPoints,
                                         maxDegree,verbose,plots);
      
      # remove the last row and column from the links matrix
      prob[[2]] = prob[[2]][-dataDim[2],-dataDim[2]];
      
      # no test and no errors
      testM = -1;          errors = -1;
    }
    else if (varType[dataDim[2]]==1)       # regression without discretization
    {
      # create the model
      prob = CreateModelNBRegression(data,dataVal,dataMin,dataMax,varType,splitPoints,
                                     maxDegree,verbose,plots);

      # remove the last row and column from the links matrix
      prob[[2]] = prob[[2]][-dataDim[2],-dataDim[2]];

      # no test and no errors
      testM = -1;          errors = -1;
    }
    else       # regression with discretization
    {
      # create a classification model
      prob = CreateModelNBClassification(data,dataVal,dataMin,dataMax,varType,splitPoints,
                                         maxDegree,verbose,plots);
      
      # remove the last row and column from the links matrix
      prob[[2]] = prob[[2]][-dataDim[2],-dataDim[2]];
      
      # no test and no errors
      testM = -1;          errors = -1;
    }
  }
  
  else if (nValidations==1)          # no cross validation, but test
  {
    # the size of dataTrain and dataTest is in train
    size = floor(dataDim[1]*train);
    
    # separate data train and data test
    set.seed(20);            selectTrain = sample(1:dataDim[1],size);
    dataTrain = data[selectTrain,];
    dataTest = data[-selectTrain,];
    
    if (varType[dataDim[2]]==0)       # classification
    {
      # create the model
      prob = CreateModelNBClassification(dataTrain,dataVal,dataMin,dataMax,varType,
                                          splitPoints,maxDegree,verbose,plots);
      
      # remove the last row and column from the links matrix
      prob[[2]] = prob[[2]][-dataDim[2],-dataDim[2]];
      
      # test the model and get the % of well classified and errors
      prediction = TestModelClassification(prob,dataTest,verbose,plots);
      testM = prediction[[1]];          errors = prediction[[2]];
    }
    else if (varType[dataDim[2]]==1)       # regression without discretization
    {
      # create the model
      prob = CreateModelNBRegression(dataTrain,dataVal,dataMin,dataMax,varType,splitPoints,
                                      maxDegree,verbose,plots);

      # remove the last row and column from the links matrix
      prob[[2]] = prob[[2]][-dataDim[2],-dataDim[2]];

      # test the model (negative errors:   higher error = better model)
      testM = -TestModelRegression(prob,dataTest,verbose,plots);
      
      # no matrix with errors in regression
      errors = 0;
    }
    else       # regression with discretization
    {
      # create a classification model
      prob = CreateModelNBClassification(dataTrain,dataVal,dataMin,dataMax,varType,
                                         splitPoints,maxDegree,verbose,plots);
      
      # remove the last row and column from the links matrix
      prob[[2]] = prob[[2]][-dataDim[2],-dataDim[2]];
      
      # predict the datatest as a classification model
      prediction = PredictDataClassification(prob,dataTest[,1:(ncol(dataTest)-1)]);
      
      # calculate the mid points of each interval of the discretization
      cutPoints = splitPoints[[length(splitPoints)]];
      midPoints = cutPoints[-length(cutPoints)] + diff(cutPoints) / 2;
      
      # make the continuous prediction
      predictionCont = midPoints[prediction];
      
      # select the test of the original class
      classTest = class[-selectTrain];
      
      # test the model (negative errors:   higher error = better model)
      testM = -mean((predictionCont - classTest) ^ 2);
      
      # no matrix with errors in regression
      errors = 0;
    }
  }
  
  else          # cross validation
  {
    # initialize the % of well classified observations in each validation and errors
    test = vector(mode='numeric',length=nValidations);
    errors = 0;
    
    # the vector with the cut points for the cross validation
    size = floor(dataDim[1] / nValidations);
    divisions = c(0,(1:(nValidations-1))*size,dataDim[1]);
    
    # for each cross validation
    for (i in 1:nValidations)
    {
      # separate the training data and the test data
      selectTest = (divisions[i]+1):divisions[i+1];
      dataTrain = data[-selectTest,];
      dataTest = data[selectTest,];
      
      # classification or regression
      if (varType[dataDim[2]]==0)       # classification
      {
        # create the model
        prob = CreateModelNBClassification(dataTrain,dataVal,dataMin,dataMax,varType,
                                            splitPoints,maxDegree,verbose,plots);
        
        # remove the last row and column from the links matrix
        prob[[2]] = prob[[2]][-dataDim[2],-dataDim[2]];
        
        # test the model and get the % of well classified and errors
        prediction = TestModelClassification(prob,dataTest,verbose,plots);
        test[i] = prediction[[1]];          errors = errors + prediction[[2]];
      }
      else if (varType[dataDim[2]]==1)       # regression without discretization
      {
        # create the model
        prob = CreateModelNBRegression(dataTrain,dataVal,dataMin,dataMax,varType,
                                        splitPoints,maxDegree,verbose,plots);

        # remove the last row and column from the links matrix
        prob[[2]] = prob[[2]][-dataDim[2],-dataDim[2]];

        # test the model (negative errors:   higher error = better model)
        test[i] = -TestModelRegression(prob,dataTest,verbose,plots);
      }
      else       # regression with discretization
      {
        # create a classification model
        prob = CreateModelNBClassification(dataTrain,dataVal,dataMin,dataMax,varType,
                                           splitPoints,maxDegree,verbose,plots);
        
        # remove the last row and column from the links matrix
        prob[[2]] = prob[[2]][-dataDim[2],-dataDim[2]];
        
        # predict the datatest as a classification model
        prediction = PredictDataClassification(prob,dataTest[,1:(ncol(dataTest)-1)]);
        
        # calculate the mid points of each interval of the discretization
        cutPoints = splitPoints[[length(splitPoints)]];
        midPoints = cutPoints[-length(cutPoints)] + diff(cutPoints) / 2;
        
        # make the continuous prediction
        predictionCont = midPoints[prediction];
        
        # select the test of the original class
        classTest = class[selectTest];
        
        # test the model (negative errors:   higher error = better model)
        test[i] = -mean((predictionCont - classTest) ^ 2);
        
        # no matrix with errors in regression
        errors = 0;
      }
    }
    
    # the mean of the percentage of correct classifications or the errors
    testM = mean(test);
  }
  
  # the final model is made using the whole dataset
  # classification or regression
  if (varType[dataDim[2]]==0)       # classification
  {
    # create the final model
    prob = CreateModelNBClassification(data,dataVal,dataMin,dataMax,varType,splitPoints,
                                        maxDegree,verbose,plots);
  }
  else if (varType[dataDim[2]]==1)       # regression without discretization
  {
    # create the model
    prob = CreateModelNBRegression(data,dataVal,dataMin,dataMax,varType,splitPoints,
                                    maxDegree,verbose,plots);
  }
  else       # regression with discretization
  {
    # create a classification model
    prob = CreateModelNBClassification(data,dataVal,dataMin,dataMax,varType,
                                       splitPoints,maxDegree,verbose,plots);
  }
  
  return(list(prob,testM,errors));
}
############################################################









