
# CARGAR DEPENDENCIAS AL PRINCIPIO
library(classInt);
library(igraph);

GeneralBN = function(data,varType,discretize=0,discretizeType=0,train=0.8,nValidations=5,
                     maxDegree=7,verbose=FALSE,plots=FALSE)
{
  # It makes the structure of a general BN restricted to two parents from a database
  # composed by a set of discrete and continuous variables
  # INPUT: - data = the database in a .RData (class is the last one)
  #        - varType (optional) = a vector with 0 if the variable is discrete, 1
  #       if it's continuous and 2 if it has to be discretized. Program calculates it
  #        - discretize (optional) = it decides when use a discretization. 
  #       0=never, 1=between 20 and 5% of the values, 2=more than 20 values
  #        - discretizeType (optional) = it decides the type of discretization. 0=dynamic,
  #       1=equal width
  #        - train (optional) = percentage of dataTrain when cross validation is not used
  #        - nValidations (optional) = number of cross validations
  #        - maxDegree (optional) = maximum degree for the polynomial fit
  #        - verbose (optional) = show the 'cats' on screen
  #        - plots (optional) = draw the graphs
  # OUTPUT: - prob: list with 3 elements per variable, the links, the type of variables
  #         and the information of the class
  #         - bic = BIC of the final model
  
  
  # number of variables
  nVar = ncol(data);
  
  # if there's not varType, calculate it
  if (missing(varType))   {   varType = DiscOrCont(data,discretize);   }
  
  # calculate the min and max of all the variables (NAs introduced)
  dataMin = as.numeric(sapply(data,min)) - 0.000000001;
  dataMax = as.numeric(sapply(data,max)) + 0.000000001;
  
  # check if there are variables to discretize
  if (length(which(varType==2))>0)
  {
    # do the discretization of the data frame 
    discretized = Discretization(data,varType,discretizeType,dataMin,dataMax);
    
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
    else if (varType[i]==2)       # discretized
    {
      dataVal[[i]] = 1:max(data[,i]);
    }
  }
  
  # separate dataTrain and dataTest
  size = floor(nrow(data)*train);
  
  # separate data train and data test randomly
  set.seed(20);           selectTrain = sample(1:nrow(data),size);
  dataTrain = data[selectTrain,];
  dataTest = data[-selectTrain,];
  
  
  
  
  # EL PROGRAMA CALCULA LA INFORMACIÓN MUTUA Y LA ORDENA AUNQUE POR AHORA NO HAGA FALTA
  
  # calculate the mutual information
  #MIs = MutualInfoGeneral(data,varType);
  #MI = MIs[[1]];     MIC = MIs[[2]];
  
  # order the mutual info by decreasing and removing zeroes (not necessary)
  #order = sort.int(t(MI),decreasing=TRUE,index.return=TRUE)$ix;
  
  # select only the values with mutual information >0
  #n = (nVar^2 - nVar) / 2;   # number of elements under the diagonal
  #order = order[1:n];
  
  # get the position in the matrix of each mutual information
  #nRow = ceiling(order/nVar);
  #nCol = order - (nRow-1) * nVar;
  
  
  
  
  # list fo the marginal distributions
  marginal = FitMarginalDistributiuons(dataTrain,dataVal,dataMin,dataMax,varType,splitPoints,
                                 maxDegree);
  fX = marginal[[1]];
  
  # list for the distributions conditioned to one parent (rows are the parents)
  fXY = vector(mode='list',length=nVar^2);
  
  # list for the distributions conditioned to two parents
  fXYZ = vector(mode='list',length=nVar^3);
  
  
  # create the links matrix and initialize the relations to 0 (row is parent of column)
  links = matrix(0,nrow=nVar,ncol=nVar);
  
  # vector for the diagonal
  diagLinks = (1:nVar - 1) * nVar + 1:nVar;
  
  # create the BIC and the partial BICs for each distribution
  BICpartX = marginal[[2]];
  BICpartXY = rep(NA,length=nVar^2);
  BICpartXYZ = rep(NA,length=nVar^3);
  BIC = -999999;
  
  
  continue = TRUE;      nChanges = 0;   # to know how many changes have been made
  
  # repeat the process while the BIC increases
  while (continue==TRUE)
  {
    cat('\n \n \n  Iteración  ',nChanges+1);
    
    # OPTION 1: add a new link
    
    
    cat('\n \n           OPCIÓN 1');
    
    
    # select the links that are not previously included
    notLinks = which(t(links)==0);
    
    # delete the diagonal
    toInclude = notLinks[! notLinks %in% diagLinks];
    
    # initialize a list to save the matrices of links and the BIC
    newLinks1 = vector(mode='list',length=length(toInclude));
    newBIC1 = rep(NA,length=length(toInclude));
    
    # for each not included arc
    for (i in 1:length(toInclude))
    {
      # get the position in the links matrix (transpose of the matrix)
      r = ceiling(toInclude[i] / nVar);
      c = toInclude[i] - (r - 1) * nVar;
      
      # check if the child variable has less than 2 parents
      if (sum(links[,c])<2)
      {
        cat('\n \n incluye enlace de la variable ',r,'(padre) a la variable ',c,'(hija)');
        
        
        # add the new link
        links1 = links;
        links1[r,c] = 1;
        
        # check if the tree is a DAG
        DAG1 = igraph :: is_dag(graph_from_adjacency_matrix(links1,mode='directed'));
        
        if (DAG1==TRUE)
        {
          # estimate the parameters of the new model
          newModel = LearnParameters(dataTrain,dataVal,dataMin,dataMax,varType,splitPoints,
                                     maxDegree,links1,fX,fXY,fXYZ,BICpartX,BICpartXY,
                                     BICpartXYZ,verbose);
          
          # upload the distributions
          fX = newModel[[2]];       fXY = newModel[[3]];       fXYZ = newModel[[4]];
          
          # upload the partial BIC
          #BICpartX = newModel[[5]];
          BICpartXY = newModel[[6]];
          BICpartXYZ = newModel[[7]];
          
          # the links and BIC of the new model
          newLinks1[[i]] = links1;
          newBIC1[i] = BICNetwork(links1,BICpartX,BICpartXY,BICpartXYZ);
          #newBIC1[i] = newModel[[8]];
        }
      }
    }
    
    
    
    
    
    # OPTION 2: change the direction of a link
    
    
    cat('\n \n           OPCIÓN 2');
    
    # select the existing links
    toChange = which(t(links)==1);
    
    # initialize a list to save the matrices of links and the BIC
    newLinks2 = vector(mode='list',length=length(toChange));
    newBIC2 = rep(NA,length=length(toChange));
    
    # check if there are existing arcs
    if (length(toChange)>0)
    {
      # for each existing arc
      for (i in 1:length(toChange))
      {
        # get the position in the links matrix (transpose of the matrix)
        r = ceiling(toChange[i] / nVar);
        c = toChange[i] - (r - 1) * nVar;
        
        # check if the child variable has less than 2 parents (change column by row)
        if (sum(links[,r])<2)
        {
          # change the direction of the link
          links1 = links;
          links1[r,c] = 0;     links1[c,r] = 1;
          
          cat('\n \n cambia enlace de ',r,'a',c,' por ',c,'a',r);
          
          # check if the tree is a DAG
          DAG1 = igraph :: is_dag(graph_from_adjacency_matrix(links1,mode='directed'));
          
          if (DAG1==TRUE)
          {
            newModel = LearnParameters(dataTrain,dataVal,dataMin,dataMax,varType,splitPoints,
                                       maxDegree,links1,fX,fXY,fXYZ,BICpartX,BICpartXY,
                                       BICpartXYZ,verbose);
            
            # upload the distributions
            fX = newModel[[2]];       fXY = newModel[[3]];       fXYZ = newModel[[4]];
            
            # upload the partial BIC
            #BICpartX = newModel[[5]];
            BICpartXY = newModel[[6]];
            BICpartXYZ = newModel[[7]];
            
            # the links and BIC of the new model
            newLinks2[[i]] = links1;
            newBIC2[i] = BICNetwork(links1,BICpartX,BICpartXY,BICpartXYZ);
            #newBIC2[i] = newModel[[8]];
          }
        }
      }
    }
    
    
    
    # OPTION 3: delete an existing link
    
    
    cat('\n \n           OPCIÓN 3');
    
    
    # select the links that are previously included
    toDelete = which(t(links)==1);
    
    # initialize a list to save the matrices of links and the BIC
    newLinks3 = vector(mode='list',length=length(toDelete));
    newBIC3 = rep(NA,length=length(toDelete));
    
    # check if there are existing arcs (more than one)
    if (length(toDelete)>1)
    {
      # for each included arc
      for (i in 1:length(toDelete))
      {
        # get the position in the links matrix (transpose of the matrix)
        r = ceiling(toDelete[i] / nVar);
        c = toDelete[i] - (r - 1) * nVar;
        
        cat('\n \n elimina enlace de la variable ',r,'(padre) a la variable ',c,'(hija)');
        
        # delete the existing link
        links1 = links;
        links1[r,c] = 0;
        
        # check if the tree is a DAG (no need to check the number of parents)
        DAG1 = igraph :: is_dag(graph_from_adjacency_matrix(links1,mode='directed'));
        
        if (DAG1==TRUE)
        {
          # estimate the parameters of the new model
          newModel = LearnParameters(dataTrain,dataVal,dataMin,dataMax,varType,splitPoints,
                                     maxDegree,links1,fX,fXY,fXYZ,BICpartX,BICpartXY,
                                     BICpartXYZ,verbose);
          
          # upload the distributions
          fX = newModel[[2]];       fXY = newModel[[3]];       fXYZ = newModel[[4]];
          
          # upload the partial BIC
          #BICpartX = newModel[[5]];
          BICpartXY = newModel[[6]];
          BICpartXYZ = newModel[[7]];
          
          # the links and BIC of the new model
          newLinks3[[i]] = links1;
          newBIC3[i] = BICNetwork(links1,BICpartX,BICpartXY,BICpartXYZ);
          #newBIC3[i] = newModel[[8]];
        }
      }
    }
    
    
    
    # calculate the maximum of the new BIC
    maxBIC = max(c(newBIC1,newBIC2,newBIC3),na.rm=TRUE);
    
    cat('\n Mejor BIC = ',maxBIC);
    
    # check if BIC has increased
    if (maxBIC>BIC)
    {
      # check what kind of change must be made and upload the BIC and the links matrix
      if (maxBIC==max(newBIC1,na.rm=TRUE))
      {
        cat('\n El mejor BIC está en la parte de añadir enlace. Se añade el ',
            toInclude[which.max(newBIC1)]);
        
        BIC = newBIC1[which.max(newBIC1)];
        links = newLinks1[[which.max(newBIC1)]];
      }
      else if (maxBIC==max(newBIC2,na.rm=TRUE))
      {
        cat('\n El mejor BIC está en la parte de cambiar sentido de enlace. Se cambia el ',
            toChange[which.max(newBIC2)]);
        
        BIC = newBIC2[which.max(newBIC2)];
        links = newLinks2[[which.max(newBIC2)]];
      }
      else
      {
        cat('\n El mejor BIC está en la parte de eliminar enlace. Se elimina el ',
            toDelete[which.max(newBIC3)]);
        
        BIC = newBIC3[which.max(newBIC3)];
        links = newLinks3[[which.max(newBIC3)]];
      }
      
      # add one change
      nChanges = nChanges + 1;
    }
    
    else   {   continue = FALSE;   }
    
    
    ShowMatrix(links);
  }
  
  
  # write the final model using the previously estimated distributions
  p = LearnParameters(dataTrain,dataVal,dataMin,dataMax,varType,splitPoints,maxDegree,links,
                             fX,fXY,fXYZ,BICpartX,BICpartXY,BICpartXYZ,verbose);
  model = p[[1]];
  
  
  cat('\n BIC marginales = \n ');
  print(BICpartX)
  cat('\n BIC condicionados a un padre = \n ');
  print(BICpartXY)
  
  
  
  
  # CALCULAR EL BIC DEL DATA TEST
  
  
  
  return(model)
}
############################################################










