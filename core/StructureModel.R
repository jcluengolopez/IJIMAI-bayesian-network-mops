

StructureModel = function(prob,namesVar,varType,discretize,discretizeType,method,verbose=FALSE)
{
  # it changes the structure of a model to access easily to each variable, probability
  # and the tree
  # INPUT  - prob = a list with the structure of the model, values and probabilities 
  #        - namesVar = vector with the names of each variable in the database
  #        - varType = a vector with 0 (discrete), 1 (continuous), 2 (discretized)
  #        - discretize = 0 = never, 1 = 20 - 5% of data, 2 = more than 20 values
  #        - discretizeType (optional) = it decides the type of discretization. 0=dynamic,
  #       1=equal width
  #        - method (optional) = a string with the method ('TAN' or 'NaiveBayes')
  #        - verbose = show the 'cats' on screen
  # OUTPUT  - TAN = new list with the structure of the model, values and probabilities 
  
  
  # number of variables
  nVar = length(varType);
  
  # create the list (one for the prob of each variable, names and tree)
  TAN = vector("list",nVar+3);
  
  # the kind of adjusted model
  TAN[[1]]  = TypeModel(discretize,discretizeType,varType[nVar],method);
  
  # create the new vector with the type of variables
  varType2 = vector('character',nVar);
  varType2[which(varType==0)] = 'discrete';
  varType2[which(varType==1)] = 'continuous';
  varType2[which(varType==2)] = 'pseudo-discrete';
  
  # create the first list for the names and type of variables
  TAN[[2]] = list(names=namesVar,varType=varType2);
  
  # convert the adjacency matrix into a matrix with the arcs (like bnlearn)
  TAN[[3]] = ChangeTreeMatrix(prob[[2]],namesVar);
  
  # for each variable except the class
  for (i in 1:(nVar-1))
  {
    # parents and children of each variable
    parents = which(prob[[2]][,i]==1);
    parentsName = namesVar[parents];
    children = namesVar[which(prob[[2]][i,]==1)];
    
    # information about the adjust and type of parents
    adjust = AdjustType(i,prob,parents,varType);
    
    # include in the model a list with all the information about that variable
    TAN[[i+3]] = list(parents=parentsName,children=children,adjust=adjust[[1]],
                      values=adjust[[3]],prob=adjust[[2]],parentValues=adjust[[4]]);
  }
  
  # the list for the class
  if (varType[nVar]==0)   {   typeClass = 'discrete';   }
  else   {   typeClass = 'continuous';   }
  TAN[[nVar+3]] = list(children=namesVar[1:(nVar-1)],adjust=typeClass,
                       values=prob[[3*nVar]],prob=prob[[3*nVar+1]]);
  
  # name all the positions of the main list
  names(TAN) = c('model','nodes','arcs',namesVar);
  
  return(TAN);
}
############################################################








StructureModelIverse = function(prob,verbose=FALSE)
{
  # it changes the structure of a model returning to the original
  # INPUT: - prob = the model with attributes prepared to the used
  #        - verbose = show the 'cats' on screen
  # OUTPUT: - TAN = list with the structure, probabilities and values of the model 
  
  
  # type and number of variables
  varType = prob$nodes$varType;
  nVar = length(varType);
  
  # create the new vector with the type of variables
  varType2 = vector('numeric',nVar);
  varType2[which(varType=='discrete')] = 0;
  varType2[which(varType=='continuous')] = 1;
  varType2[which(varType=='pseudo-discrete')] = 2;
  
  # name of variables
  namesVar = prob$nodes$names;
  
  # create the list (one for the prob of each variable, names and tree)
  TAN = vector("list",nVar+2);
  
  # the type of variable in the first position
  TAN[[1]] = varType2;
  
  # the adyacency matrix in the second position
  TAN[[2]] = ChangeTreeMatrixInverse(prob$arcs,namesVar);
  
  # for each variable except the class
  for (i in 1:(length(varType)-1))
  {
    TAN[[3*i]] = prob[[i+3]]$values;                # values
    TAN[[3*i+1]] = prob[[i+3]]$parentValues;        # cut points parent
    TAN[[3*i+2]] = prob[[i+3]]$prob;                # probabilities
    
    # when the parents are discrete, the probabilities are saved in the second position
    # and the third position is NULL
    if (length(TAN[[3*i+1]])==0)
    {
      TAN[[3*i+1]] = TAN[[3*i+2]];
      TAN[[3*i+2]] = c();
    }
  }
  
  # the class
  TAN[[3*nVar]] = prob[[nVar+3]]$values;        # values
  TAN[[3*nVar+1]] = prob[[nVar+3]]$prob;        # probabilities
  
  return(TAN);
}
############################################################







  


TypeModel = function(discretize,discretizeType,class,method)
{
  # it creates a string explaining the type of model created
  # INPUT  - discretize = 0 = never, 1 = 20 - 5% of data, 2 = more than 20 values
  #        - discretizeType (optional) = it decides the type of discretization. 0=dynamic,
  #       1=equal width
  #        - class = 0 for classification and 1 for regression
  #        - method (optional) = a string with the method ('TAN' or 'NaiveBayes')
  # OUTPUT  - model = string with the type of model
  
  
  # define the method
  if (method=='TAN')   {   method = 'TAN';   }
  else    {   method= 'Naive Bayes';   }
  
  # define the type of model
  if (class==0)   {   model = 'classification model';   }
  else    {   model = 'regression model';   }
  
  # define when it discretizes
  if (discretize==0)   {   disc = 'with no discretization';   discType = ''   }
  else
  {
    if (discretize==1)   {   disc = 'with partial';   }
    else   {   disc = 'with total';   }
    
    # define the type of discretization
    if (discretizeType==0)   {   discType = 'dynamic discretization';   }
    else   {   discType = 'equal-width discretization';   }
  }
  
  # create the string
  sol = paste(method,model,disc,discType);
  
  return(sol);
}
############################################################






ChangeTreeMatrix = function(adjacency,names)
{
  # it converts the adjacency matrix of the tree into a matrix with two columns
  # describing the arcs like the matrix used in bnlearn
  # INPUT  - adjacency = adjacency matrix
  #        - names = vector with the name of each variable
  # OUTPUT  - arcs = matrix with two columns describing the arcs
  
  
  # create a vector for the beginning and another for the ending
  fromArcs = vector('character');          toArcs = fromArcs;
  
  # for each column of the adjacency matrix
  for (j in 1:dim(adjacency)[2])
  {
    # extract the position of the 1 in the column
    pos = which(adjacency[,j]==1);
    
    # add the arc that joins both variables
    if (length(pos)>0)
    {
      fromArcs = c(fromArcs,names[pos]);      toArcs = c(toArcs,rep(names[j],length(pos)));
    }
  }
  
  # join both vectors as columns of the matrix
  arcs = cbind(fromArcs,toArcs);
  colnames(arcs) = c('from','to')
  
  return(arcs);
}
############################################################







ChangeTreeMatrixInverse = function(arcs,names)
{
  # it converts the two-columns matrix with the arcs into an adjacency matrix to use
  # during the process
  # INPUT: - arcs = matrix with two columns describing the arcs
  #        - names = vector with the name of each variable
  # OUTPUT: - adjacency = adjacency matrix
  
  
  # initialize the adyacency matrix with 0 in each position
  adyacency = matrix(0,nrow=length(names),ncol=length(names));
  
  # for each variable
  for (i in 1:length(names))
  {
    # check if it is the parent of any variable
    parents = which(names[i]==arcs[,1]);
    
    if(length(parents)>0)
    {
      # for each son of that variable
      for (j in 1:length(parents))
      {
        # each parent and children
        parent = which(names==arcs[parents[j],1]);
        children = which(names==arcs[parents[j],2]);
        
        # put a 1 in the position of the (parent,children)
        adyacency[parent,children] = 1;
      }
    }
  }
  
  return(adyacency);
}
############################################################










AdjustType = function(pos,prob,parents,varType)
{
  # it gets the type of adjust made according to the kind of variable it and its parents
  # are (continuous or discrete) and the number of parents
  # INPUT: - pos = position of the studied variable
  #        - prob = a list with the structure of the model, values and probabilities 
  #        - parents = vector with the positions of the parents
  #        - varType = a vector with 0 (discrete), 1 (continuous), 2 (discretized)
  # OUTPUT: - type = string with the type of adjust
  #         - p = a list with the structure of the model, values and probabilities
  #         - values = matrix with two rows and 3*C.val columns. It contains the
  #         min and max root, the min and max value and the probab. of the queues
  #         - cutPoints = a vector with the cut points in the discretization of C
  
  
  if (length(parents)==1)     # 1 parent
  {
    if (varType[parents]==0 && varType[pos]==0)     # discr var with 1 discr parent
    {
      type = 'discrete variable with one discrete parent';
      p = prob[[3*pos+1]];     values = prob[[3*pos]];     cutPoints = c();
    }
    else if (varType[parents]==0 && varType[pos]==1)     # cont var with 1 discr parent
    {
      type = 'continuous variable with one discrete parent';
      p = prob[[3*pos+1]];     values = prob[[3*pos]];     cutPoints = c();
    }
    else if (varType[parents]==0 && varType[pos]==2)     # disctz var with 1 discr parent
    {
      type = 'pseudo-discrete variable with one discrete parent';
      p = prob[[3*pos+1]];     values = prob[[3*pos]];     cutPoints = c();
    }
    else if (varType[parents]==1 && varType[pos]==0)     # discr var with 1 cont parent
    {
      type = 'discrete variable with one continuous parent';
      p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
    }
    else if (varType[parents]==1 && varType[pos]==1)     # cont var with 1 cont parent
    {
      type = 'continuous variable with one continuous parent';
      p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
    }
    else if (varType[parents]==1 && varType[pos]==2)     # disctz var with 1 cont parent
    {
      type = 'pseudo-discrete variable with one continuous parent';
      p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
    }
    else if (varType[parents]==2 && varType[pos]==0)     # discr var with 1 disctz parent
    {
      type = 'discrete variable with one pseudo-discrete parent';
      p = prob[[3*pos+1]];     values = prob[[3*pos]];     cutPoints = c();
    }
    else if (varType[parents]==2 && varType[pos]==1)     # cont var with 1 disctz parent
    {
      type = 'continuous variable with one pseudo-discrete parent';
      p = prob[[3*pos+1]];     values = prob[[3*pos]];     cutPoints = c();
    }
    else if (varType[parents]==2 && varType[pos]==2)     # disctz var with 1 disctz parent
    {
      type = 'pseudo-discrete variable with one pseudo-discrete parent';
      p = prob[[3*pos+1]];     values = prob[[3*pos]];     cutPoints = c();
    }
  }
  
  else     # 2 parents
  {
    if (varType[parents[1]]==0)   # one discr parent
    {
      if (varType[parents[2]]==0 && varType[pos]==0)    # discr var with 2 disc parents
      {
        type = 'discrete variable with two discrete parents';
        p = prob[[3*pos+1]];     values = prob[[3*pos]];     cutPoints = c();
      }
      else if (varType[parents[2]]==0 && varType[pos]==1)   # cont var with 2 disc parents
      {
        type = 'continuous variable with two discrete parents';
        p = prob[[3*pos+1]];     values = prob[[3*pos]];     cutPoints = c();
      }
      else if (varType[parents[2]]==0 && varType[pos]==2)    # disctz var with 2 disc parents
      {
        type = 'pseudo-discrete variable with two discrete parents';
        p = prob[[3*pos+1]];     values = prob[[3*pos]];     cutPoints = c();
      }
      else if (varType[parents[2]]==1 && varType[pos]==0)   # disc var with 1 disc & 1 cont
      {
        type = 'discrete variable with one discrete and one continuous parent';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
      }
      else if (varType[parents[2]]==1 && varType[pos]==1)   # cont var with 1 disc & 1 cont
      {
        type = 'continuous variable with one discrete and one continuous parent';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
      }
      else if (varType[parents[2]]==1 && varType[pos]==2)   # distz var with 1 disc & 1 cont
      {
        type = 'pseudo-discrete variable with one discrete and one continuous parent';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
      }
      else if (varType[parents[2]]==2 && varType[pos]==0)   # disc var with 1 disc & 1 dstz
      {
        type = 'discrete variable with one discrete and one pseudo-discrete parent';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = c();
      }
      else if (varType[parents[2]]==2 && varType[pos]==1)   # cont var with 1 disc & 1 dstz
      {
        type = 'continuous variable with one discrete and one pseudo-discrete parent';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = c();
      }
      else if (varType[parents[2]]==2 && varType[pos]==2)   # distz var with 1 disc & 1 dstz
      {
        type = 'pseudo-discrete variable with one discrete and one pseudo-discrete parent';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = c();
      }
    }
    
    if (varType[parents[1]]==1)   # one cont parent
    {
      if (varType[parents[2]]==0 && varType[pos]==0)    # discr var with 1 disc & 1 cont
      {
        type = 'discrete variable with one discrete and one continuous parent';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
      }
      else if (varType[parents[2]]==0 && varType[pos]==1)   # cont var with 1 disc & 1 cont
      {
        type = 'continuous variable with one discrete and one continuous parent';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
      }
      else if (varType[parents[2]]==0 && varType[pos]==2)    # distz var with 1 disc & 1 cont
      {
        type = 'pseudo-discrete variable with one discrete and one continuous parent';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
      }
      else if (varType[parents[2]]==1 && varType[pos]==0)   # disc var with with 2 cont
      {
        type = 'discrete variable with two continuous parents';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
      }
      else if (varType[parents[2]]==1 && varType[pos]==1)   # cont var with 2 cont
      {
        type = 'continuous variable with two continuous parents';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
      }
      else if (varType[parents[2]]==1 && varType[pos]==2)   # distz var with 2 cont
      {
        type = 'pseudo-discrete variable with two continuous parents';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
      }
      else if (varType[parents[2]]==2 && varType[pos]==0)   # disc var with 1 cont & 1 dstz
      {
        type = 'discrete variable with one continuous and one pseudo-discrete parent';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
      }
      else if (varType[parents[2]]==2 && varType[pos]==1)   # cont var with 1 cont & 1 dstz
      {
        type = 'continuous variable with one continuous and one pseudo-discrete parent';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
      }
      else if (varType[parents[2]]==2 && varType[pos]==2)   # distz var with 1 cont & 1 dstz
      {
        type = 'pseudo-discrete variable with one continuous and one pseudo-discrete parent';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
      }
    }
    
    if (varType[parents[1]]==2)   # one disctz parent
    {
      if (varType[parents[2]]==0 && varType[pos]==0)    # discr var with 1 disc & 1 dstz
      {
        type = 'discrete variable with one discrete and one pseudo-discrete paren';
        p = prob[[3*pos+1]];     values = prob[[3*pos]];     cutPoints = c();
      }
      else if (varType[parents[2]]==0 && varType[pos]==1)   # cont var with 1 disc & 1 dstz
      {
        type = 'continuous variable with one discrete and one pseudo-discrete paren';
        p = prob[[3*pos+1]];     values = prob[[3*pos]];     cutPoints = c();
      }
      else if (varType[parents[2]]==0 && varType[pos]==2)    # distz var with 1 disc & 1 dstz
      {
        type = 'pseudo-discrete variable with one discrete and one pseudo-discrete paren';
        p = prob[[3*pos+1]];     values = prob[[3*pos]];     cutPoints = c();
      }
      else if (varType[parents[2]]==1 && varType[pos]==0)   # disc var with 1 distz & 1 cont
      {
        type = 'discrete variable with one pseudo-discrete and one continuous parent';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
      }
      else if (varType[parents[2]]==1 && varType[pos]==1)   # cont var with 1 distz & 1 cont
      {
        type = 'continuous variable with one pseudo-discrete and one continuous parent';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
      }
      else if (varType[parents[2]]==1 && varType[pos]==2)   # distz var with 1 distz & 1 cont
      {
        type = 'pseudo-discrete variable with one pseudo-discrete and one continuous parent';
        p = prob[[3*pos+2]];     values = prob[[3*pos]];     cutPoints = prob[[3*pos+1]];
      }
      else if (varType[parents[2]]==2 && varType[pos]==0)   # disc var with 2 dsctz
      {
        type = 'discrete variable with two pseudo-discrete parents';
        p = prob[[3*pos+1]];     values = prob[[3*pos]];     cutPoints = c();
      }
      else if (varType[parents[2]]==2 && varType[pos]==1)   # cont var with 2 dsctz
      {
        type = 'continuous variable with two pseudo-discrete parents';
        p = prob[[3*pos+1]];     values = prob[[3*pos]];     cutPoints = c();
      }
      else if (varType[parents[2]]==2 && varType[pos]==2)   # distz var with 2 dsctz
      {
        type = 'pseudo-discrete variable with two pseudo-discrete parents';
        p = prob[[3*pos+1]];     values = prob[[3*pos]];     cutPoints = c();
      }
    }
  }
  
  return(list(type,p,values,cutPoints));
}
############################################################









