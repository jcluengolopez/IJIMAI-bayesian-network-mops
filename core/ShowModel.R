

ShowModelClassification = function(prob)
{
  # it shows the classfication model in the screen in an appropriate format
  # INPUT:  - prob = the TAN model


  cat('\n \n \n          TAN CLASSIFICATION MODEL');

  # type of variables
  varType = prob[[1]];
  cat('\n \n \n Type of variables: 0 = discrete, 1 = continuous, 2 = discretized \n');
  cat(varType);

  # tree structure
  links = prob[[2]]
  cat('\n \n \n Structure of the tree. Parents in rows and children in columns \n');
  print.table(links);

  # different values of the class
  classVal = prob[[length(prob)-1]];

  # sum the columns of links to get the number of parents of each variable
  nParents = apply(links,2,sum);

  # the columns of the discrete and discretized variables are the values of the class
  colDisc = vector(mode='character',length=length(classVal));
  for (j in 1:length(colDisc))   {   colDisc[j] = paste("C =",classVal[j]);   }

  # print variable by variable
  for (i in 1:(length(varType)-1))
  {
    # the model depends on the parents
    if (nParents[i]==0)   # 0 parent (it mus be the node variable)
    {
      cat('\n \n VARIABLE  ',i,', which is the node and it is ');

      # the model depends on the type of variable
      if (varType[i]==1)  # when the variable is continuous
      {
        # the used part of the model
        values = prob[[3*i]];    # values of the variable
        polynom = prob[[3*i+1]];   # conditioned polynomials

        # the data from the polynomials
        cat('Continuous');

        # show the information and the polynomial of each value of the class
        for (j in 1:dim(polynom)[1])
        {
          v = values[,(3*j-2):(3*j)];

          # print the polynomial (use only the coefficients, not NAs)
          cat('\n p(X[',i,'] / C =',classVal[j],') =  ');
          poly = polynom[j,];
          cat('\n poly = ');
          print(as.polynomial(poly[!is.na(poly)]),decreasing=TRUE);
          cat('\n')
          # print the matrix using row and col for the names
          rownames(v) = c("Min   ","Max   ");
          colnames(v) = c("Root","Extreme","Prob. Queue");
          print.table(as.table(v));
        }
      }
      else if (varType[i]==0)   # when the variable is discrete
      {
        # the used part of the model
        dataVal = prob[[3*i]];    # values of the variable
        p = prob[[3*i+1]];    # conditioned probabilities

        # each row is a value of the variable
        row = vector(mode='character',length=length(dataVal));
        for (j in 1:length(row))   {   row[j] = paste("X =",dataVal[j],"   ");   }

        # print the matrix using row and col for the names
        cat('Discrete \n');
        if (length(row) > 1)     {     rownames(p) = row;     }
        if (length(colDisc) > 1)     {     colnames(p) = colDisc;     }

        print.table(as.table(p));
      }
      else   # when the variable is discretized
      {
        # the used part of the model
        splitPoints = prob[[3*i]];    # split points of the variable
        dataVal = 1:(length(splitPoints)-1);    # the rows are 1, 2, 3, 4...
        p = prob[[3*i+1]];    # conditioned probabilities

        # in discretized variables it's possible that this is empty
        row = vector(mode='character',length=length(dataVal));
        for (j in 1:length(row))   {   row[j] = paste("X =",dataVal[j],"   ");   }

        # print the matrix using row and col for the names
        cat('Discretized \n');
        cat('\n Split points in the discretization = ',splitPoints);
        cat('\n');
        if (length(row) > 1 && length(colDisc) > 1)
        {
          rownames(p) = row;               colnames(p) = colDisc;
        }
        print.table(as.table(p));
      }
    }



    else   # 2 parents: a parent apart from the class
    {
      cat('\n \n VARIABLE  ',i,', which is ');

      # obtain the parent of the variable
      posParent = which(links[,i]==1);

      # the model depends on the type of parent
      if (varType[posParent]==1)   # 2 parents (the class and another continuous)
      {
        # the model depends on the type of variable
        if (varType[i]==1)  # when the variable is continuous
        {
          # the used part of the model
          polynom = prob[[3*i+2]];   # conditioned polynomials

          # the data from the polynomials
          cat('Continuous with a continuous parent, which is ',posParent);

          # show the information and the polynomial of each value of both parents
          for (n in 1:dim(polynom)[3])       # values of the class
          {
            # calculate the number of intervals in the discretization in this case
            p = polynom[,1,n];        nInt = length(p[!is.na(p)]);

            # select the used part of the model
            values = prob[[3*i]][(2*n-1):(2*n),1:(3*nInt)];    # values of the variable
            cutPoints = prob[[3*i+1]][n,1:(nInt+1)];   # split points (parent)
            poly = polynom[1:nInt,,n];      # polynomials conditioned to the other parent

            cat('\n Split points of the discretization of the continuous parent = ',
                cutPoints);

            #  for each interval in the discretization of the parent
            for (k in 1:nInt)
            {
              v = values[,(3*k-2):(3*k)];

              # print the polynomial (use only the coefficients, not NAs)
              cat('\n p(X[',i,'] / X[',posParent,'] =',k,' , C =',classVal[n],') =  ');
              #cat('\n k =',k);
              aux = poly[k,];
              cat('\n poly = ');
              print(as.polynomial(aux[!is.na(aux)]),decreasing=TRUE)
              cat('\n');

              # print the matrix using row and col for the names
              rownames(v) = c("Min   ","Max   ");
              colnames(v) = c("Root","Extreme","Prob. Queue");
              print.table(as.table(v));
            }
          }
        }


        else   # when the variable is discrete or discretized
        {
          # the values of the variable depends on its type
          if (varType[i]==0)   # when the variable is discrete
          {
            # the used part of the model
            dataVal = prob[[3*i]];    # values of the variable

            # the data from the polynomials
            cat('Discrete with a continuous parent, which is ',posParent);
          }
          else   # when the variable is discretized
          {
            # the used part of the model
            splitPoints = prob[[3*i]];    # split points of the variable
            dataVal = 1:(length(splitPoints)-1);    # the rows are 1, 2, 3, 4...

            # the data from the polynomials
            cat('Discrized with a continuous parent, which is ',posParent);
            cat('\n Split points in the discretization = ',splitPoints);
          }

          # matrix of probabilities
          dataMatrix = prob[[3*i+2]];

          # show the information and the polynomial of each value of both parents
          for (n in 1:(dim(dataMatrix)[3]))       # values of the class
          {
            # calculate the number of intervals in the discretization in this case
            p = dataMatrix[1,,n];        nInt = length(p[!is.na(p)]);

            # select the used part of the model
            p = dataMatrix[,1:nInt,n];     # probabilities conditioned to the parent
            cutPoints = prob[[3*i+1]][n,1:(nInt+1)];   # split points (parent)

            cat('\n Split points of the discretization of the continuous parent = ',
                cutPoints);

            # each row is a value of the variable
            row = vector(mode='character',length=length(dataVal));
            for (j in 1:length(row))   {   row[j] = paste("X =",dataVal[j],"   ");   }

            # print the matrix using row and col for the names
            cat('\n p(X[',i,'] / X[',posParent,'], C =',classVal[n],') =  \n');
            if (length(row) > 1 && length(cutPoints) > 1)
            {
              rownames(p) = row;               colnames(p) = cutPoints[1:nInt];
            }
            print.table(as.table(p));
          }
        }
      }



      else    # 2 parents (the class and another discrete or discretized)
      {
        # extract the values of the parent and the number of them
        if (varType[posParent]==0)   {   parentVal = prob[[3*posParent]];   }
        else          # discretized
        {
          # the used part of the model
          splitPoints = prob[[3*posParent]];    # split points of the variable
          parentVal = 1:(length(splitPoints)-1);    # the rows are 1, 2, 3, 4...
        }

        nInt = length(parentVal);

        # the model depends on the type of variable
        if (varType[i]==1)  # when the variable is continuous
        {
          # the used part of the model
          polynom = prob[[3*i+1]];   # conditioned polynomials

          # the data from the polynomials
          if (varType[posParent]==0)
          {
            cat('Continuous with a discrete parent, which is ',posParent);
            cat('\n Different values of the discrete parent = ',parentVal);
          }
          else
          {
            cat('Continuous with a discretized parent, which is ',posParent);
            cat('\n Different values of the discretized parent = ',parentVal);
          }

          # show the information and the polynomial of each value of both parents
          for (n in 1:dim(polynom)[3])       # values of the class
          {
            # select the used part of the model
            values = prob[[3*i]][,(3*n-2):(3*n)];    # values of the variable
            poly = polynom[,,n];      # polynomials conditioned to the other parent

            #  for each interval in the discretization of the parent
            for (k in 1:nInt)
            {
              v = values[(2*k-1):(2*k),];

              # print the polynomial (use only the coefficients, not NAs)
              cat('\n p(X[',i,'] / X[',posParent,'] =',k,' , C =',classVal[n],') =  ');
              aux = poly[k,];
              cat('\n poly = ');
              print(as.polynomial(aux[!is.na(aux)]),decreasing=TRUE)
              cat('\n');

              # print the matrix using row and col for the names
              rownames(v) = c("Min   ","Max   ");
              colnames(v) = c("Root","Extreme","Prob. Queue");
              print.table(as.table(v));
            }
          }
        }


        else   # when the variable is discrete or discretized
        {
          # the values of the variable depends on the type
          if (varType[i]==0)   # when the variable is discrete
          {
            # the used part of the model
            dataVal = prob[[3*i]];    # values of the variable

            if (varType[posParent]==0)
            {
              cat('Discrete with a discrete parent, which is ',posParent);
              cat('\n Different values of the discrete parent = ',parentVal);
            }
            else
            {
              cat('Discrete with a discretized parent, which is ',posParent);
              cat('\n Different values of the discretized parent = ',parentVal);
            }
          }
          else   # when the variable is discretized
          {
            # the used part of the model
            splitPoints = prob[[3*i]];    # split points of the variable
            dataVal = 1:(length(splitPoints)-1);    # the rows are 1, 2, 3, 4...

            if (varType[posParent]==0)
            {
              cat('Discretized with a discrete parent, which is ',posParent);
              cat('\n Different values of the discrete parent = ',parentVal);
            }
            else
            {
              cat('Discretized with a discretized parent, which is ',posParent);
              cat('\n Different values of the discretized parent = ',parentVal);
            }
            cat('\n Split points in the discretization = ',splitPoints);
          }

          # matrix of probabilities
          dataMatrix = prob[[3*i+1]];

          # show the information and the polynomial of each value of both parents
          for (n in 1:(dim(dataMatrix)[3]))       # values of the class
          {
            # select the used part of the model
            p = dataMatrix[,,n];     # probabilities conditioned to the parent

            # each row is a value of the variable
            row = vector(mode='character',length=length(dataVal));
            for (j in 1:length(row))   {   row[j] = paste("X =",dataVal[j],"   ");   }

            # print the matrix using row and col for the names
            cat('\n p(X[',i,'] / X[',posParent,'], C =',classVal[n],') =  \n');
            if (length(row) > 1 && length(parentVal) > 1)
            {
              rownames(p) = row;               colnames(p) = parentVal;
            }

            print.table(as.table(p));
          }
        }
      }
    }
  }



  # print the probabilities of the class
  cat('\n \n CLASS \n');
  c = t(prob[[length(prob)]]);
  rownames(c) = "p(C)   ";   colnames(c) = colDisc;
  print.table(as.table(c));
  cat('\n \n');
}
############################################################








ShowModelRegression = function(prob)
{
  # it shows the regression model in the screen in an appropriate format
  # INPUT:  - prob = the TAN model


  cat('\n \n \n          TAN REGRESSION MODEL');

  # type of variables
  varType = prob[[1]];
  cat('\n \n \n Type of variables: 0 = discrete, 1 = continuous, 2 = discretized \n');
  cat(varType);

  # tree structure
  links = prob[[2]]
  cat('\n \n \n Structure of the tree. Parents in rows and children in columns \n');
  print.table(links);

  # sum the columns of links to get the number of parents of each variable
  nParents = apply(links,2,sum);

  # extract the matrix of values and polynomial of the class
  l = length(prob);
  classVal = prob[[l-1]];          classPoly = prob[[l]];

  # minimum and maximum values of the class
  classMin = classVal[1,2];               classMax = classVal[2,2];

  # calculate all the possible cut points in the discretization of the class
  distClass = (classMax - classMin) / 8;
  nInt = 0:8;
  cutPointsClass = classMin + distClass * nInt;

  # print variable by variable
  for (i in 1:(length(varType)-1))
  {
    # the model depends on the parents
    if (nParents[i]==0)   # 0 parent (it mus be the node variable)
    {
      cat('\n \n VARIABLE  ',i,', which is the node and it is ');

      # the model depends on the type of variable
      if (varType[i]==1)  # when the variable is continuous
      {
        # the data from the polynomials
        cat('Continuous');

        # the used part of the model
        values = prob[[3*i]];        # values of the variable
        cutPoints = prob[[3*i+1]]    # split points (class)
        polynom = prob[[3*i+2]];     # conditioned polynomials

        # number of intervals of the discretization
        nInt = length(cutPoints) - 1;

        cat('\n Split points of the discretization of the class = ',cutPoints,'\n');

        #  for each interval in the discretization of the parent
        for (k in 1:nInt)
        {
          v = values[,(3*k-2):(3*k)];

          # print the polynomial (use only the coefficients, not NAs)
          cat('\n p(X[',i,'] / C_',k,') =  ');
          cat('\n poly = ');
          aux = polynom[k,];
          print(as.polynomial(aux[!is.na(aux)]),decreasing=TRUE)
          cat('\n');

          # print the matrix using row and col for the names
          rownames(v) = c("Min   ","Max   ");
          colnames(v) = c("Root","Extreme","Prob. Queue");
          print.table(as.table(v));
        }
      }

      else        # when the variable is discrete or discretized
      {
        # the values of the variable depends on the type
        if (varType[i]==0)   # when the variable is discrete
        {
          # the used part of the model
          dataVal = prob[[3*i]];    # values of the variable

          # the data from the polynomials
          cat('Discrete');
        }
        else   # when the variable is discretized
        {
          # the used part of the model
          splitPoints = prob[[3*i]];    # split points of the variable
          dataVal = 1:(length(splitPoints)-1);    # the rows are 1, 2, 3, 4...

          # the data from the polynomials
          cat('Discrized');
          cat('\n Split points in the discretization of the variable = ',splitPoints);
        }

        # the matrix with the probabilities conditioned to the class
        p = prob[[3*i+2]];
        cutPoints = prob[[3*i+1]]    # split points (class)

        cat('\n Split points of the discretization of the class = ',cutPoints,'\n');

        # each row is a value of the variable
        row = vector(mode='character',length=length(dataVal));
        for (j in 1:length(row))   {   row[j] = paste("X =",dataVal[j],"   ");   }

        # the columns are the discretization of the class
        col = vector(mode='character',length=length(cutPoints)-1);
        for (j in 1:length(col))   {   col[j] = paste("C_",j,"  ");   }

        # print the matrix using row and col for the names
        if (length(row) > 1 && length(col) > 1)
        {
          rownames(p) = row;               colnames(p) = col;
        }
        print.table(as.table(p));
      }
    }




    else   # 2 parents: a parent apart from the class
    {
      cat('\n \n VARIABLE  ',i,', which is ');

      # obtain the parent of the variable
      posParent = which(links[,i]==1);
      parentVal = prob[[3*posParent]];

      # the model depends on the type of parent
      if (varType[posParent]==1)   # 2 parents (the class and another continuous)
      {
        # matrix with the discretization of both parents
        cutPoints = prob[[3*i+1]];

        # the model depends on the type of variable
        if (varType[i]==1)  # when the variable is continuous
        {
          # the used part of the model
          polynom = prob[[3*i+2]];   # conditioned polynomials
          values = prob[[3*i]];    # values of the variable

          # the data from the polynomials
          cat('Continuous with a continuous parent, which is ',posParent);

          # show the matrix with the discretization
          cat('\n Matrix with the discretization of X[',posParent,'] and the class \n');
          print.table(as.table(cutPoints),row.names=FALSE,col.names=FALSE);
          cat('\n (X[',posParent,'] by rows and the class by columns)');

          # show the interval of each value of the variable
          for (k in 1:max(cutPoints))
          {
            # vectors with the cut points of both discretizations
            disc = CutPointsDiscr(cutPoints,parentVal[1,2],parentVal[2,2],classMin,
                                  classMax,k);

            cat('\n ',k,'covers the interval (',disc[1,1],',',disc[1,2],') in X[',
                k,'] and the interval (',disc[2,1],',',disc[2,2],') in the class');

            # matrix with the values of the adjust
            v = values[,(3*k-2):(3*k)];

            # print the polynomial (use only the coefficients, not NAs)
            aux = polynom[k,];
            cat('\n poly = ');
            print(as.polynomial(aux[!is.na(aux)]),decreasing=TRUE)
            cat('\n');

            # print the matrix using row and col for the names
            rownames(v) = c("Min   ","Max   ");
            colnames(v) = c("Root","Extreme","Prob. Queue");
            print.table(as.table(v));
          }
        }

        else   # when the variable is discrete or discretized
        {
          # the values of the variable depends on its type
          if (varType[i]==0)   # when the variable is discrete
          {
            # the used part of the model
            dataVal = prob[[3*i]];    # values of the variable

            # the data from the polynomials
            cat('Discrete with a continuous parent, which is ',posParent);
          }
          else   # when the variable is discretized
          {
            # the used part of the model
            splitPoints = prob[[3*i]];    # split points of the variable
            dataVal = 1:(length(splitPoints)-1);    # the rows are 1, 2, 3, 4...

            # the data from the polynomials
            cat('Discrized with a continuous parent, which is ',posParent);
            cat('\n Split points in the discretization of the variable = ',splitPoints);
          }

          # matrix of probabilities
          dataMatrix = prob[[3*i+2]];

          # show the matrix with the discretization
          cat('\n Matrix with the discretization of X[',posParent,'] and the class \n');
          print.table(as.table(cutPoints),row.names=FALSE,col.names=FALSE);
          cat('\n (X[',posParent,'] by rows and the class by columns)');

          # show the interval of each value of the variable
          for (k in 1:max(cutPoints))
          {
            # vectors with the cut points of both discretizations
            disc = CutPointsDiscr(cutPoints,parentVal[1,2],parentVal[2,2],classMin,
                                  classMax,k);

            cat('\n ',k,'covers the interval (',disc[1,1],',',disc[1,2],') in X[',
                posParent,'] and the interval (',disc[2,1],',',disc[2,2],') in the class');
          }

          # each row is an interval of the discretization
          row = vector(mode='character',length=max(cutPoints));
          for (j in 1:length(row))   {   row[j] = paste("disc_",j,"   ");   }

          # each column is a value of the variable
          col = vector(mode='character',length=length(dataVal));
          for (j in 1:length(col))   {   col[j] = paste("X =",dataVal[j],"  ");   }

          # print the matrix using row and col for the names
          if (length(row) > 1 && length(col) > 1)
          {
            rownames(dataMatrix) = row;               colnames(dataMatrix) = col;
          }
          print.table(as.table(dataMatrix));
        }
      }


      else    # 2 parents (the class and another discrete or discretized)
      {
        # extract the values of the parent
        if (varType[posParent]==0)   {   parentVal = prob[[3*posParent]];   }
        else          # discretized
        {
          # the used part of the model
          splitPoints = prob[[3*posParent]];    # split points of the variable
          parentVal = 1:(length(splitPoints)-1);    # the rows are 1, 2, 3, 4...
        }

        # the model depends on the type of variable
        if (varType[i]==1)  # when the variable is continuous
        {
          # the used part of the model
          polynom = prob[[3*i+2]];   # conditioned polynomials

          if (varType[posParent]==0)
          {
            cat('Continuous with a discrete parent, which is ',posParent);
            cat('\n Different values of the discrete parent = ',parentVal);
          }
          else
          {
            cat('Continuous with a discretized parent, which is ',posParent);
            cat('\n Different values of the discretized parent = ',parentVal);
          }

          # show the information and the polynomial of each value of both parents
          for (n in 1:dim(polynom)[3])       # values of the discrete parent
          {
            # calculate the number of intervals in the discretization in this case
            p = polynom[,1,n];        nInt = length(p[!is.na(p)]);

            # select the used part of the model
            values = prob[[3*i]][(2*n-1):(2*n),1:(3*nInt)];    # values of the variable
            cutPoints = prob[[3*i+1]][n,1:(nInt+1)];   # split points (parent)
            poly = polynom[1:nInt,,n];      # polynomials conditioned to the class

            cat('\n Split points of the discretization of the class = ',cutPoints,'\n');

            #  for each interval in the discretization of the class
            for (k in 1:nInt)
            {
              v = values[,(3*k-2):(3*k)];

              # print the polynomial (use only the coefficients, not NAs)
              cat('\n p(X[',i,'] / X[',posParent,'] =',parentVal[n],' , C_',k,') =  ');
              aux = poly[k,];
              cat('\n poly = ');
              print(as.polynomial(aux[!is.na(aux)]),decreasing=TRUE)
              cat('\n');

              # print the matrix using row and col for the names
              rownames(v) = c("Min   ","Max   ");
              colnames(v) = c("Root","Extreme","Prob. Queue");
              print.table(as.table(v));
            }
          }
        }


        else     # when the variable is discrete or discretized
        {
          # the values of the variable depends on its type
          if (varType[i]==0)   # when the variable is discrete
          {
            # the used part of the model
            dataVal = prob[[3*i]];    # values of the variable

            if (varType[posParent]==0)
            {
              cat('Discrete with a discrete parent, which is ',posParent);
              cat('\n Different values of the discrete parent = ',parentVal);
            }
            else
            {
              cat('Discrete with a discretized parent, which is ',posParent);
              cat('\n Different values of the discretized parent = ',parentVal);
            }
          }
          else   # when the variable is discretized
          {
            # the used part of the model
            splitPoints = prob[[3*i]];    # split points of the variable
            dataVal = 1:(length(splitPoints)-1);    # the rows are 1, 2, 3, 4...

            if (varType[posParent]==0)
            {
              cat('Discretized with a discrete parent, which is ',posParent);
              cat('\n Different values of the discrete parent = ',parentVal);
            }
            else
            {
              cat('Discretized with a discretized parent, which is ',posParent);
              cat('\n Different values of the discretized parent = ',parentVal);
            }
            cat('\n Split points in the discretization of teh variable = ',splitPoints);
          }

          # matrix of probabilities
          dataMatrix = prob[[3*i+2]];

          # show the information and the polynomial of each value of both parents
          for (n in 1:(dim(dataMatrix)[3]))       # values of the discrete parent
          {
            # calculate the number of intervals in the discretization in this case
            p = dataMatrix[1,,n];        nInt = length(p[!is.na(p)]);

            # select the used part of the model
            p = dataMatrix[,1:nInt,n];     # probabilities conditioned to the class
            cutPoints = prob[[3*i+1]][n,1:(nInt+1)];   # split points (class)

            cat('\n Split points of the discretization of the class = ',cutPoints,'\n');

            # each row is a value of the variable
            row = vector(mode='character',length=length(dataVal));
            for (j in 1:length(row))   {   row[j] = paste("X =",dataVal[j],"   ");   }

            # each row is a value of the variable
            col = vector(mode='character',length=nInt);
            for (j in 1:length(col))   {   col[j] = paste("C_",j,"  ");   }

            # print the matrix using row and col for the names
            cat('\n p(X[',i,'] / X[',posParent,'] =',parentVal[n],', C) =  \n');
            if (length(row) > 1 && length(col) > 1)
            {
              rownames(p) = row;               colnames(p) = col;
            }
            print.table(as.table(p));
          }
        }
      }
    }
  }


  # print the probabilities of the class
  cat('\n \n CLASS \n');
  cat('\n p(C) =  ');

  cat('\n poly = ');
  print(as.polynomial(classPoly[!is.na(classPoly)]),decreasing=TRUE)
  cat('\n');

  # print the matrix using row and col for the names
  rownames(classVal) = c("Min   ","Max   ");
  colnames(classVal) = c("Root","Extreme","Prob. Queue");
  print.table(as.table(classVal));

  cat('\n \n');

}
############################################################









ShowPoly = function(x,y,poly)
{
  # it plots the polynomial using the x-values
  # INPUT:  - (x,y) = coordinates of the points to adjust
  #         - poly = the coefficients of the polynomial


  # calculate the y-values
  yPoly = predict(as.polynomial(poly),x);

  # do the graph
  plot(x,yPoly,type='l',col="red",lwd=2,xaxt='n');
  points(x,y,col='blue');
  abline(h=0,col='green',lwd=2);

}
############################################################







ShowMatrix = function(x)
{
  # it prints a matrix
  # INPUT:  - x = a matrix


  dim = dim(x);

  cat('\n Matriz: \n');
  for (i in 1:nrow(x))
  {
    for (j in 1:ncol(x))   {   cat('  ',x[i,j]);   }
    cat('\n');
  }
}




PlotPolynomialFit = function(polynomial,values,yLim)
{
  # it plots a polynomial fit with the tails, in case it has got
  # IMPUT: - values = matrix with two rows and 3*C.val columns. It contains the
  #        min and max root, the min and max value and the probab. of the tails
  #        - polynomial = matrix with the coefficients of the polynomial  for
  #        each value of the class
  #        - yLim = maximum value of the y-axis


  # separate roots, extremes and probabilities of the queues
  roots = values[,1];        extremes = values[,2];        probTails = values[,3];

  # create a sequence for the central part
  x = seq(roots[1],roots[2],length.out=1000);

  # calculate the y-values
  y = predict(as.polynomial(polynomial),x);

  # calculate the maximum value of the y-axis in case it is missing
  if (missing(yLim))   {   yLim = max(y);   }

  # create a sequence for the left tail
  if (extremes[1] < roots[1])
  {
    xLeftTail = seq(extremes[1],roots[1],length.out=100);

    # add this sequence to the x
    x = append(xLeftTail,x);

    # add the y-part
    y = append(rep(probTails[1],100),y);
  }

  # create a sequence for the right tail
  if (extremes[2] > roots[2])
  {
    xRightTail = seq(roots[2],extremes[2],length.out=100);

    # add this sequence to the x
    x = append(x,xRightTail);

    # add the y-part
    y = append(y,rep(probTails[2],100));
  }

  # do the graph
  plot(x,y,type='l',col="#0072B2",lwd=3,xlab='',ylab='',ylim=c(0,yLim));
  #points(x,y,type='l',col="#a23d75",lwd=3)
  #abline(h=0,col='green',lwd=2);
}
############################################################









