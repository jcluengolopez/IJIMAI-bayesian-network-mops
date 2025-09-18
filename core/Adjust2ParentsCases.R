

Cases2x2 = function(F22,A22,poly,values,startP,startC,d,dataMin,dataMax,n,k,maxDegree,
                    verbose=FALSE,plots=FALSE)
{
  # it checks if and adjust that has to be done was previously done and stored in
  # the matrices F22 and A22. In that case, matrices are updated
  # INPUT: - F22 = polynomial adjust of a discretization 2x2 (blocks of 4x4)
  #        - A22 = values of the adjust of a discretization 2x2
  #        - startP, startC = firs index of the block in each parent
  #        - d = selected data to adjust
  #        - dataMin, dataMax = the minimum and maximum value in the variable
  #        - n = position of the variable in the dataframe
  #        - k = position of the adjust in the matrix
  #        - maxDegree = maximum degree for the polynomial fit
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - poly = vector with the coefficients of the polynomial for d
  #         - values = matrix with 2 rows and 3 columns. It contains the
  #         min and max root, the min and max value and the probab. of the queues
  #         - F22 = array 2x2x11 for the coefficients of the polynomials
  #         - A22 = array 4x6 for the values of the adjust
  
  
  # check the 4 block where the adjust belongs
  if (startP==1 && startC==1)   # first of a 4x4 block
  {
    # the adjust was previously done
    poly = F22[1,1,];                values = A22[1:2,1:3];
  }
  else if (startP==1 && startC==5)   # second of a 4x4 block
  {
    # the adjust was previously done
    poly = F22[1,2,];                values = A22[1:2,4:6];
  }
  else if (startP==5 && startC==1)   # third of a 4x4 block
  {
    # the adjust was previously done
    poly = F22[2,1,];                values = A22[3:4,1:3];
  }
  else if (startP==5 && startC==5)   # fourth of a 4x4 block
  {
    # the adjust was previously done
    poly = F22[2,2,];                values = A22[3:4,4:6];
  }
  
  return(list(poly,values,F22,A22));
}
############################################################







Cases2x4 = function(F24,A24,poly,values,startP,startC,d,dataMin,dataMax,n,k,maxDegree,
                    verbose=FALSE,plots=FALSE)
{
  # it checks if and adjust that has to be done was previously done and stored in
  # the matrices F24 and A24. In that case, matrices are updated
  # INPUT: - F24 = polynomial adjust of a discretization 2x4 (blocks of 4x2)
  #        - A24 = values of the adjust of a discretization 2x4
  #        - startP, startC = firs index of the block in each parent
  #        - d = selected data to adjust
  #        - dataMin, dataMax = the minimum and maximum value in the variable
  #        - n = position of the variable in the dataframe
  #        - k = position of the adjust in the matrix
  #        - maxDegree = maximum degree for the polynomial adjust
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - poly = vector with the coefficients of the polynomial for d
  #         - values = matrix with 2 rows and 3 columns. It contains the
  #         min and max root, the min and max value and the probab. of the queues
  #         - F24 = array 2x4x11 for the coefficients of the polynomials
  #         - A24 = array 4x12 for the values of the adjust
  
  
  # create the vector for the adjust (simple polynomial adjust)
  r = rep(0,length(d));
  
  # check the 8 block where the adjust belongs
  if (startP==1 && startC==1)   # first of a 4x2 block
  {
    # when the adjust was previously done
    if (is.na(F24[1,1,1])==FALSE)
    {
      poly = F24[1,1,];              values = A24[1:2,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F24[1,1,] = prob[[2]];         A24[1:2,1:3] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==3)   # second of a 4x2 block
  {
    # when the adjust was previously done
    if (is.na(F24[1,2,1])==FALSE)
    {
      poly = F24[1,2,];              values = A24[1:2,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F24[1,2,] = prob[[2]];         A24[1:2,4:6] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==5)   # third of a 4x2 block
  {
    # when the adjust was previously done
    if (is.na(F24[1,3,1])==FALSE)
    {
      poly = F24[1,3,];              values = A24[1:2,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F24[1,3,] = prob[[2]];         A24[1:2,7:9] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==7)   # fourth of a 4x2 block
  {
    # when the adjust was previously done
    if (is.na(F24[1,4,1])==FALSE)
    {
      poly = F24[1,4,];              values = A24[1:2,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F24[1,4,] = prob[[2]];         A24[1:2,10:12] = prob[[1]];
    }
  }
  
  
  else if (startP==5 && startC==1)   # fifth of a 4x2 block
  {
    # when the adjust was previously done
    if (is.na(F24[2,1,1])==FALSE)
    {
      poly = F24[2,1,];              values = A24[3:4,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F24[2,1,] = prob[[2]];         A24[3:4,1:3] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==3)   # sixth of a 4x2 block
  {
    # when the adjust was previously done
    if (is.na(F24[2,2,1])==FALSE)
    {
      poly = F24[2,2,];              values = A24[3:4,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F24[2,2,] = prob[[2]];         A24[3:4,4:6] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==5)   # seventh of a 4x2 block
  {
    # when the adjust was previously done
    if (is.na(F24[2,3,1])==FALSE)
    {
      poly = F24[2,3,];              values = A24[3:4,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F24[2,3,] = prob[[2]];         A24[3:4,7:9] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==7)   # eighth of a 4x2 block
  {
    # when the adjust was previously done
    if (is.na(F24[2,4,1])==FALSE)
    {
      poly = F24[2,4,];              values = A24[3:4,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F24[2,4,] = prob[[2]];         A24[3:4,10:12] = prob[[1]];
    }
  }
  
  return(list(poly,values,F24,A24));
}
############################################################








Cases4x2 = function(F42,A42,poly,values,startP,startC,d,dataMin,dataMax,n,k,maxDegree,
                    verbose=FALSE,plots=FALSE)
{
  # it checks if and adjust that has to be done was previously done and stored in
  # the matrices F42 and A42. In that case, matrices are updated
  # INPUT: - F42 = polynomial adjust of a discretization 4x2 (blocks of 2x4)
  #        - A42 = values of the adjust of a discretization 4x2
  #        - startP, startC = firs index of the block in each parent
  #        - d = selected data to adjust
  #        - dataMin, dataMax = the minimum and maximum value in the variable
  #        - n = position of the variable in the dataframe
  #        - k = position of the adjust in the matrix
  #        - maxDegree = maximum degree for the polynomial fit
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - poly = vector with the coefficients of the polynomial for d
  #         - values = matrix with 2 rows and 3 columns. It contains the
  #         min and max root, the min and max value and the probab. of the queues
  #         - F42 = array 4x2x11 for the coefficients of the polynomials
  #         - A42 = array 8x6 for the values of the adjust
  
  
  # create the vector for the adjust (simple polynomial adjust)
  r = rep(0,length(d));
  
  # check the 8 block where the adjust belongs
  if (startP==1 && startC==1)   # first of a 2x4 block
  {
    # when the adjust was previously done
    if (is.na(F42[1,1,1])==FALSE)
    {
      poly = F42[1,1,];              values = A42[1:2,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F42[1,1,] = prob[[2]];         A42[1:2,1:3] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==5)   # second of a 2x4 block
  {
    # when the adjust was previously done
    if (is.na(F42[1,2,1])==FALSE)
    {
      poly = F42[1,2,];              values = A42[1:2,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F42[1,2,] = prob[[2]];         A42[1:2,4:6] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==1)   # third of a 2x4 block
  {
    # when the adjust was previously done
    if (is.na(F42[2,1,1])==FALSE)
    {
      poly = F42[2,1,];              values = A42[3:4,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F42[2,1,] = prob[[2]];         A42[3:4,1:3] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==5)   # fourth of a 2x4 block
  {
    # when the adjust was previously done
    if (is.na(F42[2,2,1])==FALSE)
    {
      poly = F42[2,2,];              values = A42[3:4,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F42[2,2,] = prob[[2]];         A42[3:4,4:6] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==1)   # fifth of a 2x4 block
  {
    # when the adjust was previously done
    if (is.na(F42[3,1,1])==FALSE)
    {
      poly = F42[3,1,];              values = A42[5:6,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F42[3,1,] = prob[[2]];         A42[5:6,1:3] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==5)   # sixth of a 2x4 block
  {
    # when the adjust was previously done
    if (is.na(F42[3,2,1])==FALSE)
    {
      poly = F42[3,2,];              values = A42[5:6,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F42[3,2,] = prob[[2]];         A42[5:6,4:6] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==1)   # seventh of a 2x4 block
  {
    # when the adjust was previously done
    if (is.na(F42[4,1,1])==FALSE)
    {
      poly = F42[4,1,];              values = A42[7:8,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F42[4,1,] = prob[[2]];         A42[7:8,1:3] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==5)   # eighth of a 2x4 block
  {
    # when the adjust was previously done
    if (is.na(F42[4,2,1])==FALSE)
    {
      poly = F42[4,2,];              values = A42[7:8,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F42[4,2,] = prob[[2]];         A42[7:8,4:6] = prob[[1]];
    }
  }
  
  return(list(poly,values,F42,A42));
}
############################################################









Cases4x4 = function(F44,A44,poly,values,startP,startC,d,dataMin,dataMax,n,k,maxDegree,
                    verbose=FALSE,plots=FALSE)
{
  # it checks if and adjust that has to be done was previously done and stored in
  # the matrices F44 and A44. In that case, matrices are updated
  # INPUT: - F44 = polynomial adjust of a discretization 4x4 (blocks of 2x2)
  #        - A44 = values of the adjust of a discretization 4x4
  #        - startP, startC = firs index of the block in each parent
  #        - d = selected data to adjust
  #        - dataMin, dataMax = the minimum and maximum value in the variable
  #        - n = position of the variable in the dataframe
  #        - k = position of the adjust in the matrix
  #        - maxDegree = maximum degree for the polynomial fit
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - poly = vector with the coefficients of the polynomial for d
  #         - values = matrix with 2 rows and 3 columns. It contains the
  #         min and max root, the min and max value and the probab. of the queues
  #         - F44 = array 4x4x11 for the coefficients of the polynomials
  #         - A44 = array 8x12 for the values of the adjust
  
  
  # create the vector for the adjust (simple polynomial adjust)
  r = rep(0,length(d));
  
  # check the 16 block where the adjust belongs
  if (startP==1 && startC==1)   # first of a 2x2 block
  {
    # when the adjust was previously done
    if (is.na(F44[1,1,1])==FALSE)
    {
      poly = F44[1,1,];              values = A44[1:2,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F44[1,1,] = prob[[2]];         A44[1:2,1:3] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==3)   # second of a 2x2 block
  {
    # when the adjust was previously done
    if (is.na(F44[1,2,1])==FALSE)
    {
      poly = F44[1,2,];              values = A44[1:2,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F44[1,2,] = prob[[2]];         A44[1:2,4:6] = prob[[1]];
    }
  }
  
  if (startP==1 && startC==5)   # third of a 2x2 block
  {
    # when the adjust was previously done
    if (is.na(F44[1,3,1])==FALSE)
    {
      poly = F44[1,3,];              values = A44[1:2,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F44[1,3,] = prob[[2]];         A44[1:2,7:9] = prob[[1]];
    }
  }
  
  if (startP==1 && startC==7)   # fourth of a 2x2 block
  {
    # when the adjust was previously done
    if (is.na(F44[1,4,1])==FALSE)
    {
      poly = F44[1,4,];              values = A44[1:2,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F44[1,4,] = prob[[2]];         A44[1:2,10:12] = prob[[1]];
    }
  }
  
  
  else if (startP==3 && startC==1)   # fifth of a 2x2 block
  {
    # when the adjust was previously done
    if (is.na(F44[2,1,1])==FALSE)
    {
      poly = F44[2,1,];              values = A44[3:4,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F44[2,1,] = prob[[2]];         A44[3:4,1:3] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==3)   # sixth of a 2x2 block
  {
    # when the adjust was previously done
    if (is.na(F44[2,2,1])==FALSE)
    {
      poly = F44[2,2,];              values = A44[3:4,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F44[2,2,] = prob[[2]];         A44[3:4,4:6] = prob[[1]];
    }
  }
  
  if (startP==3 && startC==5)   # seventh of a 2x2 block
  {
    # when the adjust was previously done
    if (is.na(F44[2,3,1])==FALSE)
    {
      poly = F44[2,3,];              values = A44[3:4,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F44[2,3,] = prob[[2]];         A44[3:4,7:9] = prob[[1]];
    }
  }
  
  if (startP==3 && startC==7)   # eighth of a 2x2 block
  {
    # when the adjust was previously done
    if (is.na(F44[2,4,1])==FALSE)
    {
      poly = F44[2,4,];              values = A44[3:4,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F44[2,4,] = prob[[2]];         A44[3:4,10:12] = prob[[1]];
    }
  }
  
  
  else if (startP==5 && startC==1)   # ninth of a 2x2 block
  {
    # when the adjust was previously done
    if (is.na(F44[3,1,1])==FALSE)
    {
      poly = F44[3,1,];              values = A44[5:6,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F44[3,1,] = prob[[2]];         A44[5:6,1:3] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==3)   # tenth of a 2x2 block
  {
    # when the adjust was previously done
    if (is.na(F44[3,2,1])==FALSE)
    {
      poly = F44[3,2,];              values = A44[5:6,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F44[3,2,] = prob[[2]];         A44[5:6,4:6] = prob[[1]];
    }
  }
  
  if (startP==5 && startC==5)   # eleventh of a 2x2 block
  {
    # when the adjust was previously done
    if (is.na(F44[3,3,1])==FALSE)
    {
      poly = F44[3,3,];              values = A44[5:6,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F44[3,3,] = prob[[2]];         A44[5:6,7:9] = prob[[1]];
    }
  }
  
  if (startP==5 && startC==7)   # twelfth of a 2x2 block
  {
    # when the adjust was previously done
    if (is.na(F44[3,4,1])==FALSE)
    {
      poly = F44[3,4,];              values = A44[5:6,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F44[3,4,] = prob[[2]];         A44[5:6,10:12] = prob[[1]];
    }
  }
  
  
  else if (startP==7 && startC==1)   # thirteenth of a 2x2 block
  {
    # when the adjust was previously done
    if (is.na(F44[4,1,1])==FALSE)
    {
      poly = F44[4,1,];              values = A44[7:8,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F44[4,1,] = prob[[2]];         A44[7:8,1:3] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==3)   # fourtheenth of a 2x2 block
  {
    # when the adjust was previously done
    if (is.na(F44[4,2,1])==FALSE)
    {
      poly = F44[4,2,];              values = A44[7:8,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F44[4,2,] = prob[[2]];         A44[7:8,4:6] = prob[[1]];
    }
  }
  
  if (startP==7 && startC==5)   # fifteenth of a 2x2 block
  {
    # when the adjust was previously done
    if (is.na(F44[4,3,1])==FALSE)
    {
      poly = F44[4,3,];              values = A44[7:8,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F44[4,3,] = prob[[2]];         A44[7:8,7:9] = prob[[1]];
    }
  }
  
  if (startP==7 && startC==7)   # sixteenth of a 2x2 block
  {
    # when the adjust was previously done
    if (is.na(F44[4,4,1])==FALSE)
    {
      poly = F44[4,4,];              values = A44[7:8,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F44[4,4,] = prob[[2]];         A44[7:8,10:12] = prob[[1]];
    }
  }
  
  return(list(poly,values,F44,A44));
}
############################################################







Cases2x8 = function(F28,A28,poly,values,startP,startC,d,dataMin,dataMax,n,k,maxDegree,
                    verbose=FALSE,plots=FALSE)
{
  # it checks if and adjust that has to be done was previously done and stored in
  # the matrices F28 and A28. In that case, matrices are updated
  # INPUT: - F28 = polynomial adjust of a discretization 2x8 (blocks of 4x1)
  #        - A28 = values of the adjust of a discretization 2x8
  #        - startP, startC = firs index of the block in each parent
  #        - d = selected data to adjust
  #        - dataMin, dataMax = the minimum and maximum value in the variable
  #        - n = position of the variable in the dataframe
  #        - k = position of the adjust in the matrix
  #        - maxDegree = maximum degree for the polynomial fit
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - poly = vector with the coefficients of the polynomial for d
  #         - values = matrix with 2 rows and 3 columns. It contains the
  #         min and max root, the min and max value and the probab. of the queues
  #         - F28 = array 2x8x11 for the coefficients of the polynomials
  #         - A28 = array 4x24 for the values of the adjust
  
  
  # create the vector for the adjust (simple polynomial adjust)
  r = rep(0,length(d));
  
  # check the 16 block where the adjust belongs
  if (startP==1 && startC==1)   # first of a 4x1 block
  {
    # when the adjust was previously done
    if (is.na(F28[1,1,1])==FALSE)
    {
      poly = F28[1,1,];              values = A28[1:2,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F28[1,1,] = prob[[2]];         A28[1:2,1:3] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==2)   # second of a 4x1 block
  {
    # when the adjust was previously done
    if (is.na(F28[1,2,1])==FALSE)
    {
      poly = F28[1,2,];              values = A28[1:2,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F28[1,2,] = prob[[2]];         A28[1:2,4:6] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==3)   # third of a 4x1 block
  {
    # when the adjust was previously done
    if (is.na(F28[1,3,1])==FALSE)
    {
      poly = F28[1,3,];              values = A28[1:2,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F28[1,3,] = prob[[2]];         A28[1:2,7:9] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==4)   # fourth of a 4x1 block
  {
    # when the adjust was previously done
    if (is.na(F28[1,4,1])==FALSE)
    {
      poly = F28[1,4,];              values = A28[1:2,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F28[1,4,] = prob[[2]];         A28[1:2,10:12] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==5)   # fifth of a 4x1 block
  {
    # when the adjust was previously done
    if (is.na(F28[1,5,1])==FALSE)
    {
      poly = F28[1,5,];              values = A28[1:2,13:15];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F28[1,5,] = prob[[2]];         A28[1:2,13:15] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==6)   # sixth of a 4x1 block
  {
    # when the adjust was previously done
    if (is.na(F28[1,6,1])==FALSE)
    {
      poly = F28[1,6,];              values = A28[1:2,16:18];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F28[1,6,] = prob[[2]];         A28[1:2,16:18] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==7)   # seventh of a 4x1 block
  {
    # when the adjust was previously done
    if (is.na(F28[1,7,1])==FALSE)
    {
      poly = F28[1,7,];              values = A28[1:2,19:21];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F28[1,7,] = prob[[2]];         A28[1:2,19:21] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==8)   # eighth of a 4x1 block
  {
    # when the adjust was previously done
    if (is.na(F28[1,8,1])==FALSE)
    {
      poly = F28[1,8,];              values = A28[1:2,22:24];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F28[1,8,] = prob[[2]];         A28[1:2,22:24] = prob[[1]];
    }
  }
  
  
  else if (startP==5 && startC==1)   # ninth of a 4x1 block
  {
    # when the adjust was previously done
    if (is.na(F28[2,1,1])==FALSE)
    {
      poly = F28[2,1,];              values = A28[3:4,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F28[2,1,] = prob[[2]];         A28[3:4,1:3] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==2)   # tenth of a 4x1 block
  {
    # when the adjust was previously done
    if (is.na(F28[2,2,1])==FALSE)
    {
      poly = F28[2,2,];              values = A28[3:4,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F28[2,2,] = prob[[2]];         A28[3:4,4:6] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==3)   # eleventh of a 4x1 block
  {
    # when the adjust was previously done
    if (is.na(F28[2,3,1])==FALSE)
    {
      poly = F28[2,3,];              values = A28[3:4,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F28[2,3,] = prob[[2]];         A28[3:4,7:9] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==4)   # twelfth of a 4x1 block
  {
    # when the adjust was previously done
    if (is.na(F28[2,4,1])==FALSE)
    {
      poly = F28[2,4,];              values = A28[3:4,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F28[2,4,] = prob[[2]];         A28[3:4,10:12] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==5)   # thrirteenth of a 4x1 block
  {
    # when the adjust was previously done
    if (is.na(F28[2,5,1])==FALSE)
    {
      poly = F28[2,5,];              values = A28[3:4,13:15];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F28[2,5,] = prob[[2]];         A28[3:4,13:15] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==6)   # fourteenth of a 4x1 block
  {
    # when the adjust was previously done
    if (is.na(F28[2,6,1])==FALSE)
    {
      poly = F28[2,6,];              values = A28[3:4,16:18];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F28[2,6,] = prob[[2]];         A28[3:4,16:18] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==7)   # fifteenth of a 4x1 block
  {
    # when the adjust was previously done
    if (is.na(F28[2,7,1])==FALSE)
    {
      poly = F28[2,7,];              values = A28[3:4,19:21];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F28[2,7,] = prob[[2]];         A28[3:4,19:21] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==8)   # sixteenth of a 4x1 block
  {
    # when the adjust was previously done
    if (is.na(F28[2,8,1])==FALSE)
    {
      poly = F28[2,8,];              values = A28[3:4,22:24];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F28[2,8,] = prob[[2]];         A28[3:4,22:24] = prob[[1]];
    }
  }
  
  return(list(poly,values,F28,A28));
}
############################################################



  




Cases8x2 = function(F82,A82,poly,values,startP,startC,d,dataMin,dataMax,n,k,maxDegree,
                    verbose=FALSE,plots=FALSE)
{
  # it checks if and adjust that has to be done was previously done and stored in
  # the matrices F82 and A82. In that case, matrices are updated
  # INPUT: - F82 = polynomial adjust of a discretization 8x2 (blocks of 1x4)
  #        - A82 = values of the adjust of a discretization 8x2
  #        - startP, startC = firs index of the block in each parent
  #        - d = selected data to adjust
  #        - dataMin, dataMax = the minimum and maximum value in the variable
  #        - n = position of the variable in the dataframe
  #        - k = position of the adjust in the matrix
  #        - maxDegree = maximum degree for the polynomial fit
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - poly = vector with the coefficients of the polynomial for d
  #         - values = matrix with 2 rows and 3 columns. It contains the
  #         min and max root, the min and max value and the probab. of the queues
  #         - F82 = array 8x2x11 for the coefficients of the polynomials
  #         - A82 = array 16x6 for the values of the adjust
  
  
  # create the vector for the adjust (simple polynomial adjust)
  r = rep(0,length(d));
  
  # check the 16 block where the adjust belongs
  if (startP==1 && startC==1)   # first of a 1x4 block
  {
    # when the adjust was previously done
    if (is.na(F82[1,1,1])==FALSE)
    {
      poly = F82[1,1,];              values = A82[1:2,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F82[1,1,] = prob[[2]];         A82[1:2,1:3] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==5)   # second of a 1x4 block
  {
    # when the adjust was previously done
    if (is.na(F82[1,2,1])==FALSE)
    {
      poly = F82[1,2,];              values = A82[1:2,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F82[1,2,] = prob[[2]];         A82[1:2,4:6] = prob[[1]];
    }
  }
  
  
  else if (startP==2 && startC==1)   # third of a 1x4 block
  {
    # when the adjust was previously done
    if (is.na(F82[2,1,1])==FALSE)
    {
      poly = F82[2,1,];              values = A82[3:4,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F82[2,1,] = prob[[2]];         A82[3:4,1:3] = prob[[1]];
    }
  }
  
  else if (startP==2 && startC==5)   # fourth of a 1x4 block
  {
    # when the adjust was previously done
    if (is.na(F82[2,2,1])==FALSE)
    {
      poly = F82[2,2,];              values = A82[3:4,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F82[2,2,] = prob[[2]];         A82[3:4,4:6] = prob[[1]];
    }
  }
  
  
  else if (startP==3 && startC==1)   # fifth of a 1x4 block
  {
    # when the adjust was previously done
    if (is.na(F82[3,1,1])==FALSE)
    {
      poly = F82[3,1,];              values = A82[5:6,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F82[3,1,] = prob[[2]];         A82[5:6,1:3] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==5)   # sixth of a 1x4 block
  {
    # when the adjust was previously done
    if (is.na(F82[3,2,1])==FALSE)
    {
      poly = F82[3,2,];              values = A82[5:6,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F82[3,2,] = prob[[2]];         A82[5:6,4:6] = prob[[1]];
    }
  }
  
  
  else if (startP==4 && startC==1)   # seventh of a 1x4 block
  {
    # when the adjust was previously done
    if (is.na(F82[4,1,1])==FALSE)
    {
      poly = F82[4,1,];              values = A82[7:8,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F82[4,1,] = prob[[2]];         A82[7:8,1:3] = prob[[1]];
    }
  }
  
  else if (startP==4 && startC==5)   # eighth of a 1x4 block
  {
    # when the adjust was previously done
    if (is.na(F82[4,2,1])==FALSE)
    {
      poly = F82[4,2,];              values = A82[7:8,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F82[4,2,] = prob[[2]];         A82[7:8,4:6] = prob[[1]];
    }
  }
  
  
  else if (startP==5 && startC==1)   # ninth of a 1x4 block
  {
    # when the adjust was previously done
    if (is.na(F82[5,1,1])==FALSE)
    {
      poly = F82[5,1,];              values = A82[9:10,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F82[5,1,] = prob[[2]];         A82[9:10,1:3] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==5)   # tenth of a 1x4 block
  {
    # when the adjust was previously done
    if (is.na(F82[5,2,1])==FALSE)
    {
      poly = F82[5,2,];              values = A82[9:10,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F82[5,2,] = prob[[2]];         A82[9:10,4:6] = prob[[1]];
    }
  }
  
  
  else if (startP==6 && startC==1)   # eleventh of a 1x4 block
  {
    # when the adjust was previously done
    if (is.na(F82[6,1,1])==FALSE)
    {
      poly = F82[6,1,];              values = A82[11:12,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F82[6,1,] = prob[[2]];         A82[11:12,1:3] = prob[[1]];
    }
  }
  
  else if (startP==6 && startC==5)   # twelfth of a 1x4 block
  {
    # when the adjust was previously done
    if (is.na(F82[6,2,1])==FALSE)
    {
      poly = F82[6,2,];              values = A82[11:12,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F82[6,2,] = prob[[2]];         A82[11:12,4:6] = prob[[1]];
    }
  }
  
  
  else if (startP==7 && startC==1)   # thirteenth of a 1x4 block
  {
    # when the adjust was previously done
    if (is.na(F82[7,1,1])==FALSE)
    {
      poly = F82[7,1,];              values = A82[13:14,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F82[7,1,] = prob[[2]];         A82[13:14,1:3] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==5)   # fourteenth of a 1x4 block
  {
    # when the adjust was previously done
    if (is.na(F82[7,2,1])==FALSE)
    {
      poly = F82[7,2,];              values = A82[13:14,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F82[7,2,] = prob[[2]];         A82[13:14,4:6] = prob[[1]];
    }
  }
  
  
  else if (startP==8 && startC==1)   # fifteenth of a 1x4 block
  {
    # when the adjust was previously done
    if (is.na(F82[8,1,1])==FALSE)
    {
      poly = F82[8,1,];              values = A82[15:16,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F82[8,1,] = prob[[2]];         A82[15:16,1:3] = prob[[1]];
    }
  }
  
  else if (startP==8 && startC==5)   # sixteenth of a 1x4 block
  {
    # when the adjust was previously done
    if (is.na(F82[8,2,1])==FALSE)
    {
      poly = F82[8,2,];              values = A82[15:16,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F82[8,2,] = prob[[2]];         A82[15:16,4:6] = prob[[1]];
    }
  }
  
  return(list(poly,values,F82,A82));
}
############################################################



  





Cases4x8 = function(F48,A48,poly,values,startP,startC,d,dataMin,dataMax,n,k,maxDegree,
                    verbose=FALSE,plots=FALSE)
{
  # it checks if and adjust that has to be done was previously done and stored in
  # the matrices F48 and A48. In that case, matrices are updated
  # INPUT: - F48 = polynomial adjust of a discretization 4x8 (blocks of 2x1)
  #        - A48 = values of the adjust of a discretization 4x8
  #        - startP, startC = firs index of the block in each parent
  #        - d = selected data to adjust
  #        - dataMin, dataMax = the minimum and maximum value in the variable
  #        - n = position of the variable in the dataframe
  #        - k = position of the adjust in the matrix
  #        - maxDegree = maximum degree for the polynomial fit
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - poly = vector with the coefficients of the polynomial for d
  #         - values = matrix with 2 rows and 3 columns. It contains the
  #         min and max root, the min and max value and the probab. of the queues
  #         - F48 = array 4x8x11 for the coefficients of the polynomials
  #         - A48 = array 8x24 for the values of the adjust
  
  
  # create the vector for the adjust (simple polynomial adjust)
  r = rep(0,length(d));
  
  # check the 32 block where the adjust belongs
  if (startP==1 && startC==1)   # first of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[1,1,1])==FALSE)
    {
      poly = F48[1,1,];              values = A48[1:2,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[1,1,] = prob[[2]];         A48[1:2,1:3] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==2)   # second of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[1,2,1])==FALSE)
    {
      poly = F48[1,2,];              values = A48[1:2,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[1,2,] = prob[[2]];         A48[1:2,4:6] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==3)   # third of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[1,3,1])==FALSE)
    {
      poly = F48[1,3,];              values = A48[1:2,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[1,3,] = prob[[2]];         A48[1:2,7:9] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==4)   # fourth of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[1,4,1])==FALSE)
    {
      poly = F48[1,4,];              values = A48[1:2,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[1,4,] = prob[[2]];         A48[1:2,10:12] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==5)   # fifth of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[1,5,1])==FALSE)
    {
      poly = F48[1,5,];              values = A48[1:2,13:15];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[1,5,] = prob[[2]];         A48[1:2,13:15] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==6)   # sixth of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[1,6,1])==FALSE)
    {
      poly = F48[1,6,];              values = A48[1:2,16:18];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[1,6,] = prob[[2]];         A48[1:2,16:18] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==7)   # seventh of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[1,7,1])==FALSE)
    {
      poly = F48[1,7,];              values = A48[1:2,19:21];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[1,7,] = prob[[2]];         A48[1:2,19:21] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==8)   # eighth of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[1,8,1])==FALSE)
    {
      poly = F48[1,8,];              values = A48[1:2,22:24];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[1,8,] = prob[[2]];         A48[1:2,22:24] = prob[[1]];
    }
  }
  
  
  else if (startP==3 && startC==1)   # ninth of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[2,1,1])==FALSE)
    {
      poly = F48[2,1,];              values = A48[3:4,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[2,1,] = prob[[2]];         A48[3:4,1:3] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==2)   # tenth of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[2,2,1])==FALSE)
    {
      poly = F48[2,2,];              values = A48[3:4,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[2,2,] = prob[[2]];         A48[3:4,4:6] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==3)   # eleventh of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[2,3,1])==FALSE)
    {
      poly = F48[2,3,];              values = A48[3:4,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[2,3,] = prob[[2]];         A48[3:4,7:9] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==4)   # twelfth of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[2,4,1])==FALSE)
    {
      poly = F48[2,4,];              values = A48[3:4,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[2,4,] = prob[[2]];         A48[3:4,10:12] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==5)   # thirteenth of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[2,5,1])==FALSE)
    {
      poly = F48[2,5,];              values = A48[3:4,13:15];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[2,5,] = prob[[2]];         A48[3:4,13:15] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==6)   # fourteenth of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[2,6,1])==FALSE)
    {
      poly = F48[2,6,];              values = A48[3:4,16:18];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[2,6,] = prob[[2]];         A48[3:4,16:18] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==7)   # fifteenth of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[2,7,1])==FALSE)
    {
      poly = F48[2,7,];              values = A48[3:4,19:21];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[2,7,] = prob[[2]];         A48[3:4,19:21] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==8)   # sixteenth of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[2,8,1])==FALSE)
    {
      poly = F48[2,8,];              values = A48[3:4,22:24];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[2,8,] = prob[[2]];         A48[3:4,22:24] = prob[[1]];
    }
  }
  
  
  else if (startP==5 && startC==1)   # seventeenth of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[3,1,1])==FALSE)
    {
      poly = F48[3,1,];              values = A48[5:6,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[3,1,] = prob[[2]];         A48[5:6,1:3] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==2)   # eighteenth of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[3,2,1])==FALSE)
    {
      poly = F48[3,2,];              values = A48[5:6,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[3,2,] = prob[[2]];         A48[5:6,4:6] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==3)   # nineteenth of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[3,3,1])==FALSE)
    {
      poly = F48[3,3,];              values = A48[5:6,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[3,3,] = prob[[2]];         A48[5:6,7:9] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==4)   # 20th of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[3,4,1])==FALSE)
    {
      poly = F48[3,4,];              values = A48[5:6,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[3,4,] = prob[[2]];         A48[5:6,10:12] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==5)   # 21st of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[3,5,1])==FALSE)
    {
      poly = F48[3,5,];              values = A48[5:6,13:15];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[3,5,] = prob[[2]];         A48[5:6,13:15] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==6)   # 22nd of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[3,6,1])==FALSE)
    {
      poly = F48[3,6,];              values = A48[5:6,16:18];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[3,6,] = prob[[2]];         A48[5:6,16:18] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==7)   # 23th of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[3,7,1])==FALSE)
    {
      poly = F48[3,7,];              values = A48[5:6,19:21];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[3,7,] = prob[[2]];         A48[5:6,19:21] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==8)   # 24th of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[3,8,1])==FALSE)
    {
      poly = F48[3,8,];              values = A48[5:6,22:24];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[3,8,] = prob[[2]];         A48[5:6,22:24] = prob[[1]];
    }
  }
  
  
  if (startP==7 && startC==1)   # 25th of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[4,1,1])==FALSE)
    {
      poly = F48[4,1,];              values = A48[7:8,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[4,1,] = prob[[2]];         A48[7:8,1:3] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==2)   # 26th of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[4,2,1])==FALSE)
    {
      poly = F48[4,2,];              values = A48[7:8,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[4,2,] = prob[[2]];         A48[7:8,4:6] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==3)   # 27th of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[4,3,1])==FALSE)
    {
      poly = F48[4,3,];              values = A48[7:8,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[4,3,] = prob[[2]];         A48[7:8,7:9] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==4)   # 28th of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[4,4,1])==FALSE)
    {
      poly = F48[4,4,];              values = A48[7:8,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[4,4,] = prob[[2]];         A48[7:8,10:12] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==5)   # 29th of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[4,5,1])==FALSE)
    {
      poly = F48[4,5,];              values = A48[7:8,13:15];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[4,5,] = prob[[2]];         A48[7:8,13:15] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==6)   # 30th of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[4,6,1])==FALSE)
    {
      poly = F48[4,6,];              values = A48[7:8,16:18];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[4,6,] = prob[[2]];         A48[7:8,16:18] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==7)   # 31st of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[4,7,1])==FALSE)
    {
      poly = F48[4,7,];              values = A48[7:8,19:21];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[4,7,] = prob[[2]];         A48[7:8,19:21] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==8)   # 32nd of a 2x1 block
  {
    # when the adjust was previously done
    if (is.na(F48[4,8,1])==FALSE)
    {
      poly = F48[4,8,];              values = A48[7:8,22:24];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F48[4,8,] = prob[[2]];         A48[7:8,22:24] = prob[[1]];
    }
  }
  
  return(list(poly,values,F48,A48));
}
############################################################







Cases8x4 = function(F84,A84,poly,values,startP,startC,d,dataMin,dataMax,n,k,maxDegree,
                    verbose=FALSE,plots=FALSE)
{
  # it checks if and adjust that has to be done was previously done and stored in
  # the matrices F84 and A84. In that case, matrices are updated
  # INPUT: - F84 = polynomial adjust of a discretization 8x4 (blocks of 1x2)
  #        - A84 = values of the adjust of a discretization 8x4
  #        - startP, startC = firs index of the block in each parent
  #        - d = selected data to adjust
  #        - dataMin, dataMax = the minimum and maximum value in the variable
  #        - n = position of the variable in the dataframe
  #        - k = position of the adjust in the matrix
  #        - maxDegree = maximum degree for the polynomial adjust
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - poly = vector with the coefficients of the polynomial for d
  #         - values = matrix with 2 rows and 3 columns. It contains the
  #         min and max root, the min and max value and the probab. of the queues
  #         - F84 = array 8x4x11 for the coefficients of the polynomials
  #         - A84 = array 16x12 for the values of the adjust
  
  
  # create the vector for the adjust (simple polynomial adjust)
  r = rep(0,length(d));
  
  # check the 32 block where the adjust belongs
  if (startP==1 && startC==1)   # first of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[1,1,1])==FALSE)
    {
      poly = F84[1,1,];              values = A84[1:2,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[1,1,] = prob[[2]];         A84[1:2,1:3] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==3)   # second of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[1,2,1])==FALSE)
    {
      poly = F84[1,2,];              values = A84[1:2,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[1,2,] = prob[[2]];         A84[1:2,4:6] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==5)   # third of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[1,3,1])==FALSE)
    {
      poly = F84[1,3,];              values = A84[1:2,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[1,3,] = prob[[2]];         A84[1:2,7:9] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==7)   # fourth of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[1,4,1])==FALSE)
    {
      poly = F84[1,4,];              values = A84[1:2,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[1,4,] = prob[[2]];         A84[1:2,10:12] = prob[[1]];
    }
  }
  
  
  if (startP==2 && startC==1)   # fifth of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[2,1,1])==FALSE)
    {
      poly = F84[2,1,];              values = A84[3:4,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[2,1,] = prob[[2]];         A84[3:4,1:3] = prob[[1]];
    }
  }
  
  else if (startP==2 && startC==3)   # sixth of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[2,2,1])==FALSE)
    {
      poly = F84[2,2,];              values = A84[3:4,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[2,2,] = prob[[2]];         A84[3:4,4:6] = prob[[1]];
    }
  }
  
  else if (startP==2 && startC==5)   # seventh of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[2,3,1])==FALSE)
    {
      poly = F84[2,3,];              values = A84[3:4,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[2,3,] = prob[[2]];         A84[3:4,7:9] = prob[[1]];
    }
  }
  
  else if (startP==2 && startC==7)   # eighth of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[2,4,1])==FALSE)
    {
      poly = F84[2,4,];              values = A84[3:4,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[2,4,] = prob[[2]];         A84[3:4,10:12] = prob[[1]];
    }
  }
  
  
  if (startP==3 && startC==1)   # ninth of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[3,1,1])==FALSE)
    {
      poly = F84[3,1,];              values = A84[5:6,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[3,1,] = prob[[2]];         A84[5:6,1:3] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==3)   # tenth of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[3,2,1])==FALSE)
    {
      poly = F84[3,2,];              values = A84[5:6,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[3,2,] = prob[[2]];         A84[5:6,4:6] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==5)   # eleventh of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[3,3,1])==FALSE)
    {
      poly = F84[3,3,];              values = A84[5:6,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[3,3,] = prob[[2]];         A84[5:6,7:9] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==7)   # twelfth of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[3,4,1])==FALSE)
    {
      poly = F84[3,4,];              values = A84[5:6,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[3,4,] = prob[[2]];         A84[5:6,10:12] = prob[[1]];
    }
  }
  
  
  if (startP==4 && startC==1)   # thirteenth of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[4,1,1])==FALSE)
    {
      poly = F84[4,1,];              values = A84[7:8,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[4,1,] = prob[[2]];         A84[7:8,1:3] = prob[[1]];
    }
  }
  
  else if (startP==4 && startC==3)   # fourteenth of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[4,2,1])==FALSE)
    {
      poly = F84[4,2,];              values = A84[7:8,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[4,2,] = prob[[2]];         A84[7:8,4:6] = prob[[1]];
    }
  }
  
  else if (startP==4 && startC==5)   # fifteenth of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[4,3,1])==FALSE)
    {
      poly = F84[4,3,];              values = A84[7:8,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[4,3,] = prob[[2]];         A84[7:8,7:9] = prob[[1]];
    }
  }
  
  else if (startP==4 && startC==7)   # sixteenth of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[4,4,1])==FALSE)
    {
      poly = F84[4,4,];              values = A84[7:8,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[4,4,] = prob[[2]];         A84[7:8,10:12] = prob[[1]];
    }
  }
  
  
  if (startP==5 && startC==1)   # seventeenth of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[5,1,1])==FALSE)
    {
      poly = F84[5,1,];              values = A84[9:10,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[5,1,] = prob[[2]];         A84[9:10,1:3] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==3)   # eighteenth of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[5,2,1])==FALSE)
    {
      poly = F84[5,2,];              values = A84[9:10,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[5,2,] = prob[[2]];         A84[9:10,4:6] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==5)   # nineteenth of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[5,3,1])==FALSE)
    {
      poly = F84[5,3,];              values = A84[9:10,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[5,3,] = prob[[2]];         A84[9:10,7:9] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==7)   # 20th of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[5,4,1])==FALSE)
    {
      poly = F84[5,4,];              values = A84[9:10,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[5,4,] = prob[[2]];         A84[9:10,10:12] = prob[[1]];
    }
  }
  
  
  if (startP==6 && startC==1)   # 21st of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[6,1,1])==FALSE)
    {
      poly = F84[6,1,];              values = A84[11:12,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[6,1,] = prob[[2]];         A84[11:12,1:3] = prob[[1]];
    }
  }
  
  else if (startP==6 && startC==3)   # 22nd of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[6,2,1])==FALSE)
    {
      poly = F84[6,2,];              values = A84[11:12,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[6,2,] = prob[[2]];         A84[11:12,4:6] = prob[[1]];
    }
  }
  
  else if (startP==6 && startC==5)   # 23th of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[6,3,1])==FALSE)
    {
      poly = F84[6,3,];              values = A84[11:12,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[6,3,] = prob[[2]];         A84[11:12,7:9] = prob[[1]];
    }
  }
  
  else if (startP==6 && startC==7)   # 24th of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[6,4,1])==FALSE)
    {
      poly = F84[6,4,];              values = A84[11:12,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[6,4,] = prob[[2]];         A84[11:12,10:12] = prob[[1]];
    }
  }
  
  
  if (startP==7 && startC==1)   # 25th of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[7,1,1])==FALSE)
    {
      poly = F84[7,1,];              values = A84[13:14,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[7,1,] = prob[[2]];         A84[13:14,1:3] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==3)   # 26th of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[7,2,1])==FALSE)
    {
      poly = F84[7,2,];              values = A84[13:14,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[7,2,] = prob[[2]];         A84[13:14,4:6] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==5)   # 27th of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[7,3,1])==FALSE)
    {
      poly = F84[7,3,];              values = A84[13:14,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[7,3,] = prob[[2]];         A84[13:14,7:9] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==7)   # 28th of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[7,4,1])==FALSE)
    {
      poly = F84[7,4,];              values = A84[13:14,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[7,4,] = prob[[2]];         A84[13:14,10:12] = prob[[1]];
    }
  }
  
  
  if (startP==8 && startC==1)   # 29th of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[8,1,1])==FALSE)
    {
      poly = F84[8,1,];              values = A84[15:16,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[8,1,] = prob[[2]];         A84[15:16,1:3] = prob[[1]];
    }
  }
  
  else if (startP==8 && startC==3)   # 30th of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[8,2,1])==FALSE)
    {
      poly = F84[8,2,];              values = A84[15:16,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[8,2,] = prob[[2]];         A84[15:16,4:6] = prob[[1]];
    }
  }
  
  else if (startP==8 && startC==5)   # 31st of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[8,3,1])==FALSE)
    {
      poly = F84[8,3,];              values = A84[15:16,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[8,3,] = prob[[2]];         A84[15:16,7:9] = prob[[1]];
    }
  }
  
  else if (startP==8 && startC==7)   # 32nd of a 1x2 block
  {
    # when the adjust was previously done
    if (is.na(F84[8,4,1])==FALSE)
    {
      poly = F84[8,4,];              values = A84[15:16,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F84[8,4,] = prob[[2]];         A84[15:16,10:12] = prob[[1]];
    }
  }
  
  return(list(poly,values,F84,A84));
}
############################################################









Cases8x8 = function(F88,A88,poly,values,startP,startC,d,dataMin,dataMax,n,k,maxDegree,
                    verbose=FALSE,plots=FALSE)
{
  # it checks if and adjust that has to be done was previously done and stored in
  # the matrices F88 and A88. In that case, matrices are updated
  # INPUT: - F88 = polynomial adjust of a discretization 8x8 (blocks of 1x1)
  #        - A88 = values of the adjust of a discretization 8x8
  #        - startP, startC = firs index of the block in each parent
  #        - d = selected data to adjust
  #        - dataMin, dataMax = the minimum and maximum value in the variable
  #        - n = position of the variable in the dataframe
  #        - k = position of the adjust in the matrix
  #        - maxDegree = maximum degree for the polynomial fit
  #        - verbose = show the 'cats' on screen
  #        - plots = draw the graphs
  # OUTPUT: - poly = vector with the coefficients of the polynomial for d
  #         - values = matrix with 2 rows and 3 columns. It contains the
  #         min and max root, the min and max value and the probab. of the queues
  #         - F88 = array 8x8x11 for the coefficients of the polynomials
  #         - A88 = array 16x24 for the values of the adjust
  
  
  # create the vector for the adjust (simple polynomial adjust)
  r = rep(0,length(d));
  
  # check the 64 block where the adjust belongs
  if (startP==1 && startC==1)   # first of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[1,1,1])==FALSE)
    {
      poly = F88[1,1,];              values = A88[1:2,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[1,1,] = prob[[2]];         A88[1:2,1:3] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==2)   # second of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[1,2,1])==FALSE)
    {
      poly = F88[1,2,];              values = A88[1:2,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[1,2,] = prob[[2]];         A88[1:2,4:6] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==3)   # third of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[1,3,1])==FALSE)
    {
      poly = F88[1,3,];              values = A88[1:2,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[1,3,] = prob[[2]];         A88[1:2,7:9] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==4)   # fourth of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[1,4,1])==FALSE)
    {
      poly = F88[1,4,];              values = A88[1:2,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[1,4,] = prob[[2]];         A88[1:2,10:12] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==5)   # fifth of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[1,5,1])==FALSE)
    {
      poly = F88[1,5,];              values = A88[1:2,13:15];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[1,5,] = prob[[2]];         A88[1:2,13:15] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==6)   # sixth of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[1,6,1])==FALSE)
    {
      poly = F88[1,6,];              values = A88[1:2,16:18];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[1,6,] = prob[[2]];         A88[1:2,16:18] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==7)   # seventh of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[1,7,1])==FALSE)
    {
      poly = F88[1,7,];              values = A88[1:2,19:21];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[1,7,] = prob[[2]];         A88[1:2,19:21] = prob[[1]];
    }
  }
  
  else if (startP==1 && startC==8)   # eighth of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[1,8,1])==FALSE)
    {
      poly = F88[1,8,];              values = A88[1:2,22:24];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[1,8,] = prob[[2]];         A88[1:2,22:24] = prob[[1]];
    }
  }
  
    
  else if (startP==2 && startC==1)   # ninth of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[2,1,1])==FALSE)
    {
      poly = F88[2,1,];              values = A88[3:4,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[2,1,] = prob[[2]];         A88[3:4,1:3] = prob[[1]];
    }
  }
  
  else if (startP==2 && startC==2)   # tenth of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[2,2,1])==FALSE)
    {
      poly = F88[2,2,];              values = A88[3:4,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[2,2,] = prob[[2]];         A88[3:4,4:6] = prob[[1]];
    }
  }
  
  else if (startP==2 && startC==3)   # eleventh of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[2,3,1])==FALSE)
    {
      poly = F88[2,3,];              values = A88[3:4,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[2,3,] = prob[[2]];         A88[3:4,7:9] = prob[[1]];
    }
  }
  
  else if (startP==2 && startC==4)   # twelfth of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[2,4,1])==FALSE)
    {
      poly = F88[2,4,];              values = A88[3:4,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[2,4,] = prob[[2]];         A88[3:4,10:12] = prob[[1]];
    }
  }
  
  else if (startP==2 && startC==5)   # thirteenth of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[2,5,1])==FALSE)
    {
      poly = F88[2,5,];              values = A88[3:4,13:15];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[2,5,] = prob[[2]];         A88[3:4,13:15] = prob[[1]];
    }
  }
  
  else if (startP==2 && startC==6)   # fourteenth of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[2,6,1])==FALSE)
    {
      poly = F88[2,6,];              values = A88[3:4,16:18];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[2,6,] = prob[[2]];         A88[3:4,16:18] = prob[[1]];
    }
  }
  
  else if (startP==2 && startC==7)   # fifteenth of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[2,7,1])==FALSE)
    {
      poly = F88[2,7,];              values = A88[3:4,19:21];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[2,7,] = prob[[2]];         A88[3:4,19:21] = prob[[1]];
    }
  }
  
  else if (startP==2 && startC==8)   # sixteenth of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[2,8,1])==FALSE)
    {
      poly = F88[2,8,];              values = A88[3:4,22:24];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[2,8,] = prob[[2]];         A88[3:4,22:24] = prob[[1]];
    }
  }
  
  
  else if (startP==3 && startC==1)   # seventeenth of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[3,1,1])==FALSE)
    {
      poly = F88[3,1,];              values = A88[5:6,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[3,1,] = prob[[2]];         A88[5:6,1:3] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==2)   # eighteenth of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[3,2,1])==FALSE)
    {
      poly = F88[3,2,];              values = A88[5:6,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[3,2,] = prob[[2]];         A88[5:6,4:6] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==3)   # nineteenth of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[3,3,1])==FALSE)
    {
      poly = F88[3,3,];              values = A88[5:6,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[3,3,] = prob[[2]];         A88[5:6,7:9] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==4)   # 20th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[3,4,1])==FALSE)
    {
      poly = F88[3,4,];              values = A88[5:6,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[3,4,] = prob[[2]];         A88[5:6,10:12] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==5)   # 21st of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[3,5,1])==FALSE)
    {
      poly = F88[3,5,];              values = A88[5:6,13:15];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[3,5,] = prob[[2]];         A88[5:6,13:15] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==6)   # 22nd of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[3,6,1])==FALSE)
    {
      poly = F88[3,6,];              values = A88[5:6,16:18];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[3,6,] = prob[[2]];         A88[5:6,16:18] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==7)   # 23rd of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[3,7,1])==FALSE)
    {
      poly = F88[3,7,];              values = A88[5:6,19:21];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[3,7,] = prob[[2]];         A88[5:6,19:21] = prob[[1]];
    }
  }
  
  else if (startP==3 && startC==8)   # 24th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[3,8,1])==FALSE)
    {
      poly = F88[3,8,];              values = A88[5:6,22:24];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[3,8,] = prob[[2]];         A88[5:6,22:24] = prob[[1]];
    }
  }
  
  
  else if (startP==4 && startC==1)   # 25th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[4,1,1])==FALSE)
    {
      poly = F88[4,1,];              values = A88[7:8,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[4,1,] = prob[[2]];         A88[7:8,1:3] = prob[[1]];
    }
  }
  
  else if (startP==4 && startC==2)   # 26th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[4,2,1])==FALSE)
    {
      poly = F88[4,2,];              values = A88[7:8,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[4,2,] = prob[[2]];         A88[7:8,4:6] = prob[[1]];
    }
  }
  
  else if (startP==4 && startC==3)   # 27th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[4,3,1])==FALSE)
    {
      poly = F88[4,3,];              values = A88[7:8,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[4,3,] = prob[[2]];         A88[7:8,7:9] = prob[[1]];
    }
  }
  
  else if (startP==4 && startC==4)   # 28th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[4,4,1])==FALSE)
    {
      poly = F88[4,4,];              values = A88[7:8,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[4,4,] = prob[[2]];         A88[7:8,10:12] = prob[[1]];
    }
  }
  
  else if (startP==4 && startC==5)   # 29th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[4,5,1])==FALSE)
    {
      poly = F88[4,5,];              values = A88[7:8,13:15];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[4,5,] = prob[[2]];         A88[7:8,13:15] = prob[[1]];
    }
  }
  
  else if (startP==4 && startC==6)   # 30th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[4,6,1])==FALSE)
    {
      poly = F88[4,6,];              values = A88[7:8,16:18];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[4,6,] = prob[[2]];         A88[7:8,16:18] = prob[[1]];
    }
  }
  
  else if (startP==4 && startC==7)   # 31st of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[4,7,1])==FALSE)
    {
      poly = F88[4,7,];              values = A88[7:8,19:21];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[4,7,] = prob[[2]];         A88[7:8,19:21] = prob[[1]];
    }
  }
  
  else if (startP==4 && startC==8)   # 32nd of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[4,8,1])==FALSE)
    {
      poly = F88[4,8,];              values = A88[7:8,22:24];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[4,8,] = prob[[2]];         A88[7:8,22:24] = prob[[1]];
    }
  }
  
  
  else if (startP==5 && startC==1)   # 33rd of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[5,1,1])==FALSE)
    {
      poly = F88[5,1,];              values = A88[9:10,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[5,1,] = prob[[2]];         A88[9:10,1:3] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==2)   # 34th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[5,2,1])==FALSE)
    {
      poly = F88[5,2,];              values = A88[9:10,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[5,2,] = prob[[2]];         A88[9:10,4:6] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==3)   # 35th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[5,3,1])==FALSE)
    {
      poly = F88[5,3,];              values = A88[9:10,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[5,3,] = prob[[2]];         A88[9:10,7:9] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==4)   # 36th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[5,4,1])==FALSE)
    {
      poly = F88[5,4,];              values = A88[9:10,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[5,4,] = prob[[2]];         A88[9:10,10:12] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==5)   # 37th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[5,5,1])==FALSE)
    {
      poly = F88[5,5,];              values = A88[9:10,13:15];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[5,5,] = prob[[2]];         A88[9:10,13:15] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==6)   # 38th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[5,6,1])==FALSE)
    {
      poly = F88[5,6,];              values = A88[9:10,16:18];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[5,6,] = prob[[2]];         A88[9:10,16:18] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==7)   # 39th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[5,7,1])==FALSE)
    {
      poly = F88[5,7,];              values = A88[9:10,19:21];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[5,7,] = prob[[2]];         A88[9:10,19:21] = prob[[1]];
    }
  }
  
  else if (startP==5 && startC==8)   # 40th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[5,8,1])==FALSE)
    {
      poly = F88[5,8,];              values = A88[9:10,22:24];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[5,8,] = prob[[2]];         A88[9:10,22:24] = prob[[1]];
    }
  }
  
  
  else if (startP==6 && startC==1)   # 41st of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[6,1,1])==FALSE)
    {
      poly = F88[6,1,];              values = A88[11:12,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[6,1,] = prob[[2]];         A88[11:12,1:3] = prob[[1]];
    }
  }
  
  else if (startP==6 && startC==2)   # 42nd of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[6,2,1])==FALSE)
    {
      poly = F88[6,2,];              values = A88[11:12,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[6,2,] = prob[[2]];         A88[11:12,4:6] = prob[[1]];
    }
  }
  
  else if (startP==6 && startC==3)   # 43rd of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[6,3,1])==FALSE)
    {
      poly = F88[6,3,];              values = A88[11:12,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[6,3,] = prob[[2]];         A88[11:12,7:9] = prob[[1]];
    }
  }
  
  else if (startP==6 && startC==4)   # 44th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[6,4,1])==FALSE)
    {
      poly = F88[6,4,];              values = A88[11:12,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[6,4,] = prob[[2]];         A88[11:12,10:12] = prob[[1]];
    }
  }
  
  else if (startP==6 && startC==5)   # 45th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[6,5,1])==FALSE)
    {
      poly = F88[6,5,];              values = A88[11:12,13:15];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[6,5,] = prob[[2]];         A88[11:12,13:15] = prob[[1]];
    }
  }
  
  else if (startP==6 && startC==6)   # 46th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[6,6,1])==FALSE)
    {
      poly = F88[6,6,];              values = A88[11:12,16:18];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[6,6,] = prob[[2]];         A88[11:12,16:18] = prob[[1]];
    }
  }
  
  else if (startP==6 && startC==7)   # 47th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[6,7,1])==FALSE)
    {
      poly = F88[6,7,];              values = A88[11:12,19:21];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[6,7,] = prob[[2]];         A88[11:12,19:21] = prob[[1]];
    }
  }
  
  else if (startP==6 && startC==8)   # 48th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[6,8,1])==FALSE)
    {
      poly = F88[6,8,];              values = A88[11:12,22:24];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[6,8,] = prob[[2]];         A88[11:12,22:24] = prob[[1]];
    }
  }
  
  
  else if (startP==7 && startC==1)   # 49th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[7,1,1])==FALSE)
    {
      poly = F88[7,1,];              values = A88[13:14,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[7,1,] = prob[[2]];         A88[13:14,1:3] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==2)   # 50th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[7,2,1])==FALSE)
    {
      poly = F88[7,2,];              values = A88[13:14,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[7,2,] = prob[[2]];         A88[13:14,4:6] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==3)   # 51st of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[7,3,1])==FALSE)
    {
      poly = F88[7,3,];              values = A88[13:14,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[7,3,] = prob[[2]];         A88[13:14,7:9] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==4)   # 52nd of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[7,4,1])==FALSE)
    {
      poly = F88[7,4,];              values = A88[13:14,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[7,4,] = prob[[2]];         A88[13:14,10:12] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==5)   # 53rd of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[7,5,1])==FALSE)
    {
      poly = F88[7,5,];              values = A88[13:14,13:15];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[7,5,] = prob[[2]];         A88[13:14,13:15] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==6)   # 54th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[7,6,1])==FALSE)
    {
      poly = F88[7,6,];              values = A88[13:14,16:18];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[7,6,] = prob[[2]];         A88[13:14,16:18] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==7)   # 55st of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[7,7,1])==FALSE)
    {
      poly = F88[7,7,];              values = A88[13:14,19:21];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[7,7,] = prob[[2]];         A88[13:14,19:21] = prob[[1]];
    }
  }
  
  else if (startP==7 && startC==8)   # 56th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[7,8,1])==FALSE)
    {
      poly = F88[7,8,];              values = A88[13:14,22:24];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[7,8,] = prob[[2]];         A88[13:14,22:24] = prob[[1]];
    }
  }
  
  
  else if (startP==8 && startC==1)   # 57th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[8,1,1])==FALSE)
    {
      poly = F88[8,1,];              values = A88[15:16,1:3];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[8,1,] = prob[[2]];         A88[15:16,1:3] = prob[[1]];
    }
  }
  
  else if (startP==8 && startC==2)   # 58th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[8,2,1])==FALSE)
    {
      poly = F88[8,2,];              values = A88[15:16,4:6];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[8,2,] = prob[[2]];         A88[15:16,4:6] = prob[[1]];
    }
  }
  
  else if (startP==8 && startC==3)   # 59th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[8,3,1])==FALSE)
    {
      poly = F88[8,3,];              values = A88[15:16,7:9];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[8,3,] = prob[[2]];         A88[15:16,7:9] = prob[[1]];
    }
  }
  
  else if (startP==8 && startC==4)   # 60th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[8,4,1])==FALSE)
    {
      poly = F88[8,4,];              values = A88[15:16,10:12];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[8,4,] = prob[[2]];         A88[15:16,10:12] = prob[[1]];
    }
  }
  
  else if (startP==8 && startC==5)   # 61st of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[8,5,1])==FALSE)
    {
      poly = F88[8,5,];              values = A88[15:16,13:15];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[8,5,] = prob[[2]];         A88[15:16,13:15] = prob[[1]];
    }
  }
  
  else if (startP==8 && startC==6)   # 62nd of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[8,6,1])==FALSE)
    {
      poly = F88[8,6,];              values = A88[15:16,16:18];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[8,6,] = prob[[2]];         A88[15:16,16:18] = prob[[1]];
    }
  }
  
  else if (startP==8 && startC==7)   # 63rd of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[8,7,1])==FALSE)
    {
      poly = F88[8,7,];              values = A88[15:16,19:21];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[8,7,] = prob[[2]];         A88[15:16,19:21] = prob[[1]];
    }
  }
  
  else if (startP==8 && startC==8)   # 64th of a 1x1 block
  {
    # when the adjust was previously done
    if (is.na(F88[8,8,1])==FALSE)
    {
      poly = F88[8,8,];              values = A88[15:16,22:24];   
    }
    else    # the adjust is not done
    {
      prob = Adjust1ParentCD(d,dataMin,dataMax,r,0,n,maxDegree,verbose,plots);
      poly = prob[[2]];              values = prob[[1]];
      F88[8,8,] = prob[[2]];         A88[15:16,22:24] = prob[[1]];
    }
  }
  
  return(list(poly,values,F88,A88));
}
############################################################














