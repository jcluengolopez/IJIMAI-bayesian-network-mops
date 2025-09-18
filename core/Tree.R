
# CARGAR DEPENDENCIAS AL PRINCIPIO (polinomios y kernel)
library(igraph)

CreateTree = function(mutualInfoClass,verbose=FALSE)
{
  # create the maximum spanning tree using the mutual information of each variable
  # conditioned to the class
  # INPUT: - mutualInfoClass = matrix with mutual information conditioned to the class
  # OUTPUT: - links = matrix containing 1 if i is parent of j and 0 the others
  
  
  # number of variables (not including the class)
  d = dim(mutualInfoClass);
  
  # create the matrix and initialize the relations to 0
  links = matrix(0,nrow=d[1],ncol=d[1]);
  
  # order the mutual info by decreasing and removing zeroes (not necessary)
  order = sort.int(t(mutualInfoClass),decreasing=TRUE,index.return=TRUE)$ix;
  
  # select only the values with mutual information >0
  n = (d[1]^2 - d[1]) / 2;   # number of elements under the diagonal
  order = order[1:n];
  
  # get the position in the matrix of each mutual information
  nRow = ceiling(order/d[1]);
  nCol = order - (nRow-1) * d[1];
  
  # check if there are only two variables (one link)
  if (length(order)==1)
  {
    # only the first link
    links[2,1] = 1;
  }
  else
  {
    # the col that contains the row of the root must sum 0
    colRoot = nRow[1];
    
    # the first link is created between the root and the other from the first MIC
    links[nRow[1],nCol[1]] = 1;
    
    # for each pair of variables
    for (i in 2:length(order))
    {
      # sum the columns of the matrix
      sumCol = apply(links,2,sum);
      
      # get the two options with a new link
      links1 = links;        links2 = links;
      links1[nRow[i],nCol[i]] = 1;          links2[nCol[i],nRow[i]] = 1;
      
      # check if any of them is a DAG
      DAG1 = igraph :: is_dag(graph_from_adjacency_matrix(links1,mode='directed'));
      DAG2 = igraph :: is_dag(graph_from_adjacency_matrix(links2,mode='directed'));
      
      # check if the first node can be the parent
      # the second position can't be the parent neither has other parent
      if (nCol[i]!=colRoot & sumCol[nCol[i]]==0 & DAG1==TRUE)
      {
        # add the link
        links = links1;
      }
      
      # the first position can't be the parent neither has other parent
      else if (nRow[i]!=colRoot & sumCol[nRow[i]]==0 & DAG2==TRUE)
      {
        # add the link
        links = links2;
      }
    }
  }
  
  # check if the graph is a DAG
  DAG = igraph :: is_dag(graph_from_adjacency_matrix(links,mode='directed'));
  if (DAG == FALSE)   {   stop('the structure of the TAN is not a DAG');   }
  
  # add one extra row and one column for the class (the last row = 1)
  links = rbind(links,rep(1,d[1]));
  links = cbind(links,rep(0,d[1]+1));
  
  return(links);
}
############################################################

















