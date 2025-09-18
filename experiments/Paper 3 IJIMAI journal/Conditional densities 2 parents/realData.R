
# PRUEBA SOBRE APRENDIZAJEDE VARIABLES CONTINUAS CON PADRES CONTINUOS

library('squash');
library('plotly');
library('plot3D');


# ponemos el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# cambiamos el directorio para cargar las funciones
setwd("../../../core")
source("Adjust1Parent.R");          source("Adjust2Parents.R");      source("BIC.R");
source("PolynomialFit.R");          source("ShowModel.R");           source("Prediction.R");
source("AuxiliarFunctions.R");      source("Adjust2ParentsCases.R");


# ponemos el directorio de este script y lo cambiamos para cargar los datos
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../../../data/UCI")
load('appendicitis.RData');     load('hungarian.RData');
load('cloud.RData');            load('machine.RData');    load('strikes.RData');
load('adult.RData');            load('concrete.RData');   load('forestfires.RData');
load('housing.RData');          load('istanbul.RData');

# Volvemos al directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("DivergenceKL.R");




########################################################################

# DATASTRIKES

########################################################################


# Variables V4, V6 y V7 de dataStrikes

varTypeSt = DiscOrCont(dataStrikes,1);

dataMinSt = as.numeric(apply(dataStrikes,2,min));
dataMaxSt = as.numeric(apply(dataStrikes,2,max));



# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrain = sample(1:nrow(dataStrikes),floor(nrow(dataStrikes)*0.8));
Ma = c(which.max(dataStrikes[,4]),which.max(dataStrikes[,6]),which.max(dataStrikes[,7]))
mi = c(which.min(dataStrikes[,4]),which.min(dataStrikes[,6]),which.min(dataStrikes[,7]))
dataTrain = dataStrikes[c(selectTrain,Ma,mi),];
dataTest = dataStrikes[-selectTrain,];


# ajuste de variable continua con padre continuo (V6 y V7 son padres de V4)
fittMoP = Adjust2ParentsCCC(dataTrain[,4],dataMinSt[4],dataMaxSt[4],dataTrain[,6],
                     dataMinSt[6],dataMaxSt[6],dataTrain[,7],dataMinSt[7],dataMaxSt[7],n=4,7)


# DISCRETIZACIÓN

# discretización con equal width (5 intervalos)
a3 = DiscretizationEqWid(dataTrain[,c(4,6,7)],3);
dataDisc3 = a3[[1]]
breaks3 = a3[[2]]

# discretización con equal width (5 intervalos)
a5 = DiscretizationEqWid(dataTrain[,c(4,6,7)],5);
dataDisc5 = a5[[1]]
breaks5 = a5[[2]]

# ajuste discreto (en las filas de la matriz está el hijo)
fitDisc3 = Adjust2ParentsDDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3,dataDisc3[,3],1:3)
fitDisc5 = Adjust2ParentsDDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5,dataDisc5[,3],1:5)

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1];
area5 = diff(breaks5[1,])[1];
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;




# BIC

# BIC de tMoP
bictMoP = BIC2ParentsCC(fittMoP[[2]],dataTest[,4],fittMoP[[1]],dataTest[,6],dataTest[,7],
                        fittMoP[[3]],varType=1);


# BIC discretización con 3 intervalos
logDisc3 = 0;

for (i in 2:ncol(breaks3))
{
  for (j in 2:ncol(breaks3))
  {
    for (k in 2:ncol(breaks3))
    {
      aux = (dataTest[,4] >= breaks3[1,i-1] & dataTest[,4] < breaks3[1,i] &
             dataTest[,6] >= breaks3[2,j-1] & dataTest[,6] < breaks3[2,j] &
             dataTest[,7] >= breaks3[3,k-1] & dataTest[,7] < breaks3[3,k]) * 
             log10(prob3[i-1,j-1,k-1]);
      
      logDisc3 = logDisc3 + sum(aux);
    }
  }
}
bicDisc3 = logDisc3 - (ncol(breaks3)^3 - 1) / 2*log10(length(dataTest));



# BIC discretización con 5 intervalos
logDisc5 = 0;

for (i in 2:ncol(breaks5))
{
  for (j in 2:ncol(breaks5))
  {
    for (k in 2:ncol(breaks5))
    {
      aux = (dataTest[,4] >= breaks5[1,i-1] & dataTest[,4] < breaks5[1,i] &
             dataTest[,6] >= breaks5[2,j-1] & dataTest[,6] < breaks5[2,j] &
             dataTest[,7] >= breaks5[3,k-1] & dataTest[,7] < breaks5[3,k]) * 
             log10(prob5[i-1,j-1,k-1]);
      
      logDisc5 = logDisc5 + sum(aux);
    }
  }
}
bicDisc5 = logDisc5 - (ncol(breaks5)^3 - 1) / 2*log10(length(dataTest));



bictMoP
bicDisc3
bicDisc5









########################################################################


# Variables V6, V5 y V7 de dataStrikes

varTypeSt = DiscOrCont(dataStrikes,1);

dataMinSt = as.numeric(apply(dataStrikes,2,min));
dataMaxSt = as.numeric(apply(dataStrikes,2,max));



# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrain = sample(1:nrow(dataStrikes),floor(nrow(dataStrikes)*0.8));
Ma = c(which.max(dataStrikes[,6]),which.max(dataStrikes[,5]),which.max(dataStrikes[,7]))
mi = c(which.min(dataStrikes[,6]),which.min(dataStrikes[,5]),which.min(dataStrikes[,7]))
dataTrain = dataStrikes[c(selectTrain,Ma,mi),];
dataTest = dataStrikes[-selectTrain,];


# ajuste de variable continua con padre continuo (V7 y V5 son padres de V6)
fittMoP = Adjust2ParentsCCC(dataTrain[,6],dataMinSt[6],dataMaxSt[6],dataTrain[,5],
                     dataMinSt[5],dataMaxSt[5],dataTrain[,7],dataMinSt[7],dataMaxSt[7],n=4,7)


# DISCRETIZACIÓN

# discretización con equal width (5 intervalos)
a3 = DiscretizationEqWid(dataTrain[,c(6,5,7)],3);
dataDisc3 = a3[[1]]
breaks3 = a3[[2]]

# discretización con equal width (5 intervalos)
a5 = DiscretizationEqWid(dataTrain[,c(6,5,7)],5);
dataDisc5 = a5[[1]]
breaks5 = a5[[2]]

# ajuste discreto (en las filas de la matriz está el hijo)
fitDisc3 = Adjust2ParentsDDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3,dataDisc3[,3],1:3)
fitDisc5 = Adjust2ParentsDDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5,dataDisc5[,3],1:5)

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1];
area5 = diff(breaks5[1,])[1];
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;




# BIC

# BIC de tMoP
bictMoP = BIC2ParentsCC(fittMoP[[2]],dataTest[,6],fittMoP[[1]],dataTest[,5],dataTest[,7],
                        fittMoP[[3]],varType=1);


# BIC discretización con 3 intervalos
logDisc3 = 0;

for (i in 2:ncol(breaks3))
{
  for (j in 2:ncol(breaks3))
  {
    for (k in 2:ncol(breaks3))
    {
      aux = (dataTest[,6] >= breaks3[1,i-1] & dataTest[,6] < breaks3[1,i] &
             dataTest[,5] >= breaks3[2,j-1] & dataTest[,5] < breaks3[2,j] &
             dataTest[,7] >= breaks3[3,k-1] & dataTest[,7] < breaks3[3,k]) * 
             log10(prob3[i-1,j-1,k-1]);
      
      logDisc3 = logDisc3 + sum(aux);
    }
  }
}
bicDisc3 = logDisc3 - (ncol(breaks3)^3 - 1) / 2*log10(length(dataTest));



# BIC discretización con 5 intervalos
logDisc5 = 0;

for (i in 2:ncol(breaks5))
{
  for (j in 2:ncol(breaks5))
  {
    for (k in 2:ncol(breaks5))
    {
      aux = (dataTest[,6] >= breaks5[1,i-1] & dataTest[,6] < breaks5[1,i] &
             dataTest[,5] >= breaks5[2,j-1] & dataTest[,5] < breaks5[2,j] &
             dataTest[,7] >= breaks5[3,k-1] & dataTest[,7] < breaks5[3,k]) * 
             log10(prob5[i-1,j-1,k-1]);
      
      logDisc5 = logDisc5 + sum(aux);
    }
  }
}
bicDisc5 = logDisc5 - (ncol(breaks5)^3 - 1) / 2*log10(length(dataTest));



bictMoP
bicDisc3
bicDisc5






########################################################################

# DATAFORESTFIRES

########################################################################


# Variables V6, V7 y V8 de dataForestFires

varTypeFo = DiscOrCont(dataForestfires,1);

dataMinFo = as.numeric(apply(dataForestfires,2,min));
dataMaxFo = as.numeric(apply(dataForestfires,2,max));


# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrain = sample(1:nrow(dataForestfires),floor(nrow(dataForestfires)*0.8));
Ma = c(which.max(dataForestfires[,6]),which.max(dataForestfires[,7]),
       which.max(dataForestfires[,8]))
mi = c(which.min(dataForestfires[,6]),which.min(dataForestfires[,7]),
       which.min(dataForestfires[,8]))
dataTrain = dataForestfires[c(selectTrain,Ma,mi),];
dataTest = dataForestfires[-selectTrain,];


# ajuste de variable continua con padre continuo (V7 y V8 son padres de V6)
fittMoP = Adjust2ParentsCCC(dataTrain[,6],dataMinFo[6],dataMaxFo[6],dataTrain[,7],
                  dataMinFo[7],dataMaxFo[7],dataTrain[,8],dataMinFo[8],dataMaxFo[8],n=6,7)


# DISCRETIZACIÓN

# discretización con equal width (5 intervalos)
a3 = DiscretizationEqWid(dataTrain[,c(6,7,8)],3);
dataDisc3 = a3[[1]]
breaks3 = a3[[2]]

# discretización con equal width (5 intervalos)
a5 = DiscretizationEqWid(dataTrain[,c(6,7,8)],5);
dataDisc5 = a5[[1]]
breaks5 = a5[[2]]

# ajuste discreto (en las filas de la matriz está el hijo)
fitDisc3 = Adjust2ParentsDDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3,dataDisc3[,3],1:3)
fitDisc5 = Adjust2ParentsDDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5,dataDisc5[,3],1:5)

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1];
area5 = diff(breaks5[1,])[1];
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;



# BIC

# BIC de tMoP
bictMoP = BIC2ParentsCC(fittMoP[[2]],dataTest[,6],fittMoP[[1]],dataTest[,7],dataTest[,8],
                        fittMoP[[3]],varType=1);


# BIC discretización con 3 intervalos
logDisc3 = 0;

for (i in 2:ncol(breaks3))
{
  for (j in 2:ncol(breaks3))
  {
    for (k in 2:ncol(breaks3))
    {
      aux = (dataTest[,6] >= breaks3[1,i-1] & dataTest[,6] < breaks3[1,i] &
             dataTest[,7] >= breaks3[2,j-1] & dataTest[,7] < breaks3[2,j] &
             dataTest[,8] >= breaks3[3,k-1] & dataTest[,8] < breaks3[3,k]) * 
             log10(prob3[i-1,j-1,k-1]);
      
      logDisc3 = logDisc3 + sum(aux);
    }
  }
}
bicDisc3 = logDisc3 - (ncol(breaks3)^3 - 1) / 2*log10(length(dataTest));



# BIC discretización con 5 intervalos
logDisc5 = 0;

for (i in 2:ncol(breaks5))
{
  for (j in 2:ncol(breaks5))
  {
    for (k in 2:ncol(breaks5))
    {
      aux = (dataTest[,6] >= breaks5[1,i-1] & dataTest[,6] < breaks5[1,i] &
             dataTest[,7] >= breaks5[2,j-1] & dataTest[,7] < breaks5[2,j] &
             dataTest[,8] >= breaks5[3,k-1] & dataTest[,8] < breaks5[3,k]) * 
             log10(prob5[i-1,j-1,k-1]);
      
      logDisc5 = logDisc5 + sum(aux);
    }
  }
}
bicDisc5 = logDisc5 - (ncol(breaks5)^3 - 1) / 2*log10(length(dataTest));


bictMoP
bicDisc3
bicDisc5





########################################################################


# Variables V9, V10 y V13 de dataForestFires

varTypeFo = DiscOrCont(dataForestfires,1);

dataMinFo = as.numeric(apply(dataForestfires,2,min));
dataMaxFo = as.numeric(apply(dataForestfires,2,max));


# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrain = sample(1:nrow(dataForestfires),floor(nrow(dataForestfires)*0.8));
Ma = c(which.max(dataForestfires[,9]),which.max(dataForestfires[,10]),
       which.max(dataForestfires[,13]))
mi = c(which.min(dataForestfires[,9]),which.min(dataForestfires[,10]),
       which.min(dataForestfires[,13]))
dataTrain = dataForestfires[c(selectTrain,Ma,mi),];
dataTest = dataForestfires[-selectTrain,];


# ajuste de variable continua con padre continuo (V10 y V13 son padres de V9)
fittMoP = Adjust2ParentsCCC(dataTrain[,9],dataMinFo[9],dataMaxFo[9],dataTrain[,10],
                            dataMinFo[10],dataMaxFo[10],dataTrain[,13],dataMinFo[13],
                            dataMaxFo[13],n=9,7)


# DISCRETIZACIÓN

# discretización con equal width (5 intervalos)
a3 = DiscretizationEqWid(dataTrain[,c(9,10,13)],3);
dataDisc3 = a3[[1]]
breaks3 = a3[[2]]

# discretización con equal width (5 intervalos)
a5 = DiscretizationEqWid(dataTrain[,c(9,10,13)],5);
dataDisc5 = a5[[1]]
breaks5 = a5[[2]]

# ajuste discreto (en las filas de la matriz está el hijo)
fitDisc3 = Adjust2ParentsDDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3,dataDisc3[,3],1:3)
fitDisc5 = Adjust2ParentsDDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5,dataDisc5[,3],1:5)

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1];
area5 = diff(breaks5[1,])[1];
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;



# BIC

# BIC de tMoP
bictMoP = BIC2ParentsCC(fittMoP[[2]],dataTest[,9],fittMoP[[1]],dataTest[,10],dataTest[,13],
                        fittMoP[[3]],varType=1);


# BIC discretización con 3 intervalos
logDisc3 = 0;

for (i in 2:ncol(breaks3))
{
  for (j in 2:ncol(breaks3))
  {
    for (k in 2:ncol(breaks3))
    {
      aux = (dataTest[,9]  >= breaks3[1,i-1] & dataTest[,9]  < breaks3[1,i] &
             dataTest[,10] >= breaks3[2,j-1] & dataTest[,10] < breaks3[2,j] &
             dataTest[,13] >= breaks3[3,k-1] & dataTest[,13] < breaks3[3,k]) * 
             log10(prob3[i-1,j-1,k-1]);
      
      logDisc3 = logDisc3 + sum(aux);
    }
  }
}
bicDisc3 = logDisc3 - (ncol(breaks3)^3 - 1) / 2*log10(length(dataTest));



# BIC discretización con 5 intervalos
logDisc5 = 0;

for (i in 2:ncol(breaks5))
{
  for (j in 2:ncol(breaks5))
  {
    for (k in 2:ncol(breaks5))
    {
      aux = (dataTest[,9]  >= breaks5[1,i-1] & dataTest[,9]  < breaks5[1,i] &
             dataTest[,10] >= breaks5[2,j-1] & dataTest[,10] < breaks5[2,j] &
             dataTest[,13] >= breaks5[3,k-1] & dataTest[,13] < breaks5[3,k]) * 
             log10(prob5[i-1,j-1,k-1]);
      
      logDisc5 = logDisc5 + sum(aux);
    }
  }
}
bicDisc5 = logDisc5 - (ncol(breaks5)^3 - 1) / 2*log10(length(dataTest));


bictMoP
bicDisc3
bicDisc5





########################################################################

# DATAISTANBUL

########################################################################


# Variables V1, V2 y V3 de dataIstanbul

varTypeIs = DiscOrCont(dataIstanbul,1);

dataMinIs = as.numeric(apply(dataIstanbul,2,min));
dataMaxIs = as.numeric(apply(dataIstanbul,2,max));


# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrain = sample(1:nrow(dataIstanbul),floor(nrow(dataIstanbul)*0.8));
Ma = c(which.max(dataIstanbul[,1]),which.max(dataIstanbul[,2]),which.max(dataIstanbul[,3]))
mi = c(which.min(dataIstanbul[,1]),which.min(dataIstanbul[,2]),which.min(dataIstanbul[,3]))
dataTrain = dataIstanbul[c(selectTrain,Ma,mi),];
dataTest = dataIstanbul[-selectTrain,];


# ajuste de variable continua con padre continuo (V2 y V3 son padres de V1)
fittMoP = Adjust2ParentsCCC(dataTrain[,1],dataMinIs[1],dataMaxIs[1],dataTrain[,2],
                    dataMinIs[2],dataMaxIs[2],dataTrain[,3],dataMinIs[3],dataMaxIs[3],n=1,7)


# DISCRETIZACIÓN

# discretización con equal width (5 intervalos)
a3 = DiscretizationEqWid(dataTrain[,c(1,2,3)],3);
dataDisc3 = a3[[1]]
breaks3 = a3[[2]]

# discretización con equal width (5 intervalos)
a5 = DiscretizationEqWid(dataTrain[,c(1,2,3)],5);
dataDisc5 = a5[[1]]
breaks5 = a5[[2]]

# ajuste discreto (en las filas de la matriz está el hijo)
fitDisc3 = Adjust2ParentsDDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3,dataDisc3[,3],1:3)
fitDisc5 = Adjust2ParentsDDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5,dataDisc5[,3],1:5)

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1];
area5 = diff(breaks5[1,])[1];
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;



# BIC

# BIC de tMoP
bictMoP = BIC2ParentsCC(fittMoP[[2]],dataTest[,1],fittMoP[[1]],dataTest[,2],dataTest[,3],
                        fittMoP[[3]],varType=1);


# BIC discretización con 3 intervalos
logDisc3 = 0;

for (i in 2:ncol(breaks3))
{
  for (j in 2:ncol(breaks3))
  {
    for (k in 2:ncol(breaks3))
    {
      aux = (dataTest[,1] >= breaks3[1,i-1] & dataTest[,1] < breaks3[1,i] &
             dataTest[,2] >= breaks3[2,j-1] & dataTest[,2] < breaks3[2,j] &
             dataTest[,3] >= breaks3[3,k-1] & dataTest[,3] < breaks3[3,k]) * 
             log10(prob3[i-1,j-1,k-1]);
      
      logDisc3 = logDisc3 + sum(aux);
    }
  }
}
bicDisc3 = logDisc3 - (ncol(breaks3)^3 - 1) / 2*log10(length(dataTest));



# BIC discretización con 5 intervalos
logDisc5 = 0;

for (i in 2:ncol(breaks5))
{
  for (j in 2:ncol(breaks5))
  {
    for (k in 2:ncol(breaks5))
    {
      aux = (dataTest[,1] >= breaks5[1,i-1] & dataTest[,1] < breaks5[1,i] &
             dataTest[,2] >= breaks5[2,j-1] & dataTest[,2] < breaks5[2,j] &
             dataTest[,3] >= breaks5[3,k-1] & dataTest[,3] < breaks5[3,k]) * 
             log10(prob5[i-1,j-1,k-1]);
      
      logDisc5 = logDisc5 + sum(aux);
    }
  }
}
bicDisc5 = logDisc5 - (ncol(breaks5)^3 - 1) / 2*log10(length(dataTest));


bictMoP
bicDisc3
bicDisc5






########################################################################


# Variables V4, V5 y V6 de dataIstanbul

varTypeIs = DiscOrCont(dataIstanbul,1);

dataMinIs = as.numeric(apply(dataIstanbul,2,min));
dataMaxIs = as.numeric(apply(dataIstanbul,2,max));


# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrain = sample(1:nrow(dataIstanbul),floor(nrow(dataIstanbul)*0.8));
Ma = c(which.max(dataIstanbul[,4]),which.max(dataIstanbul[,5]),which.max(dataIstanbul[,6]))
mi = c(which.min(dataIstanbul[,4]),which.min(dataIstanbul[,5]),which.min(dataIstanbul[,6]))
dataTrain = dataIstanbul[c(selectTrain,Ma,mi),];
dataTest = dataIstanbul[-selectTrain,];


# ajuste de variable continua con padre continuo (V5 y V6 son padres de V4)
fittMoP = Adjust2ParentsCCC(dataTrain[,4],dataMinIs[4],dataMaxIs[4],dataTrain[,5],
                    dataMinIs[5],dataMaxIs[5],dataTrain[,6],dataMinIs[6],dataMaxIs[6],n=4,7)


# DISCRETIZACIÓN

# discretización con equal width (5 intervalos)
a3 = DiscretizationEqWid(dataTrain[,c(4,5,6)],3);
dataDisc3 = a3[[1]]
breaks3 = a3[[2]]

# discretización con equal width (5 intervalos)
a5 = DiscretizationEqWid(dataTrain[,c(4,5,6)],5);
dataDisc5 = a5[[1]]
breaks5 = a5[[2]]

# ajuste discreto (en las filas de la matriz está el hijo)
fitDisc3 = Adjust2ParentsDDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3,dataDisc3[,3],1:3)
fitDisc5 = Adjust2ParentsDDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5,dataDisc5[,3],1:5)

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1];
area5 = diff(breaks5[1,])[1];
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;



# BIC

# BIC de tMoP
bictMoP = BIC2ParentsCC(fittMoP[[2]],dataTest[,4],fittMoP[[1]],dataTest[,5],dataTest[,6],
                        fittMoP[[3]],varType=1);


# BIC discretización con 3 intervalos
logDisc3 = 0;

for (i in 2:ncol(breaks3))
{
  for (j in 2:ncol(breaks3))
  {
    for (k in 2:ncol(breaks3))
    {
      aux = (dataTest[,4] >= breaks3[1,i-1] & dataTest[,4] < breaks3[1,i] &
             dataTest[,5] >= breaks3[2,j-1] & dataTest[,5] < breaks3[2,j] &
             dataTest[,6] >= breaks3[3,k-1] & dataTest[,6] < breaks3[3,k]) * 
             log10(prob3[i-1,j-1,k-1]);
      
      logDisc3 = logDisc3 + sum(aux);
    }
  }
}
bicDisc3 = logDisc3 - (ncol(breaks3)^3 - 1) / 2*log10(length(dataTest));



# BIC discretización con 5 intervalos
logDisc5 = 0;

for (i in 2:ncol(breaks5))
{
  for (j in 2:ncol(breaks5))
  {
    for (k in 2:ncol(breaks5))
    {
      aux = (dataTest[,4] >= breaks5[1,i-1] & dataTest[,4] < breaks5[1,i] &
             dataTest[,5] >= breaks5[2,j-1] & dataTest[,5] < breaks5[2,j] &
             dataTest[,6] >= breaks5[3,k-1] & dataTest[,6] < breaks5[3,k]) * 
             log10(prob5[i-1,j-1,k-1]);
      
      logDisc5 = logDisc5 + sum(aux);
    }
  }
}
bicDisc5 = logDisc5 - (ncol(breaks5)^3 - 1) / 2*log10(length(dataTest));


bictMoP
bicDisc3
bicDisc5






































