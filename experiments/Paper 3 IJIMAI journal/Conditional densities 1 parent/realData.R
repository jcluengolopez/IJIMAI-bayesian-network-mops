
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



########################################################################

# DATASTRIKES

########################################################################


# Variables V5 y V6 de dataStrikes

varTypeSt = DiscOrCont(dataStrikes,1);

dataMinSt = as.numeric(apply(dataStrikes,2,min));
dataMaxSt = as.numeric(apply(dataStrikes,2,max));


# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrain = sample(1:nrow(dataStrikes),floor(nrow(dataStrikes)*0.8));
Ma = c(which.max(dataStrikes[,5]),which.max(dataStrikes[,6]))
mi = c(which.min(dataStrikes[,5]),which.min(dataStrikes[,6]))
dataTrain = dataStrikes[c(selectTrain,Ma,mi),];
dataTest = dataStrikes[-selectTrain,];


# ajuste de variable continua con padre continuo (V6 es padre de V5)
fittMoP = Adjust1ParentCC(dataTrain[,5],dataMinSt[5],dataMaxSt[5],dataTrain[,6],
                          dataMinSt[6],dataMaxSt[6],n=5,maxDegree=7);



# DISCRETIZACIÓN

# discretización con equal width (5 intervalos)
a3 = DiscretizationEqWid(dataTrain[,c(5,6)],3);
dataDisc3 = a3[[1]]
breaks3 = a3[[2]]

# discretización con equal width (5 intervalos)
a5 = DiscretizationEqWid(dataTrain[,c(5,6)],5);
dataDisc5 = a5[[1]]
breaks5 = a5[[2]]

# ajuste discreto (en las filas de la matriz está el hijo)
fitDisc3 = Adjust1ParentDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3);
fitDisc5 = Adjust1ParentDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5);

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1];
area5 = diff(breaks5[1,])[1];
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;



# BIC

# BIC de tMoP
bictMoP = BIC1ParentC(fittMoP[[2]],dataTest[,5],fittMoP[[1]],dataTest[,6],fittMoP[[3]],
                      varType=1);


# BIC discretización con 3 intervalos
logDisc3 = 0;

for (i in 2:ncol(breaks3))
{
  for (j in 2:ncol(breaks3))
  {
    aux = (dataTest[,5] >= breaks3[1,i-1] & dataTest[,5] < breaks3[1,i] &
        dataTest[,6] >= breaks3[2,j-1] & dataTest[,6] < breaks3[2,j]) * log10(prob3[i-1,j-1])

    logDisc3 = logDisc3 + sum(aux);
  }
}

bicDisc3 = logDisc3 - (ncol(breaks3)^2 - 1) / 2*log10(length(dataTest));




# BIC discretización con 5 intervalos
logDisc5 = 0;

for (i in 2:ncol(breaks5))
{
  for (j in 2:ncol(breaks5))
  {
    aux = (dataTest[,5] >= breaks5[1,i-1] & dataTest[,5] < breaks5[1,i] &
        dataTest[,6] >= breaks5[2,j-1] & dataTest[,6] < breaks5[2,j]) * log10(prob5[i-1,j-1])

    logDisc5 = logDisc5 + sum(aux);
  }
}

bicDisc5 = logDisc5 - (ncol(breaks5)^2 - 1) / 2*log10(length(dataTest));

bictMoP
bicDisc3
bicDisc5




# GRÁFICAS

# vectores para dibujar las gráficas
x1 = seq(dataMinSt[5],dataMaxSt[5],length.out=100);
y1 = seq(dataMinSt[6],dataMaxSt[6],length.out=100);


# gráfica original
x = seq(dataMinSt[5],dataMaxSt[5],length.out=20);   # hijo
y = seq(dataMinSt[6],dataMaxSt[6],length.out=20);   # padre
# la función Histograma3D está al final de script
z1 = Histograma3D(dataStrikes[,5],dataStrikes[,6],x,y);

plot_ly(x=x[-1],y=y[-1],z=t(z1),type='surface');


# gráfica del ajuste con tMoP
z2 = matrix(,nrow=length(x1),ncol=length(y1));

for (i in 1:length(x1))
{
  z2[i,] = PredicttMoPCC(x1,y1[i],fittMoP[[2]],fittMoP[[1]],fittMoP[[3]]);
}

plot_ly(x=x1,y=y1,z=z2,type='surface');



# gráfica del ajuste discretizado con 3 intervalos
z3 = matrix(,nrow=length(x1),ncol = length(x1));
for (i in 1:length(x1))
{
  for (j in 1:length(y1))
  {
    row = which(x1[i] < breaks3[1,])[1] - 1;
    column = which(y1[j] < breaks3[2,])[1] - 1;
    if (is.na(row))   row=3
    else if (row==0)   row=1
    if (is.na(column))   column=3
    else if (column==0)   column=1

    z3[i,j] = prob5[row,column]
  }
}

plot_ly(x=x1,y=y1,z=z3,type='surface');


# gráfica del ajuste discretizado con 5 intervalos
z4 = matrix(,nrow=length(x1),ncol = length(x1));
for (i in 1:length(x1))
{
  for (j in 1:length(y1))
  {
    row = which(x1[i] < breaks5[1,])[1] - 1;
    column = which(y1[j] < breaks5[2,])[1] - 1;
    if (is.na(row))   row=5
    else if (row==0)   row=1
    if (is.na(column))   column=5
    else if (column==0)   column=1

    z4[i,j] = prob5[row,column]
  }
}

plot_ly(x=x1,y=y1,z=z4,type='surface');





########################################################################


# Variables V4 y V7 de dataStrikes

varTypeSt = DiscOrCont(dataStrikes,1);

dataMinSt = as.numeric(apply(dataStrikes,2,min));
dataMaxSt = as.numeric(apply(dataStrikes,2,max));


# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrain = sample(1:nrow(dataStrikes),floor(nrow(dataStrikes)*0.8));
Ma = c(which.max(dataStrikes[,4]),which.max(dataStrikes[,7]))
mi = c(which.min(dataStrikes[,4]),which.min(dataStrikes[,7]))
dataTrain = dataStrikes[c(selectTrain,Ma,mi),];
dataTrain = dataStrikes[selectTrain,];
dataTest = dataStrikes[-selectTrain,];


# ajuste de variable continua con padre continuo (V7 es padre de V4)
fittMoP = Adjust1ParentCC(dataTrain[,4],dataMinSt[4],dataMaxSt[4],dataTrain[,7],
                          dataMinSt[7],dataMaxSt[7],n=5,maxDegree=7);



# DISCRETIZACIÓN

# discretización con equal width (5 intervalos)
a3 = DiscretizationEqWid(dataTrain[,c(4,7)],3);
dataDisc3 = a3[[1]]
breaks3 = a3[[2]]
a5 = DiscretizationEqWid(dataTrain[,c(4,7)],5);
dataDisc5 = a5[[1]]
breaks5 = a5[[2]]

# ajuste discreto (en las filas de la matriz está el hijo)
fitDisc3 = Adjust1ParentDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3);
fitDisc5 = Adjust1ParentDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5);

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1];
area5 = diff(breaks5[1,])[1];
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;



# BIC

# BIC de tMoP
bictMoP = BIC1ParentC(fittMoP[[2]],dataTest[,4],fittMoP[[1]],dataTest[,7],fittMoP[[3]],
                      varType=1);


# BIC discretización con 3 invertalos
logDisc3 = 0;

for (i in 2:ncol(breaks3))
{
  for (j in 2:ncol(breaks3))
  {
    aux = (dataTest[,4] >= breaks3[1,i-1] & dataTest[,4] < breaks3[1,i] &
        dataTest[,7] >= breaks3[2,j-1] & dataTest[,7] < breaks3[2,j]) * log10(prob3[i-1,j-1])

    logDisc3 = logDisc3 + sum(aux);
  }
}
bicDisc3 = logDisc3 - (ncol(breaks3)^2 - 1) / 2*log10(length(dataTest));


# BIC discretización con 5 invertalos
logDisc5 = 0;

for (i in 2:ncol(breaks5))
{
  for (j in 2:ncol(breaks5))
  {
    aux = (dataTest[,4] >= breaks5[1,i-1] & dataTest[,4] < breaks5[1,i] &
        dataTest[,7] >= breaks5[2,j-1] & dataTest[,7] < breaks5[2,j]) * log10(prob5[i-1,j-1])

    logDisc5 = logDisc5 + sum(aux);
  }
}
bicDisc5 = logDisc5 - (ncol(breaks5)^2 - 1) / 2*log10(length(dataTest));

bictMoP
bicDisc3
bicDisc5





# GRÁFICAS

# vectores para dibujar las gráficas
x1 = seq(dataMinSt[4],dataMaxSt[4],length.out=100);
y1 = seq(dataMinSt[7],dataMaxSt[7],length.out=100);


# gráfica original
x = seq(dataMinSt[4],dataMaxSt[4],length.out=20);   # hijo
y = seq(dataMinSt[7],dataMaxSt[7],length.out=20);   # padre
z1 = Histograma3D(dataStrikes[,4],dataStrikes[,7],x,y);

plot_ly(x=x[-1],y=y[-1],z=t(z1),type='surface');


# gráfica del ajuste con tMoP
z2 = matrix(,nrow=length(x1),ncol=length(y1));

for (i in 1:length(x1))
{
  z2[i,] = PredicttMoPCC(x1,y1[i],fittMoP[[2]],fittMoP[[1]],fittMoP[[3]]);
}

plot_ly(x=x1,y=y1,z=z2,type='surface');



# gráfica del ajuste discretizado con 3 intervalos
z3 = matrix(,nrow=length(x1),ncol = length(x1));
for (i in 1:length(x1))
{
  for (j in 1:length(y1))
  {
    row = which(x1[i] < breaks3[1,])[1] - 1;
    column = which(y1[j] < breaks3[2,])[1] - 1;
    if (is.na(row))   row=3
    else if (row==0)   row=1
    if (is.na(column))   column=3
    else if (column==0)   column=1

    z3[i,j] = prob3[row,column]
  }
}

plot_ly(x=x1,y=y1,z=t(z3),type='surface');




# gráfica del ajuste discretizado con 5 intervalos
z4 = matrix(,nrow=length(x1),ncol = length(x1));
for (i in 1:length(x1))
{
  for (j in 1:length(y1))
  {
    row = which(x1[i] < breaks5[1,])[1] - 1;
    column = which(y1[j] < breaks5[2,])[1] - 1;
    if (is.na(row))   row=5
    else if (row==0)   row=1
    if (is.na(column))   column=5
    else if (column==0)   column=1

    z4[i,j] = prob5[row,column]
  }
}

plot_ly(x=x1,y=y1,z=t(z4),type='surface');






########################################################################

# DATAFORESTFIRES

########################################################################


# Variables V6 y V9 de dataForestFires

varTypeFo = DiscOrCont(dataForestfires,1);

dataMinFo = as.numeric(apply(dataForestfires,2,min));
dataMaxFo = as.numeric(apply(dataForestfires,2,max));


# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrain = sample(1:nrow(dataForestfires),floor(nrow(dataForestfires)*0.8));
Ma = c(which.max(dataForestfires[,6]),which.max(dataForestfires[,9]))
mi = c(which.min(dataForestfires[,6]),which.min(dataForestfires[,9]))
dataTrain = dataForestfires[c(selectTrain,Ma,mi),];
dataTest = dataForestfires[-selectTrain,];



# ajuste de variable continua con padre continuo (V9 es padre de V6)
fittMoP = Adjust1ParentCC(dataTrain[,6],dataMinFo[6],dataMaxFo[6],dataTrain[,9],
                          dataMinFo[9],dataMaxFo[9],n=9,maxDegree=7);


# DISCRETIZACIÓN

# discretización con equal width (3 intervalos)
a3 = DiscretizationEqWid(dataTrain[,c(6,9)],3);
dataDisc3 = a3[[1]]
breaks3 = a3[[2]]

# discretización con equal width (5 intervalos)
a5 = DiscretizationEqWid(dataTrain[,c(6,9)],5);
dataDisc5 = a5[[1]]
breaks5 = a5[[2]]

# ajuste discreto (en las filas de la matriz está el hijo)
fitDisc3 = Adjust1ParentDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3);
fitDisc5 = Adjust1ParentDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5);

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1];
area5 = diff(breaks5[1,])[1];
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;



# BIC

# BIC de tMoP
bictMoP = BIC1ParentC(fittMoP[[2]],dataTest[,6],fittMoP[[1]],dataTest[,9],fittMoP[[3]],
                      varType=1);


# BIC discretización con 3 intervalos
logDisc3 = 0;

for (i in 2:ncol(breaks3))
{
  for (j in 2:ncol(breaks3))
  {
    aux = (dataTest[,6] >= breaks3[1,i-1] & dataTest[,6] < breaks3[1,i] &
        dataTest[,9] >= breaks3[2,j-1] & dataTest[,9] < breaks3[2,j]) * log10(prob3[i-1,j-1])

    logDisc3 = logDisc3 + sum(aux);
  }
}
bicDisc3 = logDisc3 - (ncol(breaks3)^2 - 1) / 2*log10(length(dataTest));



# BIC discretización con 5 intervalos
logDisc5 = 0;

for (i in 2:ncol(breaks5))
{
  for (j in 2:ncol(breaks5))
  {
    aux = (dataTest[,6] >= breaks5[1,i-1] & dataTest[,6] < breaks5[1,i] &
        dataTest[,9] >= breaks5[2,j-1] & dataTest[,9] < breaks5[2,j]) * log10(prob5[i-1,j-1])

    logDisc5 = logDisc5 + sum(aux);
  }
}
bicDisc5 = logDisc5 - (ncol(breaks5)^2 - 1) / 2*log10(length(dataTest));

bictMoP
bicDisc3
bicDisc5



# GRÁFICAS

# vectores para dibujar las gráficas
x1 = seq(dataMinFo[6],dataMaxFo[6],length.out=200);
y1 = seq(dataMinFo[9],dataMaxFo[9],length.out=200);


# gráfica original
x = seq(dataMinFo[6],dataMaxFo[6],length.out=25);   # hijo
y = seq(dataMinFo[9],dataMaxFo[9],length.out=25);   # padre
# la función Histograma3D está al final de script
z1 = Histograma3D(dataForestfires[,6],dataForestfires[,9],x,y);


## GUARDAR EN PDF
pdf('GraficaCondForestFires6-9-orig.pdf',width=7,height=5)

par(mar = c(1,0,1,2))
colors = colorRampPalette(c("darkblue", "darkgreen", "yellow"))
persp3D(x[-1],y[-1],t(z1),col=colors(256),theta=140,phi=35,ticktype="detailed",lwd=0.6)
title(main='Histogram of original data')

#plot_ly(x=x[-1],y=y[-1],z=t(z1),type='surface');
dev.off()




# gráfica del ajuste con tMoP
z2 = matrix(,nrow=length(x1),ncol=length(y1));

for (i in 1:length(x1))
{
  z2[i,] = PredicttMoPCC(x1,y1[i],fittMoP[[2]],fittMoP[[1]],fittMoP[[3]]);
}

## GUARDAR EN PDF
pdf('GraficaCondForestFires6-9-tmop.pdf',width=7,height=5)

par(mar = c(1,0,1,2))
colors = colorRampPalette(c("darkblue", "darkgreen", "yellow"))
persp3D(x1,y1,z2,col=colors(256),theta=140,phi=35,ticktype="detailed",lwd=0.6)
title(main='tMoP fit')

#plot_ly(x=x1,y=y1,z=z2,type='surface');
dev.off()




# gráfica del ajuste discretizado con 3 intervalos
z3 = matrix(,nrow=length(x1),ncol = length(x1));
for (i in 1:length(x1))
{
  for (j in 1:length(y1))
  {
    row = which(x1[i] < breaks3[1,])[1] - 1;
    column = which(y1[j] < breaks3[2,])[1] - 1;
    if (is.na(row))   row=3
    else if (row==0)   row=1
    if (is.na(column))   column=3
    else if (column==0)   column=1

    z3[i,j] = prob3[row,column]
  }
}

## GUARDAR EN PDF
pdf('GraficaCondForestFires6-9-disc3.pdf',width=7,height=5)

par(mar = c(1,0,1,2))
colors = colorRampPalette(c("darkblue", "darkgreen", "yellow"))
persp3D(x1,y1,t(z3),col=colors(256),theta=140,phi=35,ticktype="detailed",lwd=0.6)
title(main='Discretization in 3 intervals per variable')

#plot_ly(x=x1,y=y1,z=z4,type='surface');
dev.off()





# gráfica del ajuste discretizado con 5 intervalos
z4 = matrix(,nrow=length(x1),ncol = length(x1));
for (i in 1:length(x1))
{
  for (j in 1:length(y1))
  {
    row = which(x1[i] < breaks5[1,])[1] - 1;
    column = which(y1[j] < breaks5[2,])[1] - 1;
    if (is.na(row))   row=5
    else if (row==0)   row=1
    if (is.na(column))   column=5
    else if (column==0)   column=1

    z4[i,j] = prob5[row,column]
  }
}

## GUARDAR EN PDF
pdf('GraficaCondForestFires6-9-disc5.pdf',width=7,height=5)

par(mar = c(1,0,1,2))
colors = colorRampPalette(c("darkblue", "darkgreen", "yellow"))
persp3D(x1,y1,t(z4),col=colors(256),theta=140,phi=35,ticktype="detailed",lwd=0.6)
title(main='Discretization in 5 intervals per variable')

#plot_ly(x=x1,y=y1,z=z4,type='surface');
dev.off()



########################################################################



# Variables V5 y V13 de dataForestFires

varTypeFo = DiscOrCont(dataForestfires,1);

dataMinFo = as.numeric(apply(dataForestfires,2,min));
dataMaxFo = as.numeric(apply(dataForestfires,2,max));


# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrain = sample(1:nrow(dataForestfires),floor(nrow(dataForestfires)*0.8));
Ma = c(which.max(dataForestfires[,5]),which.max(dataForestfires[,13]))
mi = c(which.min(dataForestfires[,5]),which.min(dataForestfires[,13]))
dataTrain = dataForestfires[c(selectTrain,Ma,mi),];
dataTest = dataForestfires[-selectTrain,];



# ajuste de variable continua con padre continuo (V13 es padre de V5)
fittMoP = Adjust1ParentCC(dataTrain[,5],dataMinFo[5],dataMaxFo[5],dataTrain[,13],
                          dataMinFo[13],dataMaxFo[13],n=9,maxDegree=7);


# DISCRETIZACIÓN

# discretización con equal width (3 intervalos)
a3 = DiscretizationEqWid(dataTrain[,c(5,13)],3);
dataDisc3 = a3[[1]]
breaks3 = a3[[2]]

# discretización con equal width (5 intervalos)
a5 = DiscretizationEqWid(dataTrain[,c(5,13)],5);
dataDisc5 = a5[[1]]
breaks5 = a5[[2]]

# ajuste discreto (en las filas de la matriz está el hijo)
fitDisc3 = Adjust1ParentDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3);
fitDisc5 = Adjust1ParentDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5);

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1];
area5 = diff(breaks5[1,])[1];
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;



# BIC

# BIC de tMoP
bictMoP = BIC1ParentC(fittMoP[[2]],dataTest[,5],fittMoP[[1]],dataTest[,13],fittMoP[[3]],
                      varType=1);


# BIC discretización con 3 intervalos
logDisc3 = 0;

for (i in 2:ncol(breaks3))
{
  for (j in 2:ncol(breaks3))
  {
    aux = (dataTest[,5] >= breaks3[1,i-1] & dataTest[,5] < breaks3[1,i] &
      dataTest[,13] >= breaks3[2,j-1] & dataTest[,13] < breaks3[2,j]) * log10(prob3[i-1,j-1])

    logDisc3 = logDisc3 + sum(aux);
  }
}
bicDisc3 = logDisc3 - (ncol(breaks3)^2 - 1) / 2*log10(length(dataTest));


# BIC discretización con 5 intervalos
logDisc5 = 0;

for (i in 2:ncol(breaks5))
{
  for (j in 2:ncol(breaks5))
  {
    aux = (dataTest[,5] >= breaks5[1,i-1] & dataTest[,5] < breaks5[1,i] &
      dataTest[,13] >= breaks5[2,j-1] & dataTest[,13] < breaks5[2,j]) * log10(prob5[i-1,j-1])

    logDisc5 = logDisc5 + sum(aux);
  }
}
bicDisc5 = logDisc5 - (ncol(breaks5)^2 - 1) / 2*log10(length(dataTest));

bictMoP
bicDisc3
bicDisc5



# GRÁFICAS

# vectores para dibujar las gráficas
x1 = seq(dataMinFo[5],dataMaxFo[5],length.out=100);
y1 = seq(dataMinFo[13],dataMaxFo[13],length.out=100);


# gráfica original
x = seq(dataMinFo[5],dataMaxFo[5],length.out=20);   # hijo
y = seq(dataMinFo[13],dataMaxFo[13],length.out=20);   # padre
# la función Histograma3D está al final de script
z1 = Histograma3D(dataForestfires[,5],dataForestfires[,13],x,y);

plot_ly(x=x[-1],y=y[-1],z=t(z1),type='surface');


# gráfica del ajuste con tMoP
z2 = matrix(,nrow=length(x1),ncol=length(y1));

for (i in 1:length(x1))
{
  z2[i,] = PredicttMoPCC(x1,y1[i],fittMoP[[2]],fittMoP[[1]],fittMoP[[3]]);
}

plot_ly(x=x1,y=y1,z=z2,type='surface');


# gráfica del ajuste discretizado con 3 intervalos
z3 = matrix(,nrow=length(x1),ncol = length(x1));
for (i in 1:length(x1))
{
  for (j in 1:length(y1))
  {
    row = which(x1[i] < breaks3[1,])[1] - 1;
    column = which(y1[j] < breaks3[2,])[1] - 1;
    if (is.na(row))   row=3
    else if (row==0)   row=1
    if (is.na(column))   column=3
    else if (column==0)   column=1

    z3[i,j] = prob3[row,column]
  }
}

plot_ly(x=x1,y=y1,z=z3,type='surface');



# gráfica del ajuste discretizado con 5 intervalos
z4 = matrix(,nrow=length(x1),ncol = length(x1));
for (i in 1:length(x1))
{
  for (j in 1:length(y1))
  {
    row = which(x1[i] < breaks5[1,])[1] - 1;
    column = which(y1[j] < breaks5[2,])[1] - 1;
    if (is.na(row))   row=5
    else if (row==0)   row=1
    if (is.na(column))   column=5
    else if (column==0)   column=1

    z4[i,j] = prob5[row,column]
  }
}

plot_ly(x=x1,y=y1,z=z4,type='surface');





########################################################################

# DATAISTANBUL

########################################################################


# Variables V6 y V9 de dataIstanbul

varTypeIs = DiscOrCont(dataIstanbul,1);

dataMinIs = as.numeric(apply(dataIstanbul,2,min));
dataMaxIs = as.numeric(apply(dataIstanbul,2,max));


# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrain = sample(1:nrow(dataIstanbul),floor(nrow(dataIstanbul)*0.8));
Ma = c(which.max(dataIstanbul[,6]),which.max(dataIstanbul[,9]))
mi = c(which.min(dataIstanbul[,6]),which.min(dataIstanbul[,9]))
dataTrain = dataIstanbul[c(selectTrain,Ma,mi),];
dataTest = dataIstanbul[-selectTrain,];


# ajuste de variable continua con padre continuo (V9 es padre de V6)
fittMoP = Adjust1ParentCC(dataTrain[,6],dataMinIs[6],dataMaxIs[6],dataTrain[,9],
                          dataMinIs[9],dataMaxIs[9],n=6,maxDegree=7);


# DISCRETIZACIÓN

# discretización con equal width (3 intervalos)
a3 = DiscretizationEqWid(dataTrain[,c(6,9)],3);
dataDisc3 = a3[[1]]
breaks3 = a3[[2]]

# discretización con equal width (5 intervalos)
a5 = DiscretizationEqWid(dataTrain[,c(6,9)],5);
dataDisc5 = a5[[1]]
breaks5 = a5[[2]]

# ajuste discreto (en las filas de la matriz está el hijo)
fitDisc3 = Adjust1ParentDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3);
fitDisc5 = Adjust1ParentDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5);

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1];
area5 = diff(breaks5[1,])[1];
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;



# BIC

# BIC de tMoP
bictMoP = BIC1ParentC(fittMoP[[2]],dataTest[,6],fittMoP[[1]],dataTest[,9],fittMoP[[3]],
                      varType=1);


# BIC discretización con 3 intervalos
logDisc3 = 0;

for (i in 2:ncol(breaks3))
{
  for (j in 2:ncol(breaks3))
  {
    aux = (dataTest[,6] >= breaks3[1,i-1] & dataTest[,6] < breaks3[1,i] &
        dataTest[,9] >= breaks3[2,j-1] & dataTest[,9] < breaks3[2,j]) * log10(prob3[i-1,j-1])

    logDisc3 = logDisc3 + sum(aux);
  }
}
bicDisc3 = logDisc3 - (ncol(breaks3)^2 - 1) / 2*log10(length(dataTest));



# BIC discretización con 5 intervalos
logDisc5 = 0;

for (i in 2:ncol(breaks5))
{
  for (j in 2:ncol(breaks5))
  {
    aux = (dataTest[,6] >= breaks5[1,i-1] & dataTest[,6] < breaks5[1,i] &
        dataTest[,9] >= breaks5[2,j-1] & dataTest[,9] < breaks5[2,j]) * log10(prob5[i-1,j-1])

    logDisc5 = logDisc5 + sum(aux);
  }
}
bicDisc5 = logDisc5 - (ncol(breaks5)^2 - 1) / 2*log10(length(dataTest));

bictMoP
bicDisc3
bicDisc5



# GRÁFICAS

# vectores para dibujar las gráficas
x1 = seq(dataMinIs[6],dataMaxIs[6],length.out=100);
y1 = seq(dataMinIs[9],dataMaxIs[9],length.out=100);


# gráfica original
x = seq(dataMinIs[6],dataMaxIs[6],length.out=20);   # hijo
y = seq(dataMinIs[9],dataMaxIs[9],length.out=20);   # padre
# la función Histograma3D está al final de script
z1 = Histograma3D(dataIstanbul[,6],dataIstanbul[,9],x,y);

plot_ly(x=x[-1],y=y[-1],z=t(z1),type='surface');


# gráfica del ajuste con tMoP
z2 = matrix(,nrow=length(x1),ncol=length(y1));

for (i in 1:length(x1))
{
  z2[i,] = PredicttMoPCC(x1,y1[i],fittMoP[[2]],fittMoP[[1]],fittMoP[[3]]);
}

plot_ly(x=x1,y=y1,z=z2,type='surface');


# gráfica del ajuste discretizado con 5 intervalos
z3 = matrix(,nrow=length(x1),ncol = length(x1));
for (i in 1:length(x1))
{
  for (j in 1:length(y1))
  {
    row = which(x1[i] < breaks3[1,])[1] - 1;
    column = which(y1[j] < breaks3[2,])[1] - 1;
    if (is.na(row))   row=3
    else if (row==0)   row=1
    if (is.na(column))   column=3
    else if (column==0)   column=1

    z3[i,j] = prob3[row,column]
  }
}

plot_ly(x=x1,y=y1,z=z3,type='surface');


# gráfica del ajuste discretizado con 5 intervalos
z4 = matrix(,nrow=length(x1),ncol = length(x1));
for (i in 1:length(x1))
{
  for (j in 1:length(y1))
  {
    row = which(x1[i] < breaks5[1,])[1] - 1;
    column = which(y1[j] < breaks5[2,])[1] - 1;
    if (is.na(row))   row=5
    else if (row==0)   row=1
    if (is.na(column))   column=5
    else if (column==0)   column=1

    z4[i,j] = prob5[row,column]
  }
}

plot_ly(x=x1,y=y1,z=z4,type='surface');




########################################################################


# Variables V1 y V8 de dataIstanbul

varTypeIs = DiscOrCont(dataIstanbul,1);

dataMinIs = as.numeric(apply(dataIstanbul,2,min));
dataMaxIs = as.numeric(apply(dataIstanbul,2,max));


# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrain = sample(1:nrow(dataIstanbul),floor(nrow(dataIstanbul)*0.8));
Ma = c(which.max(dataIstanbul[,1]),which.max(dataIstanbul[,8]))
mi = c(which.min(dataIstanbul[,1]),which.min(dataIstanbul[,8]))
dataTrain = dataIstanbul[c(selectTrain,Ma,mi),];
dataTest = dataIstanbul[-selectTrain,];


# ajuste de variable continua con padre continuo (V8 es padre de V1)
fittMoP = Adjust1ParentCC(dataTrain[,1],dataMinIs[1],dataMaxIs[1],dataTrain[,8],
                          dataMinIs[8],dataMaxIs[8],n=6,maxDegree=7);


# DISCRETIZACIÓN

# discretización con equal width (3 intervalos)
a3 = DiscretizationEqWid(dataTrain[,c(1,8)],3);
dataDisc3 = a3[[1]]
breaks3 = a3[[2]]
# discretización con equal width (5 intervalos)
a5 = DiscretizationEqWid(dataTrain[,c(1,8)],5);
dataDisc5 = a5[[1]]
breaks5 = a5[[2]]

# ajuste discreto (en las filas de la matriz está el hijo)
fitDisc3 = Adjust1ParentDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3);
fitDisc5 = Adjust1ParentDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5);

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1];
area5 = diff(breaks5[1,])[1];
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;



# BIC

# BIC de tMoP
bictMoP = BIC1ParentC(fittMoP[[2]],dataTest[,1],fittMoP[[1]],dataTest[,8],fittMoP[[3]],
                      varType=1);

# BIC discretización con 3 intervalos
logDisc3 = 0;

for (i in 2:ncol(breaks3))
{
  for (j in 2:ncol(breaks3))
  {
    aux = (dataTest[,1] >= breaks3[1,i-1] & dataTest[,1] < breaks3[1,i] &
        dataTest[,8] >= breaks3[2,j-1] & dataTest[,8] < breaks3[2,j]) * log10(prob3[i-1,j-1])

    logDisc3 = logDisc3 + sum(aux);
  }
}
bicDisc3 = logDisc3 - (ncol(breaks3)^2 - 1) / 2*log10(length(dataTest));


# BIC discretización con 5 intervalos
logDisc5 = 0;

for (i in 2:ncol(breaks5))
{
  for (j in 2:ncol(breaks5))
  {
    aux = (dataTest[,1] >= breaks5[1,i-1] & dataTest[,1] < breaks5[1,i] &
        dataTest[,8] >= breaks5[2,j-1] & dataTest[,8] < breaks5[2,j]) * log10(prob5[i-1,j-1])

    logDisc5 = logDisc5 + sum(aux);
  }
}
bicDisc5 = logDisc5 - (ncol(breaks5)^2 - 1) / 2*log10(length(dataTest));

bictMoP
bicDisc3
bicDisc5



# GRÁFICAS

# vectores para dibujar las gráficas
x1 = seq(dataMinIs[1],dataMaxIs[1],length.out=100);
y1 = seq(dataMinIs[8],dataMaxIs[8],length.out=100);


# gráfica original
x = seq(dataMinIs[1],dataMaxIs[1],length.out=20);   # hijo
y = seq(dataMinIs[8],dataMaxIs[8],length.out=20);   # padre
# la función Histograma3D está al final de script
z1 = Histograma3D(dataIstanbul[,1],dataIstanbul[,8],x,y);

plot_ly(x=x[-1],y=y[-1],z=t(z1),type='surface');


# gráfica del ajuste con tMoP
z2 = matrix(,nrow=length(x1),ncol=length(y1));

for (i in 1:length(x1))
{
  z2[i,] = PredicttMoPCC(x1,y1[i],fittMoP[[2]],fittMoP[[1]],fittMoP[[3]]);
}

plot_ly(x=x1,y=y1,z=z2,type='surface');


# gráfica del ajuste discretizado con 3 intervalos
z3 = matrix(,nrow=length(x1),ncol = length(x1));
for (i in 1:length(x1))
{
  for (j in 1:length(y1))
  {
    row = which(x1[i] < breaks3[1,])[1] - 1;
    column = which(y1[j] < breaks3[2,])[1] - 1;
    if (is.na(row))   row=3
    else if (row==0)   row=1
    if (is.na(column))   column=3
    else if (column==0)   column=1

    z3[i,j] = prob3[row,column]
  }
}

plot_ly(x=x1,y=y1,z=z3,type='surface');



# gráfica del ajuste discretizado con 5 intervalos
z4 = matrix(,nrow=length(x1),ncol = length(x1));
for (i in 1:length(x1))
{
  for (j in 1:length(y1))
  {
    row = which(x1[i] < breaks5[1,])[1] - 1;
    column = which(y1[j] < breaks5[2,])[1] - 1;
    if (is.na(row))   row=5
    else if (row==0)   row=1
    if (is.na(column))   column=5
    else if (column==0)   column=1

    z4[i,j] = prob5[row,column]
  }
}

plot_ly(x=x1,y=y1,z=z4,type='surface');








########################################################################


Histograma3D = function(x,y,breaksX,breaksY)
{
  # obtiene una matriz con el histograma en 3d de las variables X e Y


  # crear matriz para probabilidades
  z = matrix(nrow=length(breaksX)-1,ncol=length(breaksY)-1);

  for (i in 2:length(breaksX))
  {
    for (j in 2:length(breaksY))
    {
      # contar el número de datos entre los dos puntos de corte de cada variable
      n = length(which(x>=breaksX[i-1] & x<breaksX[i] & y>=breaksY[j-1] & y<breaksY[j]));
      z[i-1,j-1] = n;
    }
  }

  # calcular la probabilidad
  for (j in 1:ncol(z))
  {
    if (sum(z[,j])>0)   z[,j] = z[,j] / sum(z[,j])
    else   z[,j] = 0
  }

  return(z)
}
############################################################





