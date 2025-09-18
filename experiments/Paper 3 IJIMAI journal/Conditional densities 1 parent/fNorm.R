
# PRUEBA SOBRE APRENDIZAJEDE VARIABLES CONTINUAS CON PADRES CONTINUOS

library(squash);
library(plotly);
library(pracma);
library(ggplot2);
library(plot3D);
library(classInt);

set.seed(2)


# ponemos el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# cambiamos el directorio para cargar las funciones
setwd("../../../core")
source("Adjust1Parent.R");          source("Adjust2Parents.R");      source("BIC.R");
source("PolynomialFit.R");          source("ShowModel.R");           source("Prediction.R");
source("AuxiliarFunctions.R");      source("Adjust2ParentsCases.R");

# ponemos el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("DivergenceKL.R");





########################################################################

# DISTRIBUCIONES NORMALES BIVARIANTES

########################################################################


# X-N(0,1) definida en el intervalo [-3,3]
# Y-N(0,1)  definida en el intervalo [-3,3],    correlación = 1/2
# X|Y-N(y/2,sqrt(3/4))

meanX = 0;         stdX = 1;
meanY = 0;         stdY = 1;
ro = 1/2;

# generar valores aleatorios para X (el padre)
x = rnorm(2000,meanX,stdX);
x = x[x>=-3 & x<=3]

# media y desviación típica de la condicionada Y|X
meanYCond = 0.5*x;
stdYCond = sqrt(3/4);

# inicializar vector con valores de Y
y = rep(0,length(x))

for (i in 1:length(x))
{
  # generar el valor de Y a partir del valor de X
  y[i] = rnorm(1,meanYCond[i],stdYCond);
}

# crear el dataframe con X e Y (hijo en columna 1, padre en columna 2)
data = as.data.frame(cbind(y[y>=-3 & y<=3],x[y>=-3 & y<=3]));
#data = as.data.frame(cbind(y,x));

# seleccionar 1000 datos
data = data[1:1000,];

# valores mínimo y máximo de cada variable
#dataMin = as.numeric(apply(data,2,min));
#dataMax = as.numeric(apply(data,2,max));


# ajuste de variable continua con padre continuo (X es padre de Y)
fittMoP = Adjust1ParentCC(data[,1],-3,3,data[,2],-3,3,n=2,maxDegree=7);


# DISCRETIZACIÓN

# discretización con equal width (3 intervalos)
a3 = DiscretizationEqWid(data,3);
dataDisc3 = a3[[1]]
#breaks3 = a3[[2]]
breaks3 = rbind(seq(-3,3,2),seq(-3,3,2))

# discretización con equal width (5 intervalos)
a5 = DiscretizationEqWid(data,5);
dataDisc5 = a5[[1]]
#breaks5 = a5[[2]]
breaks5 = rbind(seq(-3,3,1.2),seq(-3,3,1.2))

# ajuste discreto (en las filas de la matriz está el hijo)
fitDisc3 = Adjust1ParentDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3);
fitDisc5 = Adjust1ParentDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5);

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1]   # * diff(breaks3[2,])[1]
area5 = diff(breaks5[1,])[1]   # * diff(breaks5[2,])[1]
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;



# DIVERGENCIA KL

# X-N(0,1)     Y|X-N(x/2,sqrt(3/4))
fX = function(x)        exp(-(x-0)^2/(2*1^2)) / (1*sqrt(2*pi)) / 0.9973002
#fYX = function(y,x)     exp(-(y-x/2)^2 / (2*(3/4)) / (sqrt(3/4) * sqrt(2*pi))


fYX = function(y,x)
{
  s = c(seq(-3,-2,0.2),-1.5,-1,0,1,1.5,seq(2,3,0.2))
  select = vector(mode='numeric',length=length(x))
  i = c(0.9583676,0.967664,0.975176,0.9811659,0.9858787,0.9895374,0.9953052,0.9980272,
      0.999468,0.9980272,0.9953052,0.9895374,0.9858787,0.9811659,0.975176,0.967664,0.9583676)
  for (j in 1:length(x))      select[j] = which(x[j]<=s)[1]
  (exp(-(y-x/2)^2 / (2*sqrt(3/4)^2)) / (sqrt(3/4) * sqrt(2*pi))) / i[select]
}



pracma::integral(fX,-3,3)

s = c(seq(-3,-2,0.2),-1.5,-1,0,1,1.5,seq(2,3,0.2))
for (i in 1:length(s))
{
  f3 = function(y)   fYX(y,s[i])
  cat('\n s=',s[i],'  integral = ',pracma::integral(f3,-3,3))
}



DivKLtmopXY(fX,fYX,fittMoP,-3,3)

DivKLDiscretXY(fX,fYX,prob3,breaks3)
DivKLDiscretXY(fX,fYX,prob5,breaks5)



# GRÁFICAS

# vectores para dibujar las gráficas
x1 = seq(-3,3,length.out=200);     y1 = x1;

# gráfica de la distribución normal condicionada
z1 = matrix(,nrow=length(x1),ncol = length(x1));
for (i in 1:length(x1))     {     z1[i,] = dnorm(y1,0.5*x1[i],3/4);     }


## GUARDAR EN PDF
pdf('GraficaCondNorm01-01-orig.pdf',width=7,height=5)

par(mar = c(1,0,1,2))
colors = colorRampPalette(c("darkblue", "darkgreen", "yellow"))
persp3D(x1,y1,t(z1),col=colors(256),theta=130,phi=40,ticktype="detailed",lwd=0.6)
title(main='Original conditional density')

#plot_ly(x=x1,y=y1,z=z1,type='surface');
dev.off()


# gráfica del ajuste con tMoP
z2 = matrix(,nrow=length(x1),ncol = length(x1));
for (i in 1:length(x1))
{
  z2[i,] = PredicttMoPCC(x1,y1[i],fittMoP[[2]],fittMoP[[1]],fittMoP[[3]]);
}

pdf('GraficaCondNorm01-01-tmop.pdf',width=7,height=5)

par(mar = c(1,0,1,2))
colors = colorRampPalette(c("darkblue", "darkgreen", "yellow"))
persp3D(x1,y1,t(z2),col=colors(256),theta=130,phi=40,ticktype="detailed",lwd=0.6)
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


pdf('GraficaCondNorm01-01-disc3.pdf',width=7,height=5)

par(mar = c(1,0,1,2))
colors = colorRampPalette(c("darkblue", "darkgreen", "yellow"))
persp3D(x1,y1,z3,col=colors(256),theta=130,phi=40,ticktype="detailed",lwd=0.6)
title(main='Discretization in 3 intervals per variable')

#plot_ly(x=x1,y=y1,z=z3,type='surface');
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


pdf('GraficaCondNorm01-01-disc5.pdf',width=7,height=5)

par(mar = c(1,0,1,2))
colors = colorRampPalette(c("darkblue", "darkgreen", "yellow"))
persp3D(x1,y1,z4,col=colors(256),theta=130,phi=40,ticktype="detailed",lwd=0.6)
title(main='Discretization in 5 intervals per variable')

#plot_ly(x=x1,y=y1,z=z4,type='surface');
dev.off()









########################################################################

set.seed(2)

# X-N(0,1) definida en el intervalo [-3,3]
# Y-N(1,3) definida en el intervalo [-9,11],    correlación = 3/4
# Y|X-N(1+9x/4,3*sqrt(7)/4)

meanX = 0;         stdX = 1;
meanY = 1;         stdY = 3;
ro = 3/4;

# generar valores aleatorios para X
x = rnorm(2000,meanX,stdX);
x = x[x>=-3 & x<=3]

# media y desviación típica de la condicionada Y|X
meanYCond = 1+9*x/4;
stdYCond = 3*sqrt(7)/4;

# inicializar vector con valores de Y
y = rep(0,length(x))

for (i in 1:length(x))
{
  # generar el valor de Y a partir del valor de X
  y[i] = rnorm(1,meanYCond[i],stdYCond);
}

# crear el dataframe con X e Y (X es padre de Y)
data = as.data.frame(cbind(y[y>=-9 & y<=11],x[y>=-9 & y<=11]));
#data = as.data.frame(cbind(y,x));

# seleccionar 1000 datos
data = data[1:1000,];


# valores mínimo y máximo de cada variable
dataMin = as.numeric(apply(data,2,min));
dataMax = as.numeric(apply(data,2,max));


# ajuste de variable continua con padre continuo (X es padre de Y)
fittMoP = Adjust1ParentCC(data[,1],-9,11,data[,2],-3,3,n=2,maxDegree=7)




# DISCRETIZACIÓN

# discretización con equal width (3 intervalos)
a3 = DiscretizationEqWid(data,3);
dataDisc3 = a3[[1]]
#breaks3 = a3[[2]]
breaks3 = rbind(seq(-9,11,length.out=4),seq(-3,3,2))

# discretización con equal width (3 intervalos)
a5 = DiscretizationEqWid(data,5);
dataDisc5 = a5[[1]]
#breaks5 = a5[[2]]
breaks5 = rbind(seq(-9,11,4),seq(-3,3,1.2))

# ajuste discreto (en las filas de la matriz está el hijo)
fitDisc3 = Adjust1ParentDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3);
fitDisc5 = Adjust1ParentDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5);

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1]   # * diff(breaks3[2,])[1]
area5 = diff(breaks5[1,])[1]   # * diff(breaks5[2,])[1]
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;



# DIVERGENCIA KL

# X-N(0,1)     Y|X-N(1+9x/4,3*sqrt(7)/4)
fX = function(x)        exp(-(x-0)^2/(2*1^2)) / (1*sqrt(2*pi)) / 0.9973002
#fYX = function(y,x)     exp(-(y-1-9*x/4)^2 / (2*9*7/16)) / ((3/4) * sqrt(2*pi*7))

fYX = function(y,x)
{
  s = c(seq(-3,-2,0.2),-1,1,seq(2,3,0.2))
  select = vector(mode='numeric',length=length(x))
  i = c(0.9492731,0.968883,0.981754,0.9897803,0.9945355,0.9972121,0.999953,0.999953,
        0.9972121,0.9945355,0.9897803,0.981754,0.968883,0.9492731)
  for (j in 1:length(x))      select[j] = which(x[j]<=s)[1]
  (exp(-(y-1-9*x/4)^2 / (2*9*7/16)) / ((3/4) * sqrt(2*pi*7))) / i[select]
}


pracma::integral(fX,-3,3)

s = seq(-3,3,0.5)
s = c(seq(-3,-2,0.2),-1,1,seq(2,3,0.2))
for (i in 1:length(s))
{
  f3 = function(y)   fYX(y,s[i])
  cat('\n s=',s[i],'  integral = ',pracma::integral(f3,-9,11))
}



DivKLtmopXY(fX,fYX,fittMoP,-9,11)

DivKLDiscretXY(fX,fYX,prob3,breaks3)
DivKLDiscretXY(fX,fYX,prob5,breaks5)






# GRÁFICAS

# vectores para dibujar las gráficas
x1 = seq(-3,3,length.out=1000);     y1 = seq(-9,11,length.out=1000);


# gráfica de la distribución normal condicionada
z1 = matrix(,nrow=length(x1),ncol=length(y1));
for (i in 1:length(x1))     {     z1[i,] = dnorm(y1,1+9*x1[i]/4,3*sqrt(7)/4);     }

plot_ly(x=x1,y=y1,z=t(z1),type='surface');


# gráfica del ajuste con tMoP
z2 = matrix(,nrow=length(x1),ncol=length(y1));
for (i in 1:length(x1))
{
  z2[i,] = PredicttMoPCC(y1,x1[i],fittMoP[[2]],fittMoP[[1]],fittMoP[[3]]);
}

plot_ly(x=x1,y=y1,z=t(z2),type='surface');



# gráfica del ajuste discretizado con 3 intervalos
z3 = matrix(,nrow=length(x1),ncol = length(y1));
for (i in 1:length(x1))
{
  for (j in 1:length(y1))
  {
    row = which(x1[i] < breaks3[2,])[1] - 1;
    column = which(y1[j] < breaks3[1,])[1] - 1;
    if (is.na(row))   row=3
    else if (row==0)   row=1
    if (is.na(column))   column=3
    else if (column==0)   column=1

    z3[i,j] = prob3[row,column]
  }
}

plot_ly(x=x1,y=y1,z=t(z3),type='surface');



# gráfica del ajuste discretizado con 5 intervalos
z4 = matrix(,nrow=length(x1),ncol = length(y1));
for (i in 1:length(x1))
{
  for (j in 1:length(y1))
  {
    row = which(x1[i] < breaks5[2,])[1] - 1;
    column = which(y1[j] < breaks5[1,])[1] - 1;
    if (is.na(row))   row=5
    else if (row==0)   row=1
    if (is.na(column))   column=5
    else if (column==0)   column=1

    z4[i,j] = prob5[row,column]
  }
}

plot_ly(x=x1,y=y1,z=t(z4),type='surface');









########################################################################

set.seed(2)

# X-N(1,3) definida en el intervalo [-5,7]
# Y-N(0,5) definida en el intervalo [-12,12],    correlación = 2/3
# Y|X-N(10(x-1)/9,5*sqrt(5)/3)

meanX = 1;         stdX = 3;
meanY = 0;         stdY = 5;
ro = 2/3;

# generar valores aleatorios para X
x = rnorm(2000,meanX,stdX);
x = x[x>=-5 & x<=7]

# media y desviación típica de la condicionada Y|X
meanYCond = 10*(x-1)/9;
stdYCond = 5*sqrt(5)/3;

# inicializar vector con valores de Y
y = rep(0,length(x))

for (i in 1:length(x))
{
  # generar el valor de Y a partir del valor de X
  y[i] = rnorm(1,meanYCond[i],stdYCond);
}

# crear el dataframe con X e Y (X es padre de Y)
data = as.data.frame(cbind(y[y>=-12 & y<=12],x[y>=-12 & y<=12]));
#data = as.data.frame(cbind(y,x));

# seleccionar 1000 datos
data = data[1:1000,];



# valores mínimo y máximo de cada variable
dataMin = as.numeric(apply(data,2,min));
dataMax = as.numeric(apply(data,2,max));


# ajuste de variable continua con padre continuo (X es padre de Y)
fittMoP = Adjust1ParentCC(data[,1],-12,12,data[,2],-5,7,n=2,maxDegree=7);



# DISCRETIZACIÓN

# discretización con equal width (3 intervalos)
a3 = DiscretizationEqWid(data,3);
dataDisc3 = a3[[1]]
#breaks3 = a3[[2]]
breaks3 = rbind(seq(-12,12,8),seq(-5,7,4))

# discretización con equal width (3 intervalos)
a5 = DiscretizationEqWid(data,5);
dataDisc5 = a5[[1]]
#breaks5 = a5[[2]]
breaks5 = rbind(seq(-12,12,4.8),seq(-5,7,2.4))

# ajuste discreto (en las filas de la matriz está el hijo)
fitDisc3 = Adjust1ParentDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3);
fitDisc5 = Adjust1ParentDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5);

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1]   # * diff(breaks3[2,])[1]
area5 = diff(breaks5[1,])[1]   # * diff(breaks5[2,])[1]
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;



# DIVERGENCIA KL

# X-N(1,3),    Y|X-N(10(x-1)/9,5*sqrt(5)/3)
fX = function(x)        (exp(-(x-1)^2/(2*3^2)) / (3*sqrt(2*pi))) / 0.9544997
#fYX = function(y,x)     exp(-(y-10*(x-1)/9)^2 / (2*25*5/9)) / ((5/3) * sqrt(2*pi*5))

fYX = function(y,x)
{
  s = c(seq(-5,-2,0.3),seq(-1.5,3.5,0.5),seq(4,7,0.3))
  select = vector(mode='numeric',length=length(x))
  i = c(0.9237966,0.9358102,0.946297,0.9553781,0.9631792,0.9698272,0.9754471,0.9801598,
        0.9840796,0.9873132,0.9899583,0.9932937,0.9955826,0.9970979,0.9980426,0.9985556,
        0.9987178,0.9985556,0.9980426,0.9970979,0.9955826,0.9932937,0.9899583,0.9873132,
        0.9840796,0.9801598,0.9754471,0.9698272,0.9631792,0.9553781,0.946297,0.9358102,
        0.9237966);
  for (j in 1:length(x))      select[j] = which(x[j]<=s)[1]
  (exp(-(y-10*(x-1)/9)^2 / (2*25*5/9)) / ((5/3) * sqrt(2*pi*5))) / i[select]
}



pracma::integral(fX,-5,7)

s = seq(-5,7,0.5)
s = c(seq(-5,-2,0.3),seq(-1.5,3.5,0.5),seq(4,7,0.3))
for (i in 1:length(s))
{
  f3 = function(y)   fYX(y,s[i])
  cat('\n s=',s[i],'  integral = ',pracma::integral(f3,-12,12))
}



DivKLtmopXY(fX,fYX,fittMoP,-12,12)

DivKLDiscretXY(fX,fYX,prob3,breaks3)
DivKLDiscretXY(fX,fYX,prob5,breaks5)






# GRÁFICAS

# vectores para dibujar las gráficas
x1 = seq(-5,7,length.out=1000);     y1 = seq(-12,12,length.out=1000);


# gráfica de la distribución normal condicionada
z1 = matrix(,nrow=length(x1),ncol=length(y1));
for (i in 1:length(x1))     {     z1[i,] = dnorm(y1,10*(x1[i]-1)/9,5*sqrt(5)/3);     }

plot_ly(x=x1,y=y1,z=z1,type='surface');


# gráfica del ajuste con tMoP
z2 = matrix(,nrow=length(x1),ncol=length(y1));
for (i in 1:length(x1))
{
  z2[i,] = PredicttMoPCC(y1,x1[i],fittMoP[[2]],fittMoP[[1]],fittMoP[[3]]);
}

plot_ly(x=x1,y=y1,z=z2,type='surface');



# gráfica del ajuste discretizado con 3 intervalos
z3 = matrix(,nrow=length(x1),ncol = length(y1));
for (i in 1:length(x1))
{
  for (j in 1:length(y1))
  {
    row = which(x1[i] < breaks3[2,])[1] - 1;
    column = which(y1[j] < breaks3[1,])[1] - 1;
    if (is.na(row))   row=3
    else if (row==0)   row=1
    if (is.na(column))   column=3
    else if (column==0)   column=1

    z3[i,j] = prob3[row,column]
  }
}

plot_ly(x=x1,y=y1,z=t(z3),type='surface');



# gráfica del ajuste discretizado con 5 intervalos
z4 = matrix(,nrow=length(x1),ncol = length(y1));
for (i in 1:length(x1))
{
  for (j in 1:length(y1))
  {
    row = which(x1[i] < breaks5[2,])[1] - 1;
    column = which(y1[j] < breaks5[1,])[1] - 1;
    if (is.na(row))   row=5
    else if (row==0)   row=1
    if (is.na(column))   column=5
    else if (column==0)   column=1

    z4[i,j] = prob5[row,column]
  }
}

plot_ly(x=x1,y=y1,z=t(z4),type='surface');















