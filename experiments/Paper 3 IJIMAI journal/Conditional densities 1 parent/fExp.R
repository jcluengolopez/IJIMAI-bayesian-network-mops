
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

# DISTRIBUCIONES EXPONENCIALES BIVARIANTES

########################################################################


# X-exp(0.5) definida en el intervalo [0,3]
# Y-exp(1)  definida en el intervalo [0,3]         correlación = 1/2

kX = 0.5;         kY = 1;
ro = 1/2;

# crear la función de densidad marginal y la condicionada (Y condicionado a X)
fX = function(x)     kX * exp(-kX*x)
fXY = function(x,y)    (kX+ro*y) * exp(-(kX+ro*y)*x) * ((kY+ro*x)/kY-ro/(kY*(kX+ro*y)))

# función de distribución de Y condicionada a X
FYX = function(y,x)   -(kX+ro*y)/kX * exp(-(kY+ro*x)*y) + 1



# generar valores aleatorios para X y seleccionar los menores que 6
set.seed(2)
x = rexp(2000,kX);
x = x[x<=3]
x = x[x<=10]
length(x) / 2000

# inicializar vector con valores de Y
y = rep(0,length(x))

for (i in 1:length(x))
{
  # generar un valor aleatorio de la probabilidad
  p = runif(1,0,1)
  
  # función condicionada en el valor concreto de x para resolver la ecuación
  fEq = function(y)   FYX(y,x[i]) - p
  
  # resolver la ecuación usando la función de distribución (debe haber una única solución)
  y[i] = uniroot(fEq,c(0, 1E6))$root;
}

length(y[y<=3]) / length(x)

# crear el dataframe con X e Y con valores <3 (hijo en columna 1, padre en columna 2)
data = as.data.frame(cbind(y[y<=3],x[y<=3]));
#data = as.data.frame(cbind(y[y<=10],x[y<=10]));


# seleccionar 1000 datos
data = data[1:1000,];



# AJUSTE TMOP

# ajuste de variable continua con padre continuo (X es padre de Y)
fittMoP = Adjust1ParentCC(data[,1],0,3,data[,2],0,3,n=2,maxDegree=7);
#fittMoP = Adjust1ParentCC(data[,1],0,8,data[,2],0,8,n=2,maxDegree=7);




# DISCRETIZACIÓN

# discretización con equal width (3 intervalos)
a3 = DiscretizationEqWid(data,3);
dataDisc3 = a3[[1]]
breaks3 = a3[[2]]

# discretización con equal width (3 intervalos)
a4 = DiscretizationEqWid(data,4);
dataDisc4 = a4[[1]]
#breaks4 = a4[[2]]

# discretización con equal width (3 intervalos)
a5 = DiscretizationEqWid(data,5);
dataDisc5 = a5[[1]]
#breaks5 = a5[[2]]

# puntos de corte en el intervalo [0,3]
breaks3 = rbind(0:3,0:3)
br4 = rbind(seq(0,3,0.75),seq(0,3,0.75))
breaks5 = rbind(seq(0,3,0.6),seq(0,3,0.6))
#breaks5 = rbind(seq(0,10,2),seq(0,10,2))


# ajuste discreto (en las filas de la matriz está el hijo)
fitDisc3 = Adjust1ParentDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3);
fitDisc4 = Adjust1ParentDD(dataDisc4[,1],1:4,dataDisc4[,2],1:4);
fitDisc5 = Adjust1ParentDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5);

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
prob3 = fitDisc3[[1]] / diff(breaks3[1,])[1];
prob5 = fitDisc5[[1]] / diff(breaks5[1,])[1];
#prob5 = fitDisc5[[1]];
#colnames(prob3) = c('X_1','X_2','X_3')
#colnames(prob5) = c('X_1','X_2','X_3','X_4','X_5')
#rownames(prob3) = c('Y_1','Y_2','Y_3')
#rownames(prob5) = c('Y_1','Y_2','Y_3','Y_4','Y_5')





# DIVERGENCIA KL

# X-exp(0.5) padre en  [0,3],    Y-exp(1) hija en  [0,3],    correlación = 1/2

#integrate(f,0,3)$value = 0.7768698
fX = function(x)    (0.5 * exp(-0.5*x)) / 0.7768698

#fX = function(x)    0.5 * exp(-0.5*x)
#fYX = function(y,x)  (1+1/2*x) * exp(-(1+1/2*x)*y) * ((0.5+1/2*y)/0.5-(1/2)/(0.5*(1+1/2*x)))


fYX = function(y,x)
{
  #f = function(y)   (1+1/2*x) * exp(-(1+1/2*x)*y) * ((0.5+1/2*y)/0.5-(1/2)/(0.5*(1+1/2*x)))
  #i = integrate(f,0,3)$value
  #cat('\n f =',f(y))
  s = c(seq(0,0.45,0.05),seq(0.5,1.5,0.1),1.7,2,2.5,3)
  select = vector(mode='numeric',length=length(x))
  i = c(0.8008517,0.8152415,0.8285915,0.8409769,0.8524673,0.8631275,0.8730175,0.8821928,
      0.8907051,0.8986024,0.905929,0.9190324,0.9303105,0.9400177,0.9483727,0.955564,
      0.9617536,0.967081,0.9716664,0.975613,0.9790099,0.9844502,0.990085,0.9953165,0.9977877)
  for (j in 1:length(x))      select[j] = which(x[j]<=s)[1]
  
  ((1+1/2*x) * exp(-(1+1/2*x)*y)*((0.5+1/2*y)/0.5-(1/2)/(0.5*(1+1/2*x)))) / i[select]
  #cat('\n select = ',select)
}

s = seq(0,10,0.5)
for (i in 1:length(s))
{
  f3 = function(y)   fYX(y,s[i])
  cat('\n s=',s[i],'  integral = ',pracma::integral(f3,0,10))
}




s = seq(0,3,0.1)
for (i in 1:length(s))
{
  x = s[i]
  f1 = function(y)   (1+1/2*x) * exp(-(1+1/2*x)*y)*((0.5+1/2*y)/0.5-(1/2)/(0.5*(1+1/2*x)))
  #i = integral(f1,0,3)
  #f = function(y)   ((1+1/2*x)*exp(-(1+1/2*x)*y)*((0.5+1/2*y)/0.5-(1/2)/(0.5*(1+1/2*x)))) / i
  
  k = pracma::integral(f1,0,3)
  cat('\n x =',x,'   integral =',k)
}



DivKLtmopXY(fX,fYX,fittMoP,0,3)

DivKLDiscretXY(fX,fYX,prob3,breaks3)
DivKLDiscretXY(fX,fYX,prob5,breaks5)






# GRÁFICAS

# vectores para dibujar las gráficas
x1 = seq(0,3,length.out=100);     y1 = seq(0,3,length.out=100);

# gráfica de la distribución exponencial condicionada
z1 = matrix(,nrow=length(x1),ncol=length(y1));
for (i in 1:length(x1))     {     z1[i,] = fYX(y1,x1[i])     }


pdf('GraficaCondExp05-1-orig.pdf',width=7,height=5)

par(mar = c(1,0,0.5,2))
colors = colorRampPalette(c("darkblue", "darkgreen", "yellow"))
persp3D(x1,y1,t(z1),col=colors(256),theta=30,phi=30,ticktype="detailed",lwd=0.6)

#plot_ly(x=x1,y=y1,z=t(z1),type='surface');
dev.off()



# gráfica del ajuste con tMoP (X es padre de Y)
z2 = matrix(,nrow=length(x1),ncol=length(x1));
for (i in 1:length(x1))
{
  z2[i,] = PredicttMoPCC(y1,rep(x1[i],length(x1)),fittMoP[[2]],fittMoP[[1]],fittMoP[[3]]);
}

z2[,1] = z2[,2]


pdf('GraficaCondExp05-1-tmop.pdf',width=7,height=5)

par(mar = c(1,0,0.5,2))
colors = colorRampPalette(c("darkblue", "darkgreen", "yellow"))
persp3D(x1,y1,t(z2),col=colors(256),theta=30,phi=30,ticktype="detailed",lwd=0.6)

#plot_ly(x=x1,y=y1,z=t(z2),type='surface');
dev.off()



# gráfica del ajuste discretizado con 3 intervalos
z3 = matrix(,nrow=length(x1),ncol = length(y1));
for (i in 1:length(y1))
{
  for (j in 1:length(x1))
  {
    row = which(y1[i] < breaks3[1,])[1] - 1;
    column = which(x1[j] < breaks3[2,])[1] - 1;
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
for (i in 1:length(y1))
{
  for (j in 1:length(x1))
  {
    row = which(y1[i] < breaks5[1,])[1] - 1;
    column = which(x1[j] < breaks5[2,])[1] - 1;
    if (is.na(row))   row=5
    else if (row==0)   row=1
    if (is.na(column))   column=5
    else if (column==0)   column=1
    
    z4[i,j] = prob5[row,column]
  }
}


pdf('GraficaCondExp05-1-disc.pdf',width=7,height=5)

par(mar = c(1,0,0.5,2))
colors = colorRampPalette(c("darkblue", "darkgreen", "yellow"))
persp3D(x1,y1,z4,col=colors(256),theta=30,phi=30,ticktype="detailed",lwd=0.6)

#plot_ly(x=x1,y=y1,z=t(z4),type='surface');
dev.off()







