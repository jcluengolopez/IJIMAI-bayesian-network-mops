
library(MoTBFs)
library(pracma)
library(polynom)
library(NlcOptim)


# poner el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# carga la función de la divergencia y Shenoy
source('DivergenceKL.R');          source('ShenoyFit.R')


# ponemos el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# cambiamos el directorio para cargar las funciones
setwd("../../../core")
source("CreateModelTAN.R");     source("MutualInfo.R");            source("Tree.R");
source("AuxiliarFunctions.R");  source("PolynomialFit.R");         source("ShowModel.R");
source("Adjust1Parent.R");      source("Adjust2Parents.R");        source("BIC.R");
source("Prediction.R");         source("Adjust2ParentsCases.R");   source("TAN.R");
source("CreateModelNB.R");      source("NaiveBayes.R");            source("NBTAN.R");
source("StructureModel.R");     source("Simulation.R");

# ponemos el directorio de este script y lo cambiamos para cargar los datos
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


####################################################################

# DISTRIBUCIÓN WEIBULL con alpha=1 y lambda=1 en el intervalo [0,5]
fW11 = function(x)   1 * x^(1-1) / 1^1 * exp(-(x/1)^1);
iW11 = pracma::integral(fW11,0,5);
fWeibull11 = function(x)   (1 * x^(1-1) / 1^1 * exp(-(x/1)^1)) / iW11;

# DISTRIBUCIÓN WEIBULL con alpha=1.5 y lambda=1 en el intervalo [0,3]
fW151 = function(x)   1.5 * x^(1.5-1) / 1^1.5 * exp(-(x/1)^1.5);
iW151 = pracma::integral(fW151,0,3);
fWeibull151 = function(x)   (1.5 * x^(1.5-1) / 1^1.5 * exp(-(x/1)^1.5)) / iW151;

# DISTRIBUCIÓN WEIBULL con alpha=2 y lambda=2 en el intervalo [0,5]
fW22 = function(x)   2 * x^(2-1) / 2^2 * exp(-(x/2)^2);
iW22 = pracma::integral(fW22,0,5);
fWeibull22 = function(x)   (2 * x^(2-1) / 2^2 * exp(-(x/2)^2)) / iW22;

# DISTRIBUCIÓN WEIBULL con alpha=6 y lambda=3 en el intervalo [0,4]
fW63 = function(x)   6 * x^(6-1) / 3^6 * exp(-(x/3)^6);
iW63 = pracma::integral(fW63,0,4);
fWeibull63 = function(x)   (6 * x^(6-1) / 3^6 * exp(-(x/3)^6)) / iW63;

# DISTRIBUCIÓN WEIBULL con alpha=8 y lambda=4 en el intervalo [0,5]
fW84 = function(x)   8 * x^(8-1) / 4^8 * exp(-(x/4)^8);
iW84 = pracma::integral(fW84,0,5);
fWeibull84 = function(x)   (8 * x^(8-1) / 4^8 * exp(-(x/4)^8)) / iW84;




####################################################################

# Shenoy
# - Número trozos: 1,2,3,4
# - Grado: 1 trozo (7,4,3), 2 trozos (7,4,3), 3 trozos (7,4,3), 4 trozos (7,4,3,2)

# tmop y Pérez-Bernabé
# - Grados: 7, 4 y 3

# Otros métodos:
# - Tamaño muestra: 50,100,1000



# data frame para los resultados
distribucion1 = rep('weibull(1,1)',31);
intervalo1 = rep('[0,5]',31);
distribucion2 = rep('weibull(1.5,1)',31);
intervalo2 = rep('[0,3]',31);
distribucion3 = rep('weibull(2,2)',31);
intervalo3 = rep('[0,5]',31);
distribucion4 = rep('weibull(6,3)',31);
intervalo4 = rep('[0,4]',31);
distribucion5 = rep('weibull(8,4)',31);
intervalo5 = rep('[0,5]',31);
metodo = c(rep('Shenoy',13),rep('Perez-Bernabe',9),rep('tMoPs',9));
nInter = c(rep(1,3),rep(2,3),rep(3,3),rep(4,4),rep(1,9),rep('1-3',9));
grado = c(rep(c(7,4,3),4),2,rep(c(7,4,3),6));
n = c(rep(NA,13),rep(c(50,50,50,100,100,100,1000,1000,1000),2));

r1 = data.frame(distribucion1,intervalo1,metodo,nInter,grado,n,rep(0,31))
colnames(r1) = c('f','intervalo','Método','nInter','Grado','Muestra','Div-KL');
r2 = data.frame(distribucion2,intervalo2,metodo,nInter,grado,n,rep(0,31))
colnames(r2) = c('f','intervalo','Método','nInter','Grado','Muestra','Div-KL');
r3 = data.frame(distribucion3,intervalo3,metodo,nInter,grado,n,rep(0,31))
colnames(r3) = c('f','intervalo','Método','nInter','Grado','Muestra','Div-KL');
r4 = data.frame(distribucion4,intervalo4,metodo,nInter,grado,n,rep(0,31))
colnames(r4) = c('f','intervalo','Método','nInter','Grado','Muestra','Div-KL');
r5 = data.frame(distribucion5,intervalo5,metodo,nInter,grado,n,rep(0,31))
colnames(r5) = c('f','intervalo','Método','nInter','Grado','Muestra','Div-KL');

resultados = rbind(r1,r2,r3,r4,r5)




####################################################################

# MÉTODO SHENOY (WEIBULL(1,1))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fWeibull11,7,0,5,1);
fitShenoy14 = ShenoyFit(fWeibull11,4,0,5,1);
fitShenoy13 = ShenoyFit(fWeibull11,3,0,5,1);
fitShenoy27 = ShenoyFit(fWeibull11,7,0,5,2);
fitShenoy24 = ShenoyFit(fWeibull11,4,0,5,2);
fitShenoy23 = ShenoyFit(fWeibull11,3,0,5,2);
fitShenoy37 = ShenoyFit(fWeibull11,7,0,5,3);
fitShenoy34 = ShenoyFit(fWeibull11,4,0,5,3);
fitShenoy33 = ShenoyFit(fWeibull11,3,0,5,3);
fitShenoy47 = ShenoyFit(fWeibull11,7,0,5,4);
fitShenoy44 = ShenoyFit(fWeibull11,4,0,5,4);
fitShenoy43 = ShenoyFit(fWeibull11,3,0,5,4);
fitShenoy42 = ShenoyFit(fWeibull11,2,0,5,4);



divShenoy17 = DivergenceKL2(fWeibull11,fitShenoy17[[1]],0,5);
divShenoy14 = DivergenceKL2(fWeibull11,fitShenoy14[[1]],0,5);
divShenoy13 = DivergenceKL2(fWeibull11,fitShenoy13[[1]],0,5);

sp2 = seq(0,5,length.out=3);
divShenoy27 = DivergenceKL2(fWeibull11,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fWeibull11,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fWeibull11,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fWeibull11,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fWeibull11,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fWeibull11,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,5,length.out=4);
divShenoy37 = DivergenceKL2(fWeibull11,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fWeibull11,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fWeibull11,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fWeibull11,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fWeibull11,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fWeibull11,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fWeibull11,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fWeibull11,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fWeibull11,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,5,length.out=5);
divShenoy47 = DivergenceKL2(fWeibull11,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull11,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull11,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull11,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fWeibull11,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull11,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull11,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull11,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fWeibull11,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull11,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull11,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull11,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fWeibull11,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull11,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull11,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull11,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ  (WEIBULL(1,1))
# tiene un intervalo, grados 7, 4, 3 y tamaño de muestra: 50, 100, 1000


# ajustes para distintos tamaños de muestra. Repetir cada uno con 10 muestras
divPB1 = rep(0,3);        divPB2 = rep(0,3);        divPB3 = rep(0,3);
p = c(8,5,4)

# 50
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei1 = rweibull(50,1,1);           xWei1 = xWei1[xWei1<=5];
    
    fitPB1 = univMoTBF(xWei1,POTENTIAL_TYPE="MOP",evalRange=c(0,5),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fWeibull11,fitPBpoly1,0,5);
  }
}

fW = function(x)   fWeibull11(x) / integral(fWeibull11,0,5)

# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei2 = rweibull(100,1,1);           xWei2 = xWei2[xWei2<=5];
    
    fitPB2 = univMoTBF(xWei2,POTENTIAL_TYPE="MOP",evalRange=c(0,5),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fWeibull11,fitPBpoly2,0,5);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei3 = rweibull(1000,1,1);           xWei3 = xWei3[xWei3<=5];
    
    fitPB3 = univMoTBF(xWei3,POTENTIAL_TYPE="MOP",evalRange=c(0,5),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fWeibull11,fitPBpoly3,0,5);
    
    #divPB3[j] = divPB3[j] + DivergenceKL(fW,fitPBpoly1,0,5);
    #cat('\n KL = ',divPB3[1])
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (WEIBULL(1,1))
# tiene un intervalo, grados 7, 4, 3 y tamaño de muestra: 50, 100, 1000


# ajustes para distintos tamaños de muestra. Repetir cada uno con 10 muestras
divtmop1 = rep(0,3);        divtmop2 = rep(0,3);        divtmop3 = rep(0,3);
grados = c(7,4,3)

# 50
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei1 = rweibull(50,1,1);           xWei1 = xWei1[xWei1<=5];
    
    fittmop1 = PolynomialFit(xWei1,maxDegree=grados[j],0,5,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fWeibull11,fittmop1,0,5);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei2 = rweibull(100,1,1);           xWei2 = xWei2[xWei2<=5];
    
    fittmop2 = PolynomialFit(xWei2,maxDegree=grados[j],0,5,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fWeibull11,fittmop2,0,5);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei3 = rweibull(1000,1,1);           xWei3 = xWei3[xWei3<=5];
    
    fittmop3 = PolynomialFit(xWei3,maxDegree=grados[j],0,5,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fWeibull11,fittmop3,0,5);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divtmop1 = divtmop1 / 10;       divtmop2 = divtmop2 / 10;       divtmop3 = divtmop3 / 10;




# divergencias en la columna de resultados
resultados$`Div-KL`[1:31] = c(divShenoy17,divShenoy14,divShenoy13,divShenoy27,divShenoy24,
                              divShenoy23,divShenoy37,divShenoy34,divShenoy33,divShenoy47,
                              divShenoy44,divShenoy43,divShenoy42,divPB1,divPB2,divPB3,
                              divtmop1,divtmop2,divtmop3)




####################################################################

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (WEIBULL(1,1))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,1.25,length.out=100);              x2 = seq(1.25,2.5,length.out=100);
x3 = seq(2.5,3.75,length.out=100);            x4 = seq(3.75,5,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xWei3,POTENTIAL_TYPE="MOP",evalRange=c(0,5),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xWei3,maxDegree=7,0,5,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fWeibull11,c(0,5),col='#0072B2',type='l',xlim = c(0,5),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')








####################################################################

# MÉTODO SHENOY (WEIBULL(1.5,1))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fWeibull151,7,0,3,1);
fitShenoy14 = ShenoyFit(fWeibull151,4,0,3,1);
fitShenoy13 = ShenoyFit(fWeibull151,3,0,3,1);
fitShenoy27 = ShenoyFit(fWeibull151,7,0,3,2);
fitShenoy24 = ShenoyFit(fWeibull151,4,0,3,2);
fitShenoy23 = ShenoyFit(fWeibull151,3,0,3,2);
fitShenoy37 = ShenoyFit(fWeibull151,7,0,3,3);
fitShenoy34 = ShenoyFit(fWeibull151,4,0,3,3);
fitShenoy33 = ShenoyFit(fWeibull151,3,0,3,3);
fitShenoy47 = ShenoyFit(fWeibull151,7,0,3,4);
fitShenoy44 = ShenoyFit(fWeibull151,4,0,3,4);
fitShenoy43 = ShenoyFit(fWeibull151,3,0,3,4);
fitShenoy42 = ShenoyFit(fWeibull151,2,0,3,4);



divShenoy17 = DivergenceKL2(fWeibull151,fitShenoy17[[1]],0,3);
divShenoy14 = DivergenceKL2(fWeibull151,fitShenoy14[[1]],0,3);
divShenoy13 = DivergenceKL2(fWeibull151,fitShenoy13[[1]],0,3);

sp2 = seq(0,3,length.out=3);
divShenoy27 = DivergenceKL2(fWeibull151,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fWeibull151,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fWeibull151,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fWeibull151,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fWeibull151,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fWeibull151,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,3,length.out=4);
divShenoy37 = DivergenceKL2(fWeibull151,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fWeibull151,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fWeibull151,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fWeibull151,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fWeibull151,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fWeibull151,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fWeibull151,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fWeibull151,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fWeibull151,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,3,length.out=5);
divShenoy47 = DivergenceKL2(fWeibull151,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull151,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull151,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull151,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fWeibull151,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull151,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull151,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull151,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fWeibull151,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull151,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull151,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull151,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fWeibull151,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull151,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull151,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull151,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ  (WEIBULL(1.5,1))
# tiene un intervalo, grados 7, 4, 3 y tamaño de muestra: 50, 100, 1000


# ajustes para distintos tamaños de muestra. Repetir cada uno con 10 muestras
divPB1 = rep(0,3);        divPB2 = rep(0,3);        divPB3 = rep(0,3);
p = c(8,5,4)

# 50
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei1 = rweibull(50,1.5,1);           xWei1 = xWei1[xWei1<=3];
    
    fitPB1 = univMoTBF(xWei1,POTENTIAL_TYPE="MOP",evalRange=c(0,3),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fWeibull151,fitPBpoly1,0,3);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei2 = rweibull(100,1.5,1);           xWei2 = xWei2[xWei2<=3];
    
    fitPB2 = univMoTBF(xWei2,POTENTIAL_TYPE="MOP",evalRange=c(0,3),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fWeibull151,fitPBpoly2,0,3);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei3 = rweibull(1000,1.5,1);           xWei3 = xWei3[xWei3<=3];
    
    fitPB3 = univMoTBF(xWei3,POTENTIAL_TYPE="MOP",evalRange=c(0,3),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fWeibull151,fitPBpoly3,0,3);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (WEIBULL(1.5,1))
# tiene un intervalo, grados 7, 4, 3 y tamaño de muestra: 50, 100, 1000


# ajustes para distintos tamaños de muestra. Repetir cada uno con 10 muestras
divtmop1 = rep(0,3);        divtmop2 = rep(0,3);        divtmop3 = rep(0,3);
grados = c(7,4,3)

# 50
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    # HE CAMBIADO UN POCO LA SEMILLA PORQUE EN UNA DE ELLAS SALE COEFICIENTES
    # MUY GRANDES QUE DAN PROBLEMAS AL CALCULAR LA INTEGRAL DE LA DIVERGENCIA
    set.seed(i+j-4)
    xWei1 = rweibull(50,1.5,1);           xWei1 = xWei1[xWei1<=3];
    
    fittmop1 = PolynomialFit(xWei1,maxDegree=grados[j],0,3,1,2);
    #cat('\n \n \n  ajuste = ',i);  print(fittmop1)
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fWeibull151,fittmop1,0,3);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei2 = rweibull(100,1.5,1);           xWei2 = xWei2[xWei2<=3];
    
    fittmop2 = PolynomialFit(xWei2,maxDegree=grados[j],0,3,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fWeibull151,fittmop2,0,3);
  }
}



# 1000
s = c(2:4,8:15,18:20,22:24,26:33,35:38)
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    # HE CAMBIADO UN POCO LA SEMILLA PORQUE EN UNA DE ELLAS SALE COEFICIENTES
    # MUY GRANDES QUE DAN PROBLEMAS AL CALCULAR LA INTEGRAL DE LA DIVERGENCIA
    set.seed(s[i+j-1])
    xWei3 = rweibull(1000,1.5,1);           xWei3 = xWei3[xWei3<=3];
    
    fittmop3 = PolynomialFit(xWei3,maxDegree=grados[j],0,3,1,2);
    #cat('\n \n \n  semilla = ',i+j+33);
    #cat('\n \n \n  ajuste = ',i);  print(fittmop3)
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fWeibull151,fittmop3,0,3);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divtmop1 = divtmop1 / 10;       divtmop2 = divtmop2 / 10;       divtmop3 = divtmop3 / 10;




# divergencias en la columna de resultados
resultados$`Div-KL`[32:62] = c(divShenoy17,divShenoy14,divShenoy13,divShenoy27,divShenoy24,
                              divShenoy23,divShenoy37,divShenoy34,divShenoy33,divShenoy47,
                              divShenoy44,divShenoy43,divShenoy42,divPB1,divPB2,divPB3,
                              divtmop1,divtmop2,divtmop3)




####################################################################

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (WEIBULL(1.5,1))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,0.75,length.out=100);              x2 = seq(0.75,1.5,length.out=100);
x3 = seq(1.5,2.25,length.out=100);            x4 = seq(2.25,3,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xWei3,POTENTIAL_TYPE="MOP",evalRange=c(0,3),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xWei3,maxDegree=7,0,3,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fWeibull151,c(0,3),col='#0072B2',type='l',xlim = c(0,3),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')





####################################################################

# MÉTODO SHENOY (WEIBULL(2,2))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fWeibull22,7,0,5,1);
fitShenoy14 = ShenoyFit(fWeibull22,4,0,5,1);
fitShenoy13 = ShenoyFit(fWeibull22,3,0,5,1);
fitShenoy27 = ShenoyFit(fWeibull22,7,0,5,2);
fitShenoy24 = ShenoyFit(fWeibull22,4,0,5,2);
fitShenoy23 = ShenoyFit(fWeibull22,3,0,5,2);
fitShenoy37 = ShenoyFit(fWeibull22,7,0,5,3);
fitShenoy34 = ShenoyFit(fWeibull22,4,0,5,3);
fitShenoy33 = ShenoyFit(fWeibull22,3,0,5,3);
fitShenoy47 = ShenoyFit(fWeibull22,7,0,5,4);
fitShenoy44 = ShenoyFit(fWeibull22,4,0,5,4);
fitShenoy43 = ShenoyFit(fWeibull22,3,0,5,4);
fitShenoy42 = ShenoyFit(fWeibull22,2,0,5,4);



divShenoy17 = DivergenceKL2(fWeibull22,fitShenoy17[[1]],0,5);
divShenoy14 = DivergenceKL2(fWeibull22,fitShenoy14[[1]],0,5);
divShenoy13 = DivergenceKL2(fWeibull22,fitShenoy13[[1]],0,5);

sp2 = seq(0,5,length.out=3);
divShenoy27 = DivergenceKL2(fWeibull22,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fWeibull22,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fWeibull22,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fWeibull22,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fWeibull22,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fWeibull22,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,5,length.out=4);
divShenoy37 = DivergenceKL2(fWeibull22,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fWeibull22,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fWeibull22,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fWeibull22,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fWeibull22,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fWeibull22,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fWeibull22,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fWeibull22,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fWeibull22,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,5,length.out=5);
divShenoy47 = DivergenceKL2(fWeibull22,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull22,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull22,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull22,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fWeibull22,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull22,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull22,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull22,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fWeibull22,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull22,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull22,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull22,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fWeibull22,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull22,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull22,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull22,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ  (WEIBULL(2,2))
# tiene un intervalo, grados 7, 4, 3 y tamaño de muestra: 50, 100, 1000


# ajustes para distintos tamaños de muestra. Repetir cada uno con 10 muestras
divPB1 = rep(0,3);        divPB2 = rep(0,3);        divPB3 = rep(0,3);
p = c(8,5,4)

# 50
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei1 = rweibull(50,2,2);           xWei1 = xWei1[xWei1<=5];
    
    fitPB1 = univMoTBF(xWei1,POTENTIAL_TYPE="MOP",evalRange=c(0,5),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fWeibull22,fitPBpoly1,0,5);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei2 = rweibull(100,2,2);           xWei2 = xWei2[xWei2<=5];
    
    fitPB2 = univMoTBF(xWei2,POTENTIAL_TYPE="MOP",evalRange=c(0,5),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fWeibull22,fitPBpoly2,0,5);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei3 = rweibull(1000,2,2);           xWei3 = xWei3[xWei3<=5];
    
    fitPB3 = univMoTBF(xWei3,POTENTIAL_TYPE="MOP",evalRange=c(0,5),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fWeibull22,fitPBpoly3,0,5);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (WEIBULL(2,2))
# tiene un intervalo, grados 7, 4, 3 y tamaño de muestra: 50, 100, 1000


# ajustes para distintos tamaños de muestra. Repetir cada uno con 10 muestras
divtmop1 = rep(0,3);        divtmop2 = rep(0,3);        divtmop3 = rep(0,3);
grados = c(7,4,3)

# 50
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei1 = rweibull(50,2,2);           xWei1 = xWei1[xWei1<=5];
    
    fittmop1 = PolynomialFit(xWei1,maxDegree=grados[j],0,5,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fWeibull22,fittmop1,0,5);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei2 = rweibull(100,2,2);           xWei2 = xWei2[xWei2<=5];
    
    fittmop2 = PolynomialFit(xWei2,maxDegree=grados[j],0,5,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fWeibull22,fittmop2,0,5);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei3 = rweibull(1000,2,2);           xWei3 = xWei3[xWei3<=5];
    
    fittmop3 = PolynomialFit(xWei3,maxDegree=grados[j],0,5,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fWeibull22,fittmop3,0,5);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divtmop1 = divtmop1 / 10;       divtmop2 = divtmop2 / 10;       divtmop3 = divtmop3 / 10;




# divergencias en la columna de resultados
resultados$`Div-KL`[63:93] = c(divShenoy17,divShenoy14,divShenoy13,divShenoy27,divShenoy24,
                              divShenoy23,divShenoy37,divShenoy34,divShenoy33,divShenoy47,
                              divShenoy44,divShenoy43,divShenoy42,divPB1,divPB2,divPB3,
                              divtmop1,divtmop2,divtmop3)




####################################################################

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (WEIBULL(2,2))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,1.25,length.out=100);              x2 = seq(1.25,2.5,length.out=100);
x3 = seq(2.5,3.75,length.out=100);            x4 = seq(3.75,5,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xWei3,POTENTIAL_TYPE="MOP",evalRange=c(0,5),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xWei3,maxDegree=7,0,5,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fWeibull22,c(0,5),col='#0072B2',type='l',xlim = c(0,5),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')






####################################################################

# MÉTODO SHENOY (WEIBULL(6,3))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fWeibull63,7,0,4,1);
fitShenoy14 = ShenoyFit(fWeibull63,4,0,4,1);
fitShenoy13 = ShenoyFit(fWeibull63,3,0,4,1);
fitShenoy27 = ShenoyFit(fWeibull63,7,0,4,2);
fitShenoy24 = ShenoyFit(fWeibull63,4,0,4,2);
fitShenoy23 = ShenoyFit(fWeibull63,3,0,4,2);
fitShenoy37 = ShenoyFit(fWeibull63,7,0,4,3);
fitShenoy34 = ShenoyFit(fWeibull63,4,0,4,3);
fitShenoy33 = ShenoyFit(fWeibull63,3,0,4,3);
fitShenoy47 = ShenoyFit(fWeibull63,7,0,4,4);
fitShenoy44 = ShenoyFit(fWeibull63,4,0,4,4);
fitShenoy43 = ShenoyFit(fWeibull63,3,0,4,4);
fitShenoy42 = ShenoyFit(fWeibull63,2,0,4,4);



divShenoy17 = DivergenceKL2(fWeibull63,fitShenoy17[[1]],0,4);
divShenoy14 = DivergenceKL2(fWeibull63,fitShenoy14[[1]],0,4);
divShenoy13 = DivergenceKL2(fWeibull63,fitShenoy13[[1]],0,4);

sp2 = seq(0,4,length.out=3);
divShenoy27 = DivergenceKL2(fWeibull63,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fWeibull63,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fWeibull63,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fWeibull63,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fWeibull63,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fWeibull63,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,4,length.out=4);
divShenoy37 = DivergenceKL2(fWeibull63,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fWeibull63,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fWeibull63,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fWeibull63,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fWeibull63,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fWeibull63,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fWeibull63,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fWeibull63,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fWeibull63,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,4,length.out=5);
divShenoy47 = DivergenceKL2(fWeibull63,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull63,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull63,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull63,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fWeibull63,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull63,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull63,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull63,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fWeibull63,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull63,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull63,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull63,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fWeibull63,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull63,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull63,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull63,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ  (WEIBULL(6,3))
# tiene un intervalo, grados 7, 4, 3 y tamaño de muestra: 50, 100, 1000


# ajustes para distintos tamaños de muestra. Repetir cada uno con 10 muestras
divPB1 = rep(0,3);        divPB2 = rep(0,3);        divPB3 = rep(0,3);
p = c(8,5,4)

# 50
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei1 = rweibull(50,6,3);           xWei1 = xWei1[xWei1<=4];
    
    fitPB1 = univMoTBF(xWei1,POTENTIAL_TYPE="MOP",evalRange=c(0,4),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fWeibull63,fitPBpoly1,0,4);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei2 = rweibull(100,6,3);           xWei2 = xWei2[xWei2<=4];
    
    fitPB2 = univMoTBF(xWei2,POTENTIAL_TYPE="MOP",evalRange=c(0,4),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fWeibull63,fitPBpoly2,0,4);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei3 = rweibull(1000,6,3);           xWei3 = xWei3[xWei3<=4];
    
    fitPB3 = univMoTBF(xWei3,POTENTIAL_TYPE="MOP",evalRange=c(0,4),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fWeibull63,fitPBpoly3,0,4);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (WEIBULL(6,3))
# tiene un intervalo, grados 7, 4, 3 y tamaño de muestra: 50, 100, 1000


# ajustes para distintos tamaños de muestra. Repetir cada uno con 10 muestras
divtmop1 = rep(0,3);        divtmop2 = rep(0,3);        divtmop3 = rep(0,3);
grados = c(7,4,3)

# 50
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei1 = rweibull(50,6,3);           xWei1 = xWei1[xWei1<=4];
    
    fittmop1 = PolynomialFit(xWei1,maxDegree=grados[j],0,4,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fWeibull63,fittmop1,0,4);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei2 = rweibull(100,6,3);           xWei2 = xWei2[xWei2<=4];
    
    fittmop2 = PolynomialFit(xWei2,maxDegree=grados[j],0,4,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fWeibull63,fittmop2,0,4);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei3 = rweibull(1000,6,3);           xWei3 = xWei3[xWei3<=4];
    
    fittmop3 = PolynomialFit(xWei3,maxDegree=grados[j],0,4,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fWeibull63,fittmop3,0,4);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divtmop1 = divtmop1 / 10;       divtmop2 = divtmop2 / 10;       divtmop3 = divtmop3 / 10;




# divergencias en la columna de resultados
resultados$`Div-KL`[94:124] = c(divShenoy17,divShenoy14,divShenoy13,divShenoy27,divShenoy24,
                              divShenoy23,divShenoy37,divShenoy34,divShenoy33,divShenoy47,
                              divShenoy44,divShenoy43,divShenoy42,divPB1,divPB2,divPB3,
                              divtmop1,divtmop2,divtmop3)




####################################################################

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (WEIBULL(6,3))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,1,length.out=100);              x2 = seq(1,2,length.out=100);
x3 = seq(2,3,length.out=100);              x4 = seq(3,4,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xWei3,POTENTIAL_TYPE="MOP",evalRange=c(0,4),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xWei3,maxDegree=7,0,4,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];
ytmop[length(ytmop)] = ytmop[length(ytmop)-1];


## GUARDAR EN PDF
pdf("GraficaWeibull63.pdf",width=6,height=5)

# distribución original
plot(fWeibull63,c(0,4),col='#0072B2',type='l',xlim=c(0,4),lwd=4,lty=2,ylab=NA,xlab=NA)
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#a23d75',type='l',lwd=2)
legend(0,0.7, legend=c("Theoretical", "Shenoy","MoTBF","tMoP"),seg.len=2,lwd=2,
       col=c('#0072B2','#D55E00','#009E73','#a23d75'),cex=0.8,bty='n',lty=c(2,1,1,1))

dev.off()





####################################################################

# MÉTODO SHENOY (WEIBULL(8,4))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fWeibull84,7,0,5,1);
fitShenoy14 = ShenoyFit(fWeibull84,4,0,5,1);
fitShenoy13 = ShenoyFit(fWeibull84,3,0,5,1);
fitShenoy27 = ShenoyFit(fWeibull84,7,0,5,2);
fitShenoy24 = ShenoyFit(fWeibull84,4,0,5,2);
fitShenoy23 = ShenoyFit(fWeibull84,3,0,5,2);
fitShenoy37 = ShenoyFit(fWeibull84,7,0,5,3);
fitShenoy34 = ShenoyFit(fWeibull84,4,0,5,3);
fitShenoy33 = ShenoyFit(fWeibull84,3,0,5,3);
fitShenoy47 = ShenoyFit(fWeibull84,7,0,5,4);
fitShenoy44 = ShenoyFit(fWeibull84,4,0,5,4);
fitShenoy43 = ShenoyFit(fWeibull84,3,0,5,4);
fitShenoy42 = ShenoyFit(fWeibull84,2,0,5,4);



divShenoy17 = DivergenceKL2(fWeibull84,fitShenoy17[[1]],0,5);
divShenoy14 = DivergenceKL2(fWeibull84,fitShenoy14[[1]],0,5);
divShenoy13 = DivergenceKL2(fWeibull84,fitShenoy13[[1]],0,5);

sp2 = seq(0,5,length.out=3);
divShenoy27 = DivergenceKL2(fWeibull84,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fWeibull84,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fWeibull84,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fWeibull84,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fWeibull84,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fWeibull84,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,5,length.out=4);
divShenoy37 = DivergenceKL2(fWeibull84,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fWeibull84,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fWeibull84,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fWeibull84,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fWeibull84,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fWeibull84,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fWeibull84,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fWeibull84,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fWeibull84,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,5,length.out=5);
divShenoy47 = DivergenceKL2(fWeibull84,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull84,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull84,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull84,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fWeibull84,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull84,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull84,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull84,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fWeibull84,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull84,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull84,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull84,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fWeibull84,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fWeibull84,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fWeibull84,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fWeibull84,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ  (WEIBULL(8,4))
# tiene un intervalo, grados 7, 4, 3 y tamaño de muestra: 50, 100, 1000


# ajustes para distintos tamaños de muestra. Repetir cada uno con 10 muestras
divPB1 = rep(0,3);        divPB2 = rep(0,3);        divPB3 = rep(0,3);
p = c(8,5,4)

# 50
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei1 = rweibull(50,8,4);           xWei1 = xWei1[xWei1<=5];
    
    fitPB1 = univMoTBF(xWei1,POTENTIAL_TYPE="MOP",evalRange=c(0,5),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fWeibull84,fitPBpoly1,0,5);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei2 = rweibull(100,8,4);           xWei2 = xWei2[xWei2<=5];
    
    fitPB2 = univMoTBF(xWei2,POTENTIAL_TYPE="MOP",evalRange=c(0,5),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fWeibull84,fitPBpoly2,0,5);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei3 = rweibull(1000,8,4);           xWei3 = xWei3[xWei3<=5];
    
    fitPB3 = univMoTBF(xWei3,POTENTIAL_TYPE="MOP",evalRange=c(0,5),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fWeibull84,fitPBpoly3,0,5);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (WEIBULL(8,4))
# tiene un intervalo, grados 7, 4, 3 y tamaño de muestra: 50, 100, 1000


# ajustes para distintos tamaños de muestra. Repetir cada uno con 10 muestras
divtmop1 = rep(0,3);        divtmop2 = rep(0,3);        divtmop3 = rep(0,3);
grados = c(7,4,3)

# 50
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei1 = rweibull(50,8,4);           xWei1 = xWei1[xWei1<=5];
    
    fittmop1 = PolynomialFit(xWei1,maxDegree=grados[j],0,5,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fWeibull84,fittmop1,0,5);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei2 = rweibull(100,8,4);           xWei2 = xWei2[xWei2<=5];
    
    fittmop2 = PolynomialFit(xWei2,maxDegree=grados[j],0,5,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fWeibull84,fittmop2,0,5);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xWei3 = rweibull(1000,8,4);           xWei3 = xWei3[xWei3<=5];
    
    fittmop3 = PolynomialFit(xWei3,maxDegree=grados[j],0,5,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fWeibull84,fittmop3,0,5);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divtmop1 = divtmop1 / 10;       divtmop2 = divtmop2 / 10;       divtmop3 = divtmop3 / 10;




# divergencias en la columna de resultados
resultados$`Div-KL`[125:155] =c(divShenoy17,divShenoy14,divShenoy13,divShenoy27,divShenoy24,
                              divShenoy23,divShenoy37,divShenoy34,divShenoy33,divShenoy47,
                              divShenoy44,divShenoy43,divShenoy42,divPB1,divPB2,divPB3,
                              divtmop1,divtmop2,divtmop3)


resultadosWeibull = resultados;

save(resultadosWeibull,file='ResultadosWeibull.RData');






####################################################################

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (WEIBULL(8,4))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,1.25,length.out=100);              x2 = seq(1.25,2.5,length.out=100);
x3 = seq(2.5,3.75,length.out=100);            x4 = seq(3.75,5,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xWei3,POTENTIAL_TYPE="MOP",evalRange=c(0,5),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xWei3,maxDegree=7,0,5,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];
ytmop[length(ytmop)] = ytmop[length(ytmop)-1];

# distribución original
plot(fWeibull84,c(0,5),col='#0072B2',type='l',xlim = c(0,5),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')











