
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

# DISTRIBUCIÓN LOGNORMAL de parámetros (0,0.25)  en el intervalo [0,4]
fLN025 = function(x)   1/(x*sqrt(2*pi*0.25)) * exp(-(log(x)-0)^2/(2*0.25))
iLN025 = pracma::integral(fLN025,0,4);
fLogNorm025 = function(x)   (1/(x*sqrt(2*pi*0.25)) * exp(-(log(x)-0)^2/(2*0.25))) / iLN025

# DISTRIBUCIÓN LOGNORMAL de parámetros (0,0.5) en el intervalo [0,6]
fLN05 = function(x)   1/(x*sqrt(2*pi*0.5)) * exp(-(log(x)-0)^2/(2*0.5))
iLN05 = pracma::integral(fLN05,0,6);
fLogNorm05 = function(x)   (1/(x*sqrt(2*pi*0.5)) * exp(-(log(x)-0)^2/(2*0.5))) / iLN05

# DISTRIBUCIÓN LOGNORMAL de parámetros (0,1)  en el intervalo [0,11]
fLN01 = function(x)   1/(x*sqrt(2*pi*1)) * exp(-(log(x)-0)^2/(2*1))
iLN01 = pracma::integral(fLN01,0,11);
fLogNorm01 = function(x)   (1/(x*sqrt(2*pi*1)) * exp(-(log(x)-0)^2/(2*1))) / iLN01



####################################################################

# Shenoy
# - Número trozos: 1,2,3,4
# - Grado: 1 trozo (7,4,3), 2 trozos (7,4,3), 3 trozos (7,4,3), 4 trozos (7,4,3,2)

# tmop y Pérez-Bernabé
# - Grados: 7, 4 y 3

# Otros métodos:
# - Tamaño muestra: 50,100,1000


# data frame para los resultados
distribucion1 = rep('lognorm(0,0.25)',31);
intervalo1 = rep('[0,4]',31)
distribucion2 = rep('lognorm(0,0.5)',31);
intervalo2 = rep('[0,6]',31)
distribucion3 = rep('lognorm(0,0.1)',31);
intervalo3 = rep('[0,11]',31)
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

resultados = rbind(r1,r2,r3)





####################################################################

# MÉTODO SHENOY  (LOGNORM(0,0.25))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fLogNorm025,7,0,4,1);
fitShenoy14 = ShenoyFit(fLogNorm025,4,0,4,1);
fitShenoy13 = ShenoyFit(fLogNorm025,3,0,4,1);
fitShenoy27 = ShenoyFit(fLogNorm025,7,0,4,2);
fitShenoy24 = ShenoyFit(fLogNorm025,4,0,4,2);
fitShenoy23 = ShenoyFit(fLogNorm025,3,0,4,2);
fitShenoy37 = ShenoyFit(fLogNorm025,7,0,4,3);
fitShenoy34 = ShenoyFit(fLogNorm025,4,0,4,3);
fitShenoy33 = ShenoyFit(fLogNorm025,3,0,4,3);
fitShenoy47 = ShenoyFit(fLogNorm025,7,0,4,4);
fitShenoy44 = ShenoyFit(fLogNorm025,4,0,4,4);
fitShenoy43 = ShenoyFit(fLogNorm025,3,0,4,4);
fitShenoy42 = ShenoyFit(fLogNorm025,2,0,4,4);



divShenoy17 = DivergenceKL2(fLogNorm025,fitShenoy17[[1]],0,4);
divShenoy14 = DivergenceKL2(fLogNorm025,fitShenoy14[[1]],0,4);
divShenoy13 = DivergenceKL2(fLogNorm025,fitShenoy13[[1]],0,4);

sp2 = seq(0,4,length.out=3);
divShenoy27 = DivergenceKL2(fLogNorm025,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fLogNorm025,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fLogNorm025,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fLogNorm025,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fLogNorm025,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fLogNorm025,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,4,length.out=4);
divShenoy37 = DivergenceKL2(fLogNorm025,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fLogNorm025,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fLogNorm025,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fLogNorm025,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fLogNorm025,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fLogNorm025,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fLogNorm025,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fLogNorm025,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fLogNorm025,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,4,length.out=5);
divShenoy47 = DivergenceKL2(fLogNorm025,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fLogNorm025,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fLogNorm025,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fLogNorm025,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fLogNorm025,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fLogNorm025,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fLogNorm025,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fLogNorm025,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fLogNorm025,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fLogNorm025,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fLogNorm025,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fLogNorm025,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fLogNorm025,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fLogNorm025,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fLogNorm025,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fLogNorm025,fitShenoy42[[4]],sp4[4],sp4[5]);





####################################################################

# MÉTODO PÉREZ-BERNABÉ  (LOGNORM(0,0.25))
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
    xLogNorm1 = rlnorm(50,0,0.25);            xLogNorm1 = xLogNorm1[xLogNorm1<=4];
    
    fitPB1 = univMoTBF(xLogNorm1,POTENTIAL_TYPE="MOP",evalRange=c(0,4),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fLogNorm025,fitPBpoly1,0,4);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xLogNorm2 = rlnorm(100,0,0.25);            xLogNorm2 = xLogNorm2[xLogNorm2<=4];
    
    fitPB2 = univMoTBF(xLogNorm2,POTENTIAL_TYPE="MOP",evalRange=c(0,4),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fLogNorm025,fitPBpoly2,0,4);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xLogNorm3 = rlnorm(1000,0,0.25);            xLogNorm3 = xLogNorm3[xLogNorm3<=4];
    
    fitPB3 = univMoTBF(xLogNorm3,POTENTIAL_TYPE="MOP",evalRange=c(0,4),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fLogNorm025,fitPBpoly3,0,4);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (LOGNORM(0,0.25))
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
    xLogNorm1 = rlnorm(50,0,0.25);            xLogNorm1 = xLogNorm1[xLogNorm1<=4];
    
    fittmop1 = PolynomialFit(xLogNorm1,maxDegree=grados[j],0,4,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fLogNorm025,fittmop1,0,4);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xLogNorm2 = rlnorm(100,0,0.25);            xLogNorm2 = xLogNorm2[xLogNorm2<=4];
    
    fittmop2 = PolynomialFit(xLogNorm2,maxDegree=grados[j],0,4,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fLogNorm025,fittmop2,0,4);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xLogNorm3 = rlnorm(1000,0,0.25);            xLogNorm3 = xLogNorm3[xLogNorm3<=4];
    
    fittmop3 = PolynomialFit(xLogNorm3,maxDegree=grados[j],0,4,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fLogNorm025,fittmop3,0,4);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO PARA LOGNORMAL (0,0.25)
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,1,length.out=100);           x2 = seq(1,2,length.out=100);
x3 = seq(2,3,length.out=100);           x4 = seq(3,4,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xLogNorm3,POTENTIAL_TYPE="MOP",evalRange=c(0,4),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xLogNorm3,maxDegree=7,0,4,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fLogNorm025,c(0,4),col='#0072B2',type='l',xlim = c(0,4),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')







####################################################################

# MÉTODO SHENOY  (LOGNORM(0,0.5))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fLogNorm05,7,0,6,1);
fitShenoy14 = ShenoyFit(fLogNorm05,4,0,6,1);
fitShenoy13 = ShenoyFit(fLogNorm05,3,0,6,1);
fitShenoy27 = ShenoyFit(fLogNorm05,7,0,6,2);
fitShenoy24 = ShenoyFit(fLogNorm05,4,0,6,2);
fitShenoy23 = ShenoyFit(fLogNorm05,3,0,6,2);
fitShenoy37 = ShenoyFit(fLogNorm05,7,0,6,3);
fitShenoy34 = ShenoyFit(fLogNorm05,4,0,6,3);
fitShenoy33 = ShenoyFit(fLogNorm05,3,0,6,3);
fitShenoy47 = ShenoyFit(fLogNorm05,7,0,6,4);
fitShenoy44 = ShenoyFit(fLogNorm05,4,0,6,4);
fitShenoy43 = ShenoyFit(fLogNorm05,3,0,6,4);
fitShenoy42 = ShenoyFit(fLogNorm05,2,0,6,4);



divShenoy17 = DivergenceKL2(fLogNorm05,fitShenoy17[[1]],0,6);
divShenoy14 = DivergenceKL2(fLogNorm05,fitShenoy14[[1]],0,6);
divShenoy13 = DivergenceKL2(fLogNorm05,fitShenoy13[[1]],0,6);

sp2 = seq(0,6,length.out=3);
divShenoy27 = DivergenceKL2(fLogNorm05,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fLogNorm05,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fLogNorm05,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fLogNorm05,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fLogNorm05,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fLogNorm05,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,6,length.out=4);
divShenoy37 = DivergenceKL2(fLogNorm05,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fLogNorm05,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fLogNorm05,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fLogNorm05,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fLogNorm05,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fLogNorm05,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fLogNorm05,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fLogNorm05,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fLogNorm05,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,6,length.out=5);
divShenoy47 = DivergenceKL2(fLogNorm05,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fLogNorm05,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fLogNorm05,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fLogNorm05,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fLogNorm05,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fLogNorm05,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fLogNorm05,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fLogNorm05,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fLogNorm05,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fLogNorm05,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fLogNorm05,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fLogNorm05,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fLogNorm05,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fLogNorm05,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fLogNorm05,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fLogNorm05,fitShenoy42[[4]],sp4[4],sp4[5]);





####################################################################

# MÉTODO PÉREZ-BERNABÉ  (LOGNORM(0,0.5))
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
    xLogNorm1 = rlnorm(50,0,0.5);            xLogNorm1 = xLogNorm1[xLogNorm1<=6];
    
    fitPB1 = univMoTBF(xLogNorm1,POTENTIAL_TYPE="MOP",evalRange=c(0,6),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fLogNorm05,fitPBpoly1,0,6);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xLogNorm2 = rlnorm(100,0,0.5);            xLogNorm2 = xLogNorm2[xLogNorm2<=6];
    
    fitPB2 = univMoTBF(xLogNorm2,POTENTIAL_TYPE="MOP",evalRange=c(0,6),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fLogNorm05,fitPBpoly2,0,6);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xLogNorm3 = rlnorm(1000,0,0.5);            xLogNorm3 = xLogNorm3[xLogNorm3<=6];
    
    fitPB3 = univMoTBF(xLogNorm3,POTENTIAL_TYPE="MOP",evalRange=c(0,6),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fLogNorm05,fitPBpoly3,0,6);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (LOGNORM(0,0.5))
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
    xLogNorm1 = rlnorm(50,0,0.5);            xLogNorm1 = xLogNorm1[xLogNorm1<=6];
    
    fittmop1 = PolynomialFit(xLogNorm1,maxDegree=grados[j],0,6,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fLogNorm05,fittmop1,0,6);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xLogNorm2 = rlnorm(100,0,0.5);            xLogNorm2 = xLogNorm2[xLogNorm2<=6];
    
    fittmop2 = PolynomialFit(xLogNorm2,maxDegree=grados[j],0,6,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fLogNorm05,fittmop2,0,6);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xLogNorm3 = rlnorm(1000,0,0.5);            xLogNorm3 = xLogNorm3[xLogNorm3<=6];
    
    fittmop3 = PolynomialFit(xLogNorm3,maxDegree=grados[j],0,6,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fLogNorm05,fittmop3,0,6);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO PARA LOGNORMAL (0,0.5)
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,1.5,length.out=100);           x2 = seq(1.5,3,length.out=100);
x3 = seq(3,4.5,length.out=100);           x4 = seq(4.5,6,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xLogNorm3,POTENTIAL_TYPE="MOP",evalRange=c(0,6),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xLogNorm3,maxDegree=7,0,6,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];


## GUARDAR EN PDF
pdf("GraficaLogNorm005.pdf",width=6,height=5)

# distribución original
plot(fLogNorm05,c(0,6),col='#0072B2',type='l',xlim=c(0,6),ylim=c(0,0.8),lwd=4,lty=2,
     ylab=NA,xlab=NA)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,0.75, legend=c("Original", "Shenoy","MoTBF","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.8,bty='n')

dev.off()









####################################################################

# MÉTODO SHENOY  (LOGNORM(0,1))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fLogNorm01,7,0,11,1);
fitShenoy14 = ShenoyFit(fLogNorm01,4,0,11,1);
fitShenoy13 = ShenoyFit(fLogNorm01,3,0,11,1);
fitShenoy27 = ShenoyFit(fLogNorm01,7,0,11,2);
fitShenoy24 = ShenoyFit(fLogNorm01,4,0,11,2);
fitShenoy23 = ShenoyFit(fLogNorm01,3,0,11,2);
fitShenoy37 = ShenoyFit(fLogNorm01,7,0,11,3);
fitShenoy34 = ShenoyFit(fLogNorm01,4,0,11,3);
fitShenoy33 = ShenoyFit(fLogNorm01,3,0,11,3);
fitShenoy47 = ShenoyFit(fLogNorm01,7,0,11,4);
fitShenoy44 = ShenoyFit(fLogNorm01,4,0,11,4);
fitShenoy43 = ShenoyFit(fLogNorm01,3,0,11,4);
fitShenoy42 = ShenoyFit(fLogNorm01,2,0,11,4);



divShenoy17 = DivergenceKL2(fLogNorm01,fitShenoy17[[1]],0,11);
divShenoy14 = DivergenceKL2(fLogNorm01,fitShenoy14[[1]],0,11);
divShenoy13 = DivergenceKL2(fLogNorm01,fitShenoy13[[1]],0,11);

sp2 = seq(0,11,length.out=3);
divShenoy27 = DivergenceKL2(fLogNorm01,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fLogNorm01,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fLogNorm01,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fLogNorm01,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fLogNorm01,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fLogNorm01,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,11,length.out=4);
divShenoy37 = DivergenceKL2(fLogNorm01,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fLogNorm01,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fLogNorm01,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fLogNorm01,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fLogNorm01,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fLogNorm01,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fLogNorm01,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fLogNorm01,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fLogNorm01,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,11,length.out=5);
divShenoy47 = DivergenceKL2(fLogNorm01,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fLogNorm01,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fLogNorm01,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fLogNorm01,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fLogNorm01,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fLogNorm01,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fLogNorm01,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fLogNorm01,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fLogNorm01,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fLogNorm01,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fLogNorm01,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fLogNorm01,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fLogNorm01,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fLogNorm01,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fLogNorm01,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fLogNorm01,fitShenoy42[[4]],sp4[4],sp4[5]);





####################################################################

# MÉTODO PÉREZ-BERNABÉ  (LOGNORM(0,1))
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
    xLogNorm1 = rlnorm(50,0,1);            xLogNorm1 = xLogNorm1[xLogNorm1<=11];
    
    fitPB1 = univMoTBF(xLogNorm1,POTENTIAL_TYPE="MOP",evalRange=c(0,11),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fLogNorm01,fitPBpoly1,0,11);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xLogNorm2 = rlnorm(100,0,1);            xLogNorm2 = xLogNorm2[xLogNorm2<=11];
    
    fitPB2 = univMoTBF(xLogNorm2,POTENTIAL_TYPE="MOP",evalRange=c(0,11),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fLogNorm01,fitPBpoly2,0,11);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xLogNorm3 = rlnorm(1000,0,1);            xLogNorm3 = xLogNorm3[xLogNorm3<=11];
    
    fitPB3 = univMoTBF(xLogNorm3,POTENTIAL_TYPE="MOP",evalRange=c(0,11),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fLogNorm01,fitPBpoly3,0,11);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (LOGNORM(0,1))
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
    xLogNorm1 = rlnorm(50,0,1);            xLogNorm1 = xLogNorm1[xLogNorm1<=11];
    
    fittmop1 = PolynomialFit(xLogNorm1,maxDegree=grados[j],0,11,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fLogNorm01,fittmop1,0,11);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xLogNorm2 = rlnorm(100,0,1);            xLogNorm2 = xLogNorm2[xLogNorm2<=11];
    
    fittmop2 = PolynomialFit(xLogNorm2,maxDegree=grados[j],0,11,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fLogNorm01,fittmop2,0,11);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xLogNorm3 = rlnorm(1000,0,1);            xLogNorm3 = xLogNorm3[xLogNorm3<=11];
    
    fittmop3 = PolynomialFit(xLogNorm3,maxDegree=grados[j],0,11,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fLogNorm01,fittmop3,0,11);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divtmop1 = divtmop1 / 10;       divtmop2 = divtmop2 / 10;       divtmop3 = divtmop3 / 10;




# divergencias en la columna de resultados
resultados$`Div-KL`[63:93] = c(divShenoy17,divShenoy14,divShenoy13,divShenoy27,divShenoy24,
                              divShenoy23,divShenoy37,divShenoy34,divShenoy33,divShenoy47,
                              divShenoy44,divShenoy43,divShenoy42,divPB1,divPB2,divPB3,
                              divtmop1,divtmop2,divtmop3)



resultadosLogNorm = resultados;

save(resultadosLogNorm,file='ResultadosLogNorm.RData');






####################################################################

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO PARA LOGNORMAL (0,1)
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,2.75,length.out=100);            x2 = seq(2.75,5.5,length.out=100);
x3 = seq(5.5,8.25,length.out=100);          x4 = seq(8.25,11,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xLogNorm3,POTENTIAL_TYPE="MOP",evalRange=c(0,11),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xLogNorm3,maxDegree=7,0,11,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fLogNorm01,c(0,11),col='#0072B2',type='l',xlim = c(0,11),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')






