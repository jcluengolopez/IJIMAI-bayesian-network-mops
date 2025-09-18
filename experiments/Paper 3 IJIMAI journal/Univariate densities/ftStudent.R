
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

# DISTRIBUCIÓN T-STUDENT con v=1 en el intervalo [-30,30]
ftS1 = function(x)   gamma((1+1)/2) / (sqrt(1*pi) * gamma(1/2)) * (1+x^2/1)^(-(1+1)/2);
itS1 = pracma::integral(ftS1,-30,30);
ftStudent1 = function(x)   (gamma((1+1)/2)/(sqrt(1*pi)*gamma(1/2))*(1+x^2/1)^(-(1+1)/2))/itS1

# DISTRIBUCIÓN T-STUDENT con v=2 en el intervalo [-8,8]
ftS2 = function(x)   gamma((2+1)/2) / (sqrt(2*pi) * gamma(2/2)) * (1+x^2/2)^(-(2+1)/2);
itS2 = pracma::integral(ftS2,-8,8);
ftStudent2 = function(x)   (gamma((2+1)/2)/(sqrt(2*pi)*gamma(2/2))*(1+x^2/2)^(-(2+1)/2))/itS2

# DISTRIBUCIÓN T-STUDENT con v=4 en el intervalo [-5,5]
ftS4 = function(x)   gamma((4+1)/2) / (sqrt(4*pi) * gamma(4/2)) * (1+x^2/4)^(-(4+1)/2);
itS4 = pracma::integral(ftS4,-5,5);
ftStudent4 = function(x)   (gamma((4+1)/2)/(sqrt(4*pi)*gamma(4/2))*(1+x^2/4)^(-(4+1)/2))/itS4




####################################################################

# Shenoy
# - Número trozos: 1,2,3,4
# - Grado: 1 trozo (7,4,3), 2 trozos (7,4,3), 3 trozos (7,4,3), 4 trozos (7,4,3,2)

# tmop y Pérez-Bernabé
# - Grados: 7, 4 y 3

# Otros métodos:
# - Tamaño muestra: 50,100,1000



# data frame para los resultados
distribucion1 = rep('t-student(1)',31);
intervalo1 = rep('[-30,30]',31);
distribucion2 = rep('t-student(2)',31);
intervalo2 = rep('[-8,8]',31);
distribucion3 = rep('t-student(4)',31);
intervalo3 = rep('[-5,5]',31);
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

# MÉTODO SHENOY (T-STUDENT(1))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(ftStudent1,7,-30,30,1);
fitShenoy14 = ShenoyFit(ftStudent1,4,-30,30,1);
fitShenoy13 = ShenoyFit(ftStudent1,3,-30,30,1);
fitShenoy27 = ShenoyFit(ftStudent1,7,-30,30,2);
fitShenoy24 = ShenoyFit(ftStudent1,4,-30,30,2);
fitShenoy23 = ShenoyFit(ftStudent1,3,-30,30,2);
fitShenoy37 = ShenoyFit(ftStudent1,7,-30,30,3);
fitShenoy34 = ShenoyFit(ftStudent1,4,-30,30,3);
fitShenoy33 = ShenoyFit(ftStudent1,3,-30,30,3);
fitShenoy47 = ShenoyFit(ftStudent1,7,-30,30,4);
fitShenoy44 = ShenoyFit(ftStudent1,4,-30,30,4);
fitShenoy43 = ShenoyFit(ftStudent1,3,-30,30,4);
fitShenoy42 = ShenoyFit(ftStudent1,2,-30,30,4);



divShenoy17 = DivergenceKL2(ftStudent1,fitShenoy17[[1]],-30,30);
divShenoy14 = DivergenceKL2(ftStudent1,fitShenoy14[[1]],-30,30);
divShenoy13 = DivergenceKL2(ftStudent1,fitShenoy13[[1]],-30,30);

sp2 = seq(-30,30,length.out=3);
divShenoy27 = DivergenceKL2(ftStudent1,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(ftStudent1,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(ftStudent1,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(ftStudent1,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(ftStudent1,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(ftStudent1,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(-30,30,length.out=4);
divShenoy37 = DivergenceKL2(ftStudent1,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(ftStudent1,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(ftStudent1,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(ftStudent1,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(ftStudent1,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(ftStudent1,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(ftStudent1,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(ftStudent1,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(ftStudent1,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(-30,30,length.out=5);
divShenoy47 = DivergenceKL2(ftStudent1,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(ftStudent1,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(ftStudent1,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(ftStudent1,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(ftStudent1,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(ftStudent1,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(ftStudent1,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(ftStudent1,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(ftStudent1,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(ftStudent1,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(ftStudent1,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(ftStudent1,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(ftStudent1,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(ftStudent1,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(ftStudent1,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(ftStudent1,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ  (T-STUDENT(1))
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
    xtStu1 = rt(50,1);            xtStu1 = xtStu1[xtStu1>=-30 & xtStu1<=30];
    
    fitPB1 = univMoTBF(xtStu1,POTENTIAL_TYPE="MOP",evalRange=c(-30,30),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(ftStudent1,fitPBpoly1,-30,30);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xtStu2 = rt(100,1);            xtStu2 = xtStu2[xtStu2>=-30 & xtStu2<=30];
    
    fitPB2 = univMoTBF(xtStu2,POTENTIAL_TYPE="MOP",evalRange=c(-30,30),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(ftStudent1,fitPBpoly2,-30,30);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xtStu3 = rt(1000,1);            xtStu3 = xtStu3[xtStu3>=-30 & xtStu3<=30];
    
    fitPB3 = univMoTBF(xtStu3,POTENTIAL_TYPE="MOP",evalRange=c(-30,30),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(ftStudent1,fitPBpoly3,-30,30);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (T-STUDENT(1))
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
    xtStu1 = rt(50,1);            xtStu1 = xtStu1[xtStu1>=-30 & xtStu1<=30];
    
    fittmop1 = PolynomialFit(xtStu1,maxDegree=grados[j],-30,30,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(ftStudent1,fittmop1,-30,30);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xtStu2 = rt(100,1);            xtStu2 = xtStu2[xtStu2>=-30 & xtStu2<=30];
    
    fittmop2 = PolynomialFit(xtStu2,maxDegree=grados[j],-30,30,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(ftStudent1,fittmop2,-30,30);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xtStu3 = rt(1000,1);            xtStu3 = xtStu3[xtStu3>=-30 & xtStu3<=30];
    
    fittmop3 = PolynomialFit(xtStu3,maxDegree=grados[j],-30,30,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(ftStudent1,fittmop3,-30,30);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (T-STUDENT(1))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(-30,-15,length.out=100);           x2 = seq(-15,0,length.out=100);
x3 = seq(0,15,length.out=100);              x4 = seq(15,30,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xtStu3,POTENTIAL_TYPE="MOP",evalRange=c(-30,30),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xtStu3,maxDegree=7,-30,30,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(ftStudent1,c(-30,30),col='#0072B2',type='l',xlim = c(-30,30),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')







####################################################################

# MÉTODO SHENOY (T-STUDENT(2))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(ftStudent2,7,-8,8,1);
fitShenoy14 = ShenoyFit(ftStudent2,4,-8,8,1);
fitShenoy13 = ShenoyFit(ftStudent2,3,-8,8,1);
fitShenoy27 = ShenoyFit(ftStudent2,7,-8,8,2);
fitShenoy24 = ShenoyFit(ftStudent2,4,-8,8,2);
fitShenoy23 = ShenoyFit(ftStudent2,3,-8,8,2);
fitShenoy37 = ShenoyFit(ftStudent2,7,-8,8,3);
fitShenoy34 = ShenoyFit(ftStudent2,4,-8,8,3);
fitShenoy33 = ShenoyFit(ftStudent2,3,-8,8,3);
fitShenoy47 = ShenoyFit(ftStudent2,7,-8,8,4);
fitShenoy44 = ShenoyFit(ftStudent2,4,-8,8,4);
fitShenoy43 = ShenoyFit(ftStudent2,3,-8,8,4);
fitShenoy42 = ShenoyFit(ftStudent2,2,-8,8,4);



divShenoy17 = DivergenceKL2(ftStudent2,fitShenoy17[[1]],-8,8);
divShenoy14 = DivergenceKL2(ftStudent2,fitShenoy14[[1]],-8,8);
divShenoy13 = DivergenceKL2(ftStudent2,fitShenoy13[[1]],-8,8);

sp2 = seq(-8,8,length.out=3);
divShenoy27 = DivergenceKL2(ftStudent2,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(ftStudent2,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(ftStudent2,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(ftStudent2,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(ftStudent2,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(ftStudent2,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(-8,8,length.out=4);
divShenoy37 = DivergenceKL2(ftStudent2,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(ftStudent2,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(ftStudent2,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(ftStudent2,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(ftStudent2,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(ftStudent2,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(ftStudent2,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(ftStudent2,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(ftStudent2,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(-8,8,length.out=5);
divShenoy47 = DivergenceKL2(ftStudent2,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(ftStudent2,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(ftStudent2,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(ftStudent2,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(ftStudent2,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(ftStudent2,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(ftStudent2,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(ftStudent2,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(ftStudent2,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(ftStudent2,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(ftStudent2,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(ftStudent2,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(ftStudent2,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(ftStudent2,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(ftStudent2,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(ftStudent2,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ  (T-STUDENT(2))
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
    xtStu1 = rt(50,2);            xtStu1 = xtStu1[xtStu1>=-8 & xtStu1<=8];
    
    fitPB1 = univMoTBF(xtStu1,POTENTIAL_TYPE="MOP",evalRange=c(-8,8),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(ftStudent2,fitPBpoly1,-8,8);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xtStu2 = rt(100,2);            xtStu2 = xtStu2[xtStu2>=-8 & xtStu2<=8];
    
    fitPB2 = univMoTBF(xtStu2,POTENTIAL_TYPE="MOP",evalRange=c(-8,8),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(ftStudent2,fitPBpoly2,-8,8);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xtStu3 = rt(1000,2);            xtStu3 = xtStu3[xtStu3>=-8 & xtStu3<=8];
    
    fitPB3 = univMoTBF(xtStu3,POTENTIAL_TYPE="MOP",evalRange=c(-8,8),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(ftStudent2,fitPBpoly3,-8,8);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (T-STUDENT(2))
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
    xtStu1 = rt(50,2);            xtStu1 = xtStu1[xtStu1>=-8 & xtStu1<=8];
    
    fittmop1 = PolynomialFit(xtStu1,maxDegree=grados[j],-8,8,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(ftStudent2,fittmop1,-8,8);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xtStu2 = rt(100,2);            xtStu2 = xtStu2[xtStu2>=-8 & xtStu2<=8];
    
    fittmop2 = PolynomialFit(xtStu2,maxDegree=grados[j],-8,8,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(ftStudent2,fittmop2,-8,8);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xtStu3 = rt(1000,2);            xtStu3 = xtStu3[xtStu3>=-8 & xtStu3<=8];
    
    fittmop3 = PolynomialFit(xtStu3,maxDegree=grados[j],-8,8,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(ftStudent2,fittmop3,-8,8);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (T-STUDENT(2))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(-8,-4,length.out=100);            x2 = seq(-4,0,length.out=100);
x3 = seq(0,4,length.out=100);              x4 = seq(4,8,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xtStu3,POTENTIAL_TYPE="MOP",evalRange=c(-8,8),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xtStu3,maxDegree=7,-8,8,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(ftStudent2,c(-8,8),col='#0072B2',type='l',xlim = c(-8,8),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')









####################################################################

# MÉTODO SHENOY (T-STUDENT(4))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(ftStudent4,7,-5,5,1);
fitShenoy14 = ShenoyFit(ftStudent4,4,-5,5,1);
fitShenoy13 = ShenoyFit(ftStudent4,3,-5,5,1);
fitShenoy27 = ShenoyFit(ftStudent4,7,-5,5,2);
fitShenoy24 = ShenoyFit(ftStudent4,4,-5,5,2);
fitShenoy23 = ShenoyFit(ftStudent4,3,-5,5,2);
fitShenoy37 = ShenoyFit(ftStudent4,7,-5,5,3);
fitShenoy34 = ShenoyFit(ftStudent4,4,-5,5,3);
fitShenoy33 = ShenoyFit(ftStudent4,3,-5,5,3);
fitShenoy47 = ShenoyFit(ftStudent4,7,-5,5,4);
fitShenoy44 = ShenoyFit(ftStudent4,4,-5,5,4);
fitShenoy43 = ShenoyFit(ftStudent4,3,-5,5,4);
fitShenoy42 = ShenoyFit(ftStudent4,2,-5,5,4);



divShenoy17 = DivergenceKL2(ftStudent4,fitShenoy17[[1]],-5,5);
divShenoy14 = DivergenceKL2(ftStudent4,fitShenoy14[[1]],-5,5);
divShenoy13 = DivergenceKL2(ftStudent4,fitShenoy13[[1]],-5,5);

sp2 = seq(-5,5,length.out=3);
divShenoy27 = DivergenceKL2(ftStudent4,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(ftStudent4,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(ftStudent4,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(ftStudent4,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(ftStudent4,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(ftStudent4,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(-5,5,length.out=4);
divShenoy37 = DivergenceKL2(ftStudent4,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(ftStudent4,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(ftStudent4,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(ftStudent4,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(ftStudent4,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(ftStudent4,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(ftStudent4,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(ftStudent4,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(ftStudent4,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(-5,5,length.out=5);
divShenoy47 = DivergenceKL2(ftStudent4,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(ftStudent4,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(ftStudent4,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(ftStudent4,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(ftStudent4,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(ftStudent4,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(ftStudent4,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(ftStudent4,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(ftStudent4,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(ftStudent4,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(ftStudent4,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(ftStudent4,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(ftStudent4,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(ftStudent4,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(ftStudent4,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(ftStudent4,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ  (T-STUDENT(4))
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
    xtStu1 = rt(50,4);            xtStu1 = xtStu1[xtStu1>=-5 & xtStu1<=5];
    
    fitPB1 = univMoTBF(xtStu1,POTENTIAL_TYPE="MOP",evalRange=c(-5,5),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(ftStudent4,fitPBpoly1,-5,5);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xtStu2 = rt(100,4);            xtStu2 = xtStu2[xtStu2>=-5 & xtStu2<=5];
    
    fitPB2 = univMoTBF(xtStu2,POTENTIAL_TYPE="MOP",evalRange=c(-5,5),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(ftStudent4,fitPBpoly2,-5,5);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xtStu3 = rt(1000,4);            xtStu3 = xtStu3[xtStu3>=-5 & xtStu3<=5];
    
    fitPB3 = univMoTBF(xtStu3,POTENTIAL_TYPE="MOP",evalRange=c(-5,5),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(ftStudent4,fitPBpoly3,-5,5);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (T-STUDENT(4))
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
    xtStu1 = rt(50,4);            xtStu1 = xtStu1[xtStu1>=-5 & xtStu1<=5];
    
    fittmop1 = PolynomialFit(xtStu1,maxDegree=grados[j],-5,5,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(ftStudent4,fittmop1,-5,5);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xtStu2 = rt(100,4);            xtStu2 = xtStu2[xtStu2>=-5 & xtStu2<=5];
    
    fittmop2 = PolynomialFit(xtStu2,maxDegree=grados[j],-5,5,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(ftStudent4,fittmop2,-5,5);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xtStu3 = rt(1000,4);            xtStu3 = xtStu3[xtStu3>=-5 & xtStu3<=5];
    
    fittmop3 = PolynomialFit(xtStu3,maxDegree=grados[j],-5,5,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(ftStudent4,fittmop3,-5,5);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divtmop1 = divtmop1 / 10;       divtmop2 = divtmop2 / 10;       divtmop3 = divtmop3 / 10;




# divergencias en la columna de resultados
resultados$`Div-KL`[63:93] = c(divShenoy17,divShenoy14,divShenoy13,divShenoy27,divShenoy24,
                              divShenoy23,divShenoy37,divShenoy34,divShenoy33,divShenoy47,
                              divShenoy44,divShenoy43,divShenoy42,divPB1,divPB2,divPB3,
                              divtmop1,divtmop2,divtmop3)



resultadostStudent = resultados;

save(resultadostStudent,file='ResultadostStudent.RData');







####################################################################

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (T-STUDENT(4))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(-5,-2.5,length.out=100);            x2 = seq(-2.5,0,length.out=100);
x3 = seq(0,2.5,length.out=100);              x4 = seq(2.5,5,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xtStu3,POTENTIAL_TYPE="MOP",evalRange=c(-5,5),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xtStu3,maxDegree=7,-5,5,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];


## GUARDAR EN PDF
pdf("GraficatStudent4.pdf",width=6,height=5)

# distribución original
plot(ftStudent4,c(-5,5),col='#0072B2',type='l',xlim = c(-5,5),lwd=4,lty=2,ylab=NA,xlab=NA)
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(3,0.35, legend=c("Original", "Shenoy","MoTBF","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.8,bty='n')

dev.off()







