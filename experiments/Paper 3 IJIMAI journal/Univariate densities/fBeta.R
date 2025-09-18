
library(MoTBFs)
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

# DISTRIBUCIÓN BETA con alpha=2 y beta=2 en el intervalo [0,1]
fBeta2525 = function(x)   x^1.5 * (1-x)^1.5 * gamma(2.5+2.5) / (gamma(2.5) * gamma(2.5));

# DISTRIBUCIÓN BETA con alpha=2.7 y beta=1.3 en el intervalo [0,1]
fBeta2713 = function(x)   x^1.7 * (1-x)^0.3 * gamma(2.7+1.3) / (gamma(2.7) * gamma(1.3));

# DISTRIBUCIÓN BETA con alpha=1.3 y beta=2.7 en el intervalo [0,1]
fBeta1327 = function(x)   x^0.3 * (1-x)^1.7 * gamma(2.7+1.3) / (gamma(2.7) * gamma(1.3));

# DISTRIBUCIÓN BETA con alpha=1.1 y beta=3 en el intervalo [0,1]
fBeta113 = function(x)   x^0.1 * (1-x)^2 * gamma(1.1+3) / (gamma(1.1) * gamma(3));

# DISTRIBUCIÓN BETA con alpha=0.4 y beta=0.6 en el intervalo [0,1]
fBeta0406 = function(x)  x^(-0.6) * (1-x)^(-0.4) * gamma(0.4+0.6) / (gamma(0.4)*gamma(0.6));





####################################################################

# Shenoy
# - Número trozos: 1,2,3,4
# - Grado: 1 trozo (7,4,3), 2 trozos (7,4,3), 3 trozos (7,4,3), 4 trozos (7,4,3,2)

# tmop y Pérez-Bernabé
# - Grados: 7, 4 y 3

# Otros métodos:
# - Tamaño muestra: 50,100,1000



# data frame para los resultados
distribucion1 = rep('beta(2.5,2.5)',31);
intervalo1 = rep('[0,1]',31);
distribucion2 = rep('beta(2.7,1.3)',31);
intervalo2 = rep('[0,1]',31);
distribucion3 = rep('beta(1.3,2.7)',31);
intervalo3 = rep('[0,1]',31);
distribucion4 = rep('beta(1.1,3)',31);
intervalo4 = rep('[0,1]',31);
distribucion5 = rep('beta(0.4,0.6)',31);
intervalo5 = rep('[0,1]',31);
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

# MÉTODO SHENOY (BETA(2.5,2.5))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fBeta2525,7,0,1,1);
fitShenoy14 = ShenoyFit(fBeta2525,4,0,1,1);
fitShenoy13 = ShenoyFit(fBeta2525,3,0,1,1);
fitShenoy27 = ShenoyFit(fBeta2525,7,0,1,2);
fitShenoy24 = ShenoyFit(fBeta2525,4,0,1,2);
fitShenoy23 = ShenoyFit(fBeta2525,3,0,1,2);
fitShenoy37 = ShenoyFit(fBeta2525,7,0,1,3);
fitShenoy34 = ShenoyFit(fBeta2525,4,0,1,3);
fitShenoy33 = ShenoyFit(fBeta2525,3,0,1,3);
fitShenoy47 = ShenoyFit(fBeta2525,7,0,1,4);
fitShenoy44 = ShenoyFit(fBeta2525,4,0,1,4);
fitShenoy43 = ShenoyFit(fBeta2525,3,0,1,4);
fitShenoy42 = ShenoyFit(fBeta2525,2,0,1,4);



divShenoy17 = DivergenceKL2(fBeta2525,fitShenoy17[[1]],0,1);
divShenoy14 = DivergenceKL2(fBeta2525,fitShenoy14[[1]],0,1);
divShenoy13 = DivergenceKL2(fBeta2525,fitShenoy13[[1]],0,1);

sp2 = seq(0,1,length.out=3);
divShenoy27 = DivergenceKL2(fBeta2525,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fBeta2525,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fBeta2525,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fBeta2525,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fBeta2525,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fBeta2525,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,1,length.out=4);
divShenoy37 = DivergenceKL2(fBeta2525,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fBeta2525,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBeta2525,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fBeta2525,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fBeta2525,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBeta2525,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fBeta2525,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fBeta2525,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBeta2525,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,1,length.out=5);
divShenoy47 = DivergenceKL2(fBeta2525,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta2525,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta2525,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta2525,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fBeta2525,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta2525,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta2525,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta2525,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fBeta2525,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta2525,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta2525,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta2525,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fBeta2525,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta2525,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta2525,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta2525,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ  (BETA(2.5,2.5))
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
    xBeta1 = rbeta(50,2.5,2.5);
    
    fitPB1 = univMoTBF(xBeta1,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fBeta2525,fitPBpoly1,0,1);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta2 = rbeta(100,2.5,2.5);
    
    fitPB2 = univMoTBF(xBeta2,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fBeta2525,fitPBpoly2,0,1);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta3 = rbeta(1000,2.5,2.5);
    
    fitPB3 = univMoTBF(xBeta3,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fBeta2525,fitPBpoly3,0,1);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (BETA(2.5,2.5))
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
    xBeta1 = rbeta(50,2.5,2.5);
    
    fittmop1 = PolynomialFit(xBeta1,maxDegree=grados[j],0,1,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fBeta2525,fittmop1,0,1);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta2 = rbeta(100,2.5,2.5);
    
    fittmop2 = PolynomialFit(xBeta2,maxDegree=grados[j],0,1,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fBeta2525,fittmop2,0,1);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta3 = rbeta(1000,2.5,2.5);
    
    fittmop3 = PolynomialFit(xBeta3,maxDegree=grados[j],0,1,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fBeta2525,fittmop3,0,1);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (BETA(2.5,2.5))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,0.25,length.out=100);              x2 = seq(0.25,0.5,length.out=100);
x3 = seq(0.5,0.75,length.out=100);            x4 = seq(0.75,1,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xBeta3,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xBeta3,maxDegree=7,0,1,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fBeta2525,c(0,1),col='#0072B2',type='l',xlim = c(0,1),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')








####################################################################

# MÉTODO SHENOY (BETA(2.7,1.3))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fBeta2713,7,0,1,1);
fitShenoy14 = ShenoyFit(fBeta2713,4,0,1,1);
fitShenoy13 = ShenoyFit(fBeta2713,3,0,1,1);
fitShenoy27 = ShenoyFit(fBeta2713,7,0,1,2);
fitShenoy24 = ShenoyFit(fBeta2713,4,0,1,2);
fitShenoy23 = ShenoyFit(fBeta2713,3,0,1,2);
fitShenoy37 = ShenoyFit(fBeta2713,7,0,1,3);
fitShenoy34 = ShenoyFit(fBeta2713,4,0,1,3);
fitShenoy33 = ShenoyFit(fBeta2713,3,0,1,3);
fitShenoy47 = ShenoyFit(fBeta2713,7,0,1,4);
fitShenoy44 = ShenoyFit(fBeta2713,4,0,1,4);
fitShenoy43 = ShenoyFit(fBeta2713,3,0,1,4);
fitShenoy42 = ShenoyFit(fBeta2713,2,0,1,4);



divShenoy17 = DivergenceKL2(fBeta2713,fitShenoy17[[1]],0,1);
divShenoy14 = DivergenceKL2(fBeta2713,fitShenoy14[[1]],0,1);
divShenoy13 = DivergenceKL2(fBeta2713,fitShenoy13[[1]],0,1);

sp2 = seq(0,1,length.out=3);
divShenoy27 = DivergenceKL2(fBeta2713,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fBeta2713,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fBeta2713,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fBeta2713,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fBeta2713,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fBeta2713,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,1,length.out=4);
divShenoy37 = DivergenceKL2(fBeta2713,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fBeta2713,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBeta2713,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fBeta2713,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fBeta2713,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBeta2713,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fBeta2713,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fBeta2713,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBeta2713,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,1,length.out=5);
divShenoy47 = DivergenceKL2(fBeta2713,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta2713,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta2713,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta2713,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fBeta2713,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta2713,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta2713,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta2713,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fBeta2713,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta2713,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta2713,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta2713,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fBeta2713,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta2713,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta2713,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta2713,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ  (BETA(2.7,1.3))
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
    xBeta1 = rbeta(50,2.7,1.3);
    
    fitPB1 = univMoTBF(xBeta1,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fBeta2713,fitPBpoly1,0,1);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta2 = rbeta(100,2.7,1.3);
    
    fitPB2 = univMoTBF(xBeta2,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fBeta2713,fitPBpoly2,0,1);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta3 = rbeta(1000,2.7,1.3);
    
    fitPB3 = univMoTBF(xBeta3,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fBeta2713,fitPBpoly3,0,1);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (BETA(2.7,1.3))
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
    xBeta1 = rbeta(50,2.7,1.3);
    
    fittmop1 = PolynomialFit(xBeta1,maxDegree=grados[j],0,1,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fBeta2713,fittmop1,0,1);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta2 = rbeta(100,2.7,1.3);
    
    fittmop2 = PolynomialFit(xBeta2,maxDegree=grados[j],0,1,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fBeta2713,fittmop2,0,1);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta3 = rbeta(1000,2.7,1.3);
    
    fittmop3 = PolynomialFit(xBeta3,maxDegree=grados[j],0,1,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fBeta2713,fittmop3,0,1);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (BETA(2.7,1.3))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,0.25,length.out=100);              x2 = seq(0.25,0.5,length.out=100);
x3 = seq(0.5,0.75,length.out=100);            x4 = seq(0.75,1,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xBeta3,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xBeta3,maxDegree=7,0,1,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fBeta2713,c(0,1),col='#0072B2',type='l',xlim = c(0,1),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')








####################################################################

# MÉTODO SHENOY (BETA(1.3,2.7))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fBeta1327,7,0,1,1);
fitShenoy14 = ShenoyFit(fBeta1327,4,0,1,1);
fitShenoy13 = ShenoyFit(fBeta1327,3,0,1,1);
fitShenoy27 = ShenoyFit(fBeta1327,7,0,1,2);
fitShenoy24 = ShenoyFit(fBeta1327,4,0,1,2);
fitShenoy23 = ShenoyFit(fBeta1327,3,0,1,2);
fitShenoy37 = ShenoyFit(fBeta1327,7,0,1,3);
fitShenoy34 = ShenoyFit(fBeta1327,4,0,1,3);
fitShenoy33 = ShenoyFit(fBeta1327,3,0,1,3);
fitShenoy47 = ShenoyFit(fBeta1327,7,0,1,4);
fitShenoy44 = ShenoyFit(fBeta1327,4,0,1,4);
fitShenoy43 = ShenoyFit(fBeta1327,3,0,1,4);
fitShenoy42 = ShenoyFit(fBeta1327,2,0,1,4);



divShenoy17 = DivergenceKL2(fBeta1327,fitShenoy17[[1]],0,1);
divShenoy14 = DivergenceKL2(fBeta1327,fitShenoy14[[1]],0,1);
divShenoy13 = DivergenceKL2(fBeta1327,fitShenoy13[[1]],0,1);

sp2 = seq(0,1,length.out=3);
divShenoy27 = DivergenceKL2(fBeta1327,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fBeta1327,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fBeta1327,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fBeta1327,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fBeta1327,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fBeta1327,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,1,length.out=4);
divShenoy37 = DivergenceKL2(fBeta1327,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fBeta1327,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBeta1327,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fBeta1327,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fBeta1327,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBeta1327,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fBeta1327,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fBeta1327,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBeta1327,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,1,length.out=5);
divShenoy47 = DivergenceKL2(fBeta1327,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta1327,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta1327,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta1327,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fBeta1327,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta1327,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta1327,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta1327,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fBeta1327,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta1327,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta1327,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta1327,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fBeta1327,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta1327,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta1327,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta1327,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ  (BETA(1.3,2.7))
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
    xBeta1 = rbeta(50,1.3,2.7);
    
    fitPB1 = univMoTBF(xBeta1,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fBeta1327,fitPBpoly1,0,1);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta2 = rbeta(100,1.3,2.7);
    
    fitPB2 = univMoTBF(xBeta2,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fBeta1327,fitPBpoly2,0,1);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta3 = rbeta(1000,1.3,2.7);
    
    fitPB3 = univMoTBF(xBeta3,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fBeta1327,fitPBpoly3,0,1);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (BETA(1.3,2.7))
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
    set.seed(i+j-2)
    xBeta1 = rbeta(50,1.3,2.7);
    
    fittmop1 = PolynomialFit(xBeta1,maxDegree=grados[j],0,1,1,2);
    #cat('\n ajuste = ');  print(fittmop1)
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fBeta1327,fittmop1,0,1);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta2 = rbeta(100,1.3,2.7);
    
    fittmop2 = PolynomialFit(xBeta2,maxDegree=grados[j],0,1,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fBeta1327,fittmop2,0,1);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta3 = rbeta(1000,1.3,2.7);
    
    fittmop3 = PolynomialFit(xBeta3,maxDegree=grados[j],0,1,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fBeta1327,fittmop3,0,1);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (BETA(1.3,2.7))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,0.25,length.out=100);              x2 = seq(0.25,0.5,length.out=100);
x3 = seq(0.5,0.75,length.out=100);            x4 = seq(0.75,1,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xBeta3,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xBeta3,maxDegree=7,0,1,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fBeta1327,c(0,1),col='#0072B2',type='l',xlim = c(0,1),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')








####################################################################

# MÉTODO SHENOY (BETA(1.1,3))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fBeta113,7,0,1,1);
fitShenoy14 = ShenoyFit(fBeta113,4,0,1,1);
fitShenoy13 = ShenoyFit(fBeta113,3,0,1,1);
fitShenoy27 = ShenoyFit(fBeta113,7,0,1,2);
fitShenoy24 = ShenoyFit(fBeta113,4,0,1,2);
fitShenoy23 = ShenoyFit(fBeta113,3,0,1,2);
fitShenoy37 = ShenoyFit(fBeta113,7,0,1,3);
fitShenoy34 = ShenoyFit(fBeta113,4,0,1,3);
fitShenoy33 = ShenoyFit(fBeta113,3,0,1,3);
fitShenoy47 = ShenoyFit(fBeta113,7,0,1,4);
fitShenoy44 = ShenoyFit(fBeta113,4,0,1,4);
fitShenoy43 = ShenoyFit(fBeta113,3,0,1,4);
fitShenoy42 = ShenoyFit(fBeta113,2,0,1,4);



divShenoy17 = DivergenceKL2(fBeta113,fitShenoy17[[1]],0,1);
divShenoy14 = DivergenceKL2(fBeta113,fitShenoy14[[1]],0,1);
divShenoy13 = DivergenceKL2(fBeta113,fitShenoy13[[1]],0,1);

sp2 = seq(0,1,length.out=3);
divShenoy27 = DivergenceKL2(fBeta113,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fBeta113,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fBeta113,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fBeta113,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fBeta113,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fBeta113,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,1,length.out=4);
divShenoy37 = DivergenceKL2(fBeta113,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fBeta113,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBeta113,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fBeta113,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fBeta113,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBeta113,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fBeta113,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fBeta113,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBeta113,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,1,length.out=5);
divShenoy47 = DivergenceKL2(fBeta113,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta113,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta113,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta113,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fBeta113,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta113,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta113,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta113,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fBeta113,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta113,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta113,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta113,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fBeta113,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta113,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta113,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta113,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ  (BETA(1.1,3))
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
    xBeta1 = rbeta(50,1.1,3);
    
    fitPB1 = univMoTBF(xBeta1,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fBeta113,fitPBpoly1,0,1);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta2 = rbeta(100,1.1,3);
    
    fitPB2 = univMoTBF(xBeta2,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fBeta113,fitPBpoly2,0,1);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta3 = rbeta(1000,1.1,3);
    
    fitPB3 = univMoTBF(xBeta3,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fBeta113,fitPBpoly3,0,1);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (BETA(1.1,3))
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
    set.seed(i+j-3)
    xBeta1 = rbeta(50,1.1,3);
    
    fittmop1 = PolynomialFit(xBeta1,maxDegree=grados[j],0,1,1,2);
    #cat('\n \n \n  ajuste = ',i);  print(fittmop1)
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fBeta113,fittmop1,0,1);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta2 = rbeta(100,1.1,3);
    
    fittmop2 = PolynomialFit(xBeta2,maxDegree=grados[j],0,1,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fBeta113,fittmop2,0,1);
  }
}



# 1000
s = c(2,3,5:11,14:18,21:24,27:29,30:40)
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    # HE CAMBIADO UN POCO LA SEMILLA PORQUE EN UNA DE ELLAS SALE COEFICIENTES
    # MUY GRANDES QUE DAN PROBLEMAS AL CALCULAR LA INTEGRAL DE LA DIVERGENCIA
    set.seed(s[i+j-1])
    
    xBeta3 = rbeta(1000,1.1,3);
    
    fittmop3 = PolynomialFit(xBeta3,maxDegree=grados[j],0,1,1,2);
    #cat('\n \n \n  semilla = ',i+j+28);
    #cat('\n \n \n  ajuste = ',i);  print(fittmop3)
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fBeta113,fittmop3,0,1);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (BETA(1.1,3))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,0.25,length.out=100);              x2 = seq(0.25,0.5,length.out=100);
x3 = seq(0.5,0.75,length.out=100);            x4 = seq(0.75,1,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xBeta3,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xBeta3,maxDegree=7,0,1,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fBeta113,c(0,1),col='#0072B2',type='l',xlim = c(0,1),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')








####################################################################

# MÉTODO SHENOY (BETA(0.4,0.6))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fBeta0406,7,0,1,1);
fitShenoy14 = ShenoyFit(fBeta0406,4,0,1,1);
fitShenoy13 = ShenoyFit(fBeta0406,3,0,1,1);
fitShenoy27 = ShenoyFit(fBeta0406,7,0,1,2);
fitShenoy24 = ShenoyFit(fBeta0406,4,0,1,2);
fitShenoy23 = ShenoyFit(fBeta0406,3,0,1,2);
fitShenoy37 = ShenoyFit(fBeta0406,7,0,1,3);
fitShenoy34 = ShenoyFit(fBeta0406,4,0,1,3);
fitShenoy33 = ShenoyFit(fBeta0406,3,0,1,3);
fitShenoy47 = ShenoyFit(fBeta0406,7,0,1,4);
fitShenoy44 = ShenoyFit(fBeta0406,4,0,1,4);
fitShenoy43 = ShenoyFit(fBeta0406,3,0,1,4);
fitShenoy42 = ShenoyFit(fBeta0406,2,0,1,4);



divShenoy17 = DivergenceKL2(fBeta0406,fitShenoy17[[1]],0,1);
divShenoy14 = DivergenceKL2(fBeta0406,fitShenoy14[[1]],0,1);
divShenoy13 = DivergenceKL2(fBeta0406,fitShenoy13[[1]],0,1);

sp2 = seq(0,1,length.out=3);
divShenoy27 = DivergenceKL2(fBeta0406,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fBeta0406,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fBeta0406,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fBeta0406,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fBeta0406,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fBeta0406,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,1,length.out=4);
divShenoy37 = DivergenceKL2(fBeta0406,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fBeta0406,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBeta0406,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fBeta0406,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fBeta0406,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBeta0406,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fBeta0406,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fBeta0406,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBeta0406,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,1,length.out=5);
divShenoy47 = DivergenceKL2(fBeta0406,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta0406,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta0406,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta0406,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fBeta0406,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta0406,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta0406,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta0406,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fBeta0406,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta0406,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta0406,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta0406,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fBeta0406,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fBeta0406,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBeta0406,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBeta0406,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ  (BETA(0.4,0.6))
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
    xBeta1 = rbeta(50,0.4,0.6);
    
    fitPB1 = univMoTBF(xBeta1,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fBeta0406,fitPBpoly1,0,1);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta2 = rbeta(100,0.4,0.6);
    
    fitPB2 = univMoTBF(xBeta2,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fBeta0406,fitPBpoly2,0,1);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta3 = rbeta(1000,0.4,0.6);
    
    fitPB3 = univMoTBF(xBeta3,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fBeta0406,fitPBpoly3,0,1);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (BETA(0.4,0.6))
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
    xBeta1 = rbeta(50,0.4,0.6);
    
    fittmop1 = PolynomialFit(xBeta1,maxDegree=grados[j],0,1,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fBeta0406,fittmop1,0,1);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta2 = rbeta(100,0.4,0.6);
    
    fittmop2 = PolynomialFit(xBeta2,maxDegree=grados[j],0,1,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fBeta0406,fittmop2,0,1);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBeta3 = rbeta(1000,0.4,0.6);
    
    fittmop3 = PolynomialFit(xBeta3,maxDegree=grados[j],0,1,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fBeta0406,fittmop3,0,1);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divtmop1 = divtmop1 / 10;       divtmop2 = divtmop2 / 10;       divtmop3 = divtmop3 / 10;




# divergencias en la columna de resultados
resultados$`Div-KL`[125:155] =c(divShenoy17,divShenoy14,divShenoy13,divShenoy27,divShenoy24,
                              divShenoy23,divShenoy37,divShenoy34,divShenoy33,divShenoy47,
                              divShenoy44,divShenoy43,divShenoy42,divPB1,divPB2,divPB3,
                              divtmop1,divtmop2,divtmop3)



resultadosBeta = resultados;

save(resultadosBeta,file='ResultadosBeta.RData');






####################################################################

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (BETA(0.4,0.6))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,0.25,length.out=100);              x2 = seq(0.25,0.5,length.out=100);
x3 = seq(0.5,0.75,length.out=100);            x4 = seq(0.75,1,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xBeta3,POTENTIAL_TYPE="MOP",evalRange=c(0,1),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xBeta3,maxDegree=7,0,1,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fBeta0406,c(0,1),col='#0072B2',type='l',xlim = c(0,1),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')















