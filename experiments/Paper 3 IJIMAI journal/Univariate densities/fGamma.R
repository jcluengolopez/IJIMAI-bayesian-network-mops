
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

# DISTRIBUCIÓN GAMMA con lambda=1 y r=6 en el intervalo [0,15]
fG16 = function(x)   x^5 * exp(-x) / 120;
iG16 = pracma::integral(fG16,0,15);
fGam16 = function(x)   (x^5 * exp(-x) / 120) / iG16;

# DISTRIBUCIÓN GAMMA con lambda=1 y r=8 en el intervalo [0,17]
fG18 = function(x)   x^7 * exp(-x) / 5040;
iG18 = pracma::integral(fG18,0,17);
fGam18 = function(x)   (x^7 * exp(-x) / 5040) / iG18;

# DISTRIBUCIÓN GAMMA con lambda=1 y r=11 en el intervalo [0,20]
fG111 = function(x)   x^10 * exp(-x) / 3628800;
iG111 = pracma::integral(fG111,0,20);
fGam111 = function(x)   (x^10 * exp(-x) / 3628800) / iG111;




####################################################################

# Shenoy
# - Número trozos: 1,2,3,4
# - Grado: 1 trozo (7,4,3), 2 trozos (7,4,3), 3 trozos (7,4,3), 4 trozos (7,4,3,2)

# tmop y Pérez-Bernabé
# - Grados: 7, 4 y 3

# Otros métodos:
# - Tamaño muestra: 50,100,1000



# data frame para los resultados
distribucion1 = rep('gamma(6,1)',31);
intervalo1 = rep('[0,15]',31)
distribucion2 = rep('gamma(8,1)',31);
intervalo2 = rep('[0,17]',31)
distribucion3 = rep('gamma(11,1)',31);
intervalo3 = rep('[0,20]',31)
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

# MÉTODO SHENOY (GAMMA(1,6))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fGam16,7,0,15,1);
fitShenoy14 = ShenoyFit(fGam16,4,0,15,1);
fitShenoy13 = ShenoyFit(fGam16,3,0,15,1);
fitShenoy27 = ShenoyFit(fGam16,7,0,15,2);
fitShenoy24 = ShenoyFit(fGam16,4,0,15,2);
fitShenoy23 = ShenoyFit(fGam16,3,0,15,2);
fitShenoy37 = ShenoyFit(fGam16,7,0,15,3);
fitShenoy34 = ShenoyFit(fGam16,4,0,15,3);
fitShenoy33 = ShenoyFit(fGam16,3,0,15,3);
fitShenoy47 = ShenoyFit(fGam16,7,0,15,4);
fitShenoy44 = ShenoyFit(fGam16,4,0,15,4);
fitShenoy43 = ShenoyFit(fGam16,3,0,15,4);
fitShenoy42 = ShenoyFit(fGam16,2,0,15,4);



divShenoy17 = DivergenceKL2(fGam16,fitShenoy17[[1]],0,15);
divShenoy14 = DivergenceKL2(fGam16,fitShenoy14[[1]],0,15);
divShenoy13 = DivergenceKL2(fGam16,fitShenoy13[[1]],0,15);

sp2 = seq(0,15,length.out=3);
divShenoy27 = DivergenceKL2(fGam16,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fGam16,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fGam16,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fGam16,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fGam16,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fGam16,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,15,length.out=4);
divShenoy37 = DivergenceKL2(fGam16,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fGam16,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fGam16,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fGam16,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fGam16,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fGam16,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fGam16,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fGam16,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fGam16,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,15,length.out=5);
divShenoy47 = DivergenceKL2(fGam16,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fGam16,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fGam16,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fGam16,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fGam16,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fGam16,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fGam16,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fGam16,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fGam16,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fGam16,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fGam16,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fGam16,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fGam16,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fGam16,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fGam16,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fGam16,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ  (GAMMA(1,6))
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
    xGam1 = rgamma(50,shape=6,scale=1);           xGam1 = xGam1[xGam1<=15];
    
    fitPB1 = univMoTBF(xGam1,POTENTIAL_TYPE="MOP",evalRange=c(0,15),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fGam16,fitPBpoly1,0,15);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xGam2 = rgamma(100,shape=6,scale=1);           xGam2 = xGam2[xGam2<=15];
    
    fitPB2 = univMoTBF(xGam2,POTENTIAL_TYPE="MOP",evalRange=c(0,15),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fGam16,fitPBpoly2,0,15);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xGam3 = rgamma(1000,shape=6,scale=1);           xGam3 = xGam3[xGam3<=15];
    
    fitPB3 = univMoTBF(xGam3,POTENTIAL_TYPE="MOP",evalRange=c(0,15),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fGam16,fitPBpoly3,0,15);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (GAMMA(1,6))
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
    xGam1 = rgamma(50,shape=6,scale=1);           xGam1 = xGam1[xGam1<=15];
    
    fittmop1 = PolynomialFit(xGam1,maxDegree=grados[j],0,15,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fGam16,fittmop1,0,15);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xGam2 = rgamma(100,shape=6,scale=1);           xGam2 = xGam2[xGam2<=15];
    
    fittmop2 = PolynomialFit(xGam2,maxDegree=grados[j],0,15,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fGam16,fittmop2,0,15);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xGam3 = rgamma(1000,shape=6,scale=1);           xGam3 = xGam3[xGam3<=15];
    
    fittmop3 = PolynomialFit(xGam3,maxDegree=grados[j],0,15,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fGam16,fittmop3,0,15);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (GAMMA(1,6))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,3.75,length.out=100);               x2 = seq(3.75,7.5,length.out=100);
x3 = seq(7.5,11.25,length.out=100);            x4 = seq(11.25,15,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xGam3,POTENTIAL_TYPE="MOP",evalRange=c(0,15),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xGam3,maxDegree=7,0,15,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fGam16,c(0,15),col='#0072B2',type='l',xlim = c(0,15),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')










####################################################################

# MÉTODO SHENOY (GAMMA(1,8))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fGam18,7,0,17,1);
fitShenoy14 = ShenoyFit(fGam18,4,0,17,1);
fitShenoy13 = ShenoyFit(fGam18,3,0,17,1);
fitShenoy27 = ShenoyFit(fGam18,7,0,17,2);
fitShenoy24 = ShenoyFit(fGam18,4,0,17,2);
fitShenoy23 = ShenoyFit(fGam18,3,0,17,2);
fitShenoy37 = ShenoyFit(fGam18,7,0,17,3);
fitShenoy34 = ShenoyFit(fGam18,4,0,17,3);
fitShenoy33 = ShenoyFit(fGam18,3,0,17,3);
fitShenoy47 = ShenoyFit(fGam18,7,0,17,4);
fitShenoy44 = ShenoyFit(fGam18,4,0,17,4);
fitShenoy43 = ShenoyFit(fGam18,3,0,17,4);
fitShenoy42 = ShenoyFit(fGam18,2,0,17,4);



divShenoy17 = DivergenceKL2(fGam18,fitShenoy17[[1]],0,17);
divShenoy14 = DivergenceKL2(fGam18,fitShenoy14[[1]],0,17);
divShenoy13 = DivergenceKL2(fGam18,fitShenoy13[[1]],0,17);

sp2 = seq(0,17,length.out=3);
divShenoy27 = DivergenceKL2(fGam18,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fGam18,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fGam18,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fGam18,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fGam18,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fGam18,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,17,length.out=4);
divShenoy37 = DivergenceKL2(fGam18,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fGam18,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fGam18,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fGam18,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fGam18,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fGam18,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fGam18,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fGam18,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fGam18,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,17,length.out=5);
divShenoy47 = DivergenceKL2(fGam18,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fGam18,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fGam18,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fGam18,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fGam18,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fGam18,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fGam18,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fGam18,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fGam18,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fGam18,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fGam18,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fGam18,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fGam18,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fGam18,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fGam18,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fGam18,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ  (GAMMA(1,8))
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
    xGam1 = rgamma(50,shape=8,scale=1);           xGam1 = xGam1[xGam1<=17];
    
    fitPB1 = univMoTBF(xGam1,POTENTIAL_TYPE="MOP",evalRange=c(0,17),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fGam18,fitPBpoly1,0,17);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xGam2 = rgamma(100,shape=8,scale=1);           xGam2 = xGam2[xGam2<=17];
    
    fitPB2 = univMoTBF(xGam2,POTENTIAL_TYPE="MOP",evalRange=c(0,17),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fGam18,fitPBpoly2,0,17);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xGam3 = rgamma(1000,shape=8,scale=1);           xGam3 = xGam3[xGam3<=17];
    
    fitPB3 = univMoTBF(xGam3,POTENTIAL_TYPE="MOP",evalRange=c(0,17),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fGam18,fitPBpoly3,0,17);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (GAMMA(1,8))
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
    xGam1 = rgamma(50,shape=8,scale=1);           xGam1 = xGam1[xGam1<=17];
    
    fittmop1 = PolynomialFit(xGam1,maxDegree=grados[j],0,17,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fGam18,fittmop1,0,17);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xGam2 = rgamma(100,shape=8,scale=1);           xGam2 = xGam2[xGam2<=17];
    
    fittmop2 = PolynomialFit(xGam2,maxDegree=grados[j],0,17,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fGam18,fittmop2,0,17);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xGam3 = rgamma(1000,shape=8,scale=1);           xGam3 = xGam3[xGam3<=17];
    
    fittmop3 = PolynomialFit(xGam3,maxDegree=grados[j],0,17,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fGam18,fittmop3,0,17);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (GAMMA(1,8))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,4.25,length.out=100);               x2 = seq(4.25,8.5,length.out=100);
x3 = seq(8.5,12.75,length.out=100);            x4 = seq(12.75,17,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xGam3,POTENTIAL_TYPE="MOP",evalRange=c(0,17),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xGam3,maxDegree=7,0,17,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fGam18,c(0,17),col='#0072B2',type='l',xlim = c(0,17),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')









####################################################################

# MÉTODO SHENOY (GAMMA(1,11))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fGam111,7,0,20,1);
fitShenoy14 = ShenoyFit(fGam111,4,0,20,1);
fitShenoy13 = ShenoyFit(fGam111,3,0,20,1);
fitShenoy27 = ShenoyFit(fGam111,7,0,20,2);
fitShenoy24 = ShenoyFit(fGam111,4,0,20,2);
fitShenoy23 = ShenoyFit(fGam111,3,0,20,2);
fitShenoy37 = ShenoyFit(fGam111,7,0,20,3);
fitShenoy34 = ShenoyFit(fGam111,4,0,20,3);
fitShenoy33 = ShenoyFit(fGam111,3,0,20,3);
fitShenoy47 = ShenoyFit(fGam111,7,0,20,4);
fitShenoy44 = ShenoyFit(fGam111,4,0,20,4);
fitShenoy43 = ShenoyFit(fGam111,3,0,20,4);
fitShenoy42 = ShenoyFit(fGam111,2,0,20,4);



divShenoy17 = DivergenceKL2(fGam111,fitShenoy17[[1]],0,20);
divShenoy14 = DivergenceKL2(fGam111,fitShenoy14[[1]],0,20);
divShenoy13 = DivergenceKL2(fGam111,fitShenoy13[[1]],0,20);

sp2 = seq(0,20,length.out=3);
divShenoy27 = DivergenceKL2(fGam111,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fGam111,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fGam111,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fGam111,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fGam111,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fGam111,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,20,length.out=4);
divShenoy37 = DivergenceKL2(fGam111,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fGam111,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fGam111,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fGam111,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fGam111,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fGam111,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fGam111,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fGam111,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fGam111,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,20,length.out=5);
divShenoy47 = DivergenceKL2(fGam111,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fGam111,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fGam111,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fGam111,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fGam111,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fGam111,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fGam111,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fGam111,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fGam111,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fGam111,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fGam111,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fGam111,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fGam111,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fGam111,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fGam111,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fGam111,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ  (GAMMA(1,11))
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
    xGam1 = rgamma(50,shape=11,scale=1);           xGam1 = xGam1[xGam1<=20];
    
    fitPB1 = univMoTBF(xGam1,POTENTIAL_TYPE="MOP",evalRange=c(0,20),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fGam111,fitPBpoly1,0,20);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xGam2 = rgamma(100,shape=11,scale=1);           xGam2 = xGam2[xGam2<=20];
    
    fitPB2 = univMoTBF(xGam2,POTENTIAL_TYPE="MOP",evalRange=c(0,20),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fGam111,fitPBpoly2,0,20);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xGam3 = rgamma(1000,shape=11,scale=1);           xGam3 = xGam3[xGam3<=20];
    
    fitPB3 = univMoTBF(xGam3,POTENTIAL_TYPE="MOP",evalRange=c(0,20),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fGam111,fitPBpoly3,0,20);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (GAMMA(1,11))
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
    xGam1 = rgamma(50,shape=11,scale=1);           xGam1 = xGam1[xGam1<=20];
    
    fittmop1 = PolynomialFit(xGam1,maxDegree=grados[j],0,20,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fGam111,fittmop1,0,20);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xGam2 = rgamma(100,shape=11,scale=1);           xGam2 = xGam2[xGam2<=20];
    
    fittmop2 = PolynomialFit(xGam2,maxDegree=grados[j],0,20,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fGam111,fittmop2,0,20);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xGam3 = rgamma(1000,shape=11,scale=1);           xGam3 = xGam3[xGam3<=20];
    
    fittmop3 = PolynomialFit(xGam3,maxDegree=grados[j],0,20,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fGam111,fittmop3,0,20);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divtmop1 = divtmop1 / 10;       divtmop2 = divtmop2 / 10;       divtmop3 = divtmop3 / 10;




# divergencias en la columna de resultados
resultados$`Div-KL`[63:93] = c(divShenoy17,divShenoy14,divShenoy13,divShenoy27,divShenoy24,
                              divShenoy23,divShenoy37,divShenoy34,divShenoy33,divShenoy47,
                              divShenoy44,divShenoy43,divShenoy42,divPB1,divPB2,divPB3,
                              divtmop1,divtmop2,divtmop3)



resultadosGam = resultados;

save(resultadosGam,file='ResultadosGam.RData');





####################################################################

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (GAMMA(1,11))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,5,length.out=100);              x2 = seq(5,10,length.out=100);
x3 = seq(10,15,length.out=100);            x4 = seq(15,20,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xGam3,POTENTIAL_TYPE="MOP",evalRange=c(0,20),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xGam3,maxDegree=7,0,20,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fGam111,c(0,20),col='#0072B2',type='l',xlim = c(0,20),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')















