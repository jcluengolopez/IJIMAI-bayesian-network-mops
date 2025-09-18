
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

# DISTRIBUCIÓN CHI-2 con k=2 en el intervalo [0,13]
fC2 = function(x)   (1/2)^(2/2)*x^(2/2-1) * exp(-x/2) / gamma(2/2);
iC2 = pracma::integral(fC2,0,13);
fChi2 = function(x)   ((1/2)^(2/2)*x^(2/2-1) * exp(-x/2) / gamma(2/2)) / iC2;

# DISTRIBUCIÓN CHI-2 con k=4 en el intervalo [0,17]
fC4 = function(x)   (1/2)^(4/2)*x^(4/2-1) * exp(-x/2) / gamma(4/2);
iC4 = pracma::integral(fC4,0,17);
fChi4 = function(x)   ((1/2)^(4/2)*x^(4/2-1) * exp(-x/2) / gamma(4/2)) / iC4;

# DISTRIBUCIÓN CHI-2 con k=6 en el intervalo [0,22]
fC6 = function(x)   (1/2)^(6/2)*x^(6/2-1) * exp(-x/2) / gamma(6/2);
iC6 = pracma::integral(fC6,0,22);
fChi6 = function(x)   ((1/2)^(6/2)*x^(6/2-1) * exp(-x/2) / gamma(6/2)) / iC6;

# DISTRIBUCIÓN CHI-2 con k=8 en el intervalo [0,25]
fC8 = function(x)   (1/2)^(8/2)*x^(8/2-1) * exp(-x/2) / gamma(8/2);
iC8 = pracma::integral(fC8,0,25);
fChi8 = function(x)   ((1/2)^(8/2)*x^(8/2-1) * exp(-x/2) / gamma(8/2)) / iC8;

# DISTRIBUCIÓN CHI-2 con k=10 en el intervalo [0,27]
fC10 = function(x)   (1/2)^(10/2)*x^(10/2-1) * exp(-x/2) / gamma(10/2);
iC10 = pracma::integral(fC10,0,27);
fChi10 = function(x)   ((1/2)^(10/2)*x^(10/2-1) * exp(-x/2) / gamma(10/2)) / iC10;



####################################################################

# Shenoy
# - Número trozos: 1,2,3,4
# - Grado: 1 trozo (7,4,3), 2 trozos (7,4,3), 3 trozos (7,4,3), 4 trozos (7,4,3,2)

# tmop y Pérez-Bernabé
# - Grados: 7, 4 y 3

# Otros métodos:
# - Tamaño muestra: 50,100,1000



# data frame para los resultados
distribucion1 = rep('chi-2(2)',31);
intervalo1 = rep('[0,13]',31)
distribucion2 = rep('chi-2(4)',31);
intervalo2 = rep('[0,17]',31)
distribucion3 = rep('chi-2(6)',31);
intervalo3 = rep('[0,22]',31)
distribucion4 = rep('chi-2(8)',31);
intervalo4 = rep('[0,25]',31)
distribucion5 = rep('chi-2(10)',31);
intervalo5 = rep('[0,27]',31)
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

# MÉTODO SHENOY (CHI-2(2))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fChi2,7,0,13,1);
fitShenoy14 = ShenoyFit(fChi2,4,0,13,1);
fitShenoy13 = ShenoyFit(fChi2,3,0,13,1);
fitShenoy27 = ShenoyFit(fChi2,7,0,13,2);
fitShenoy24 = ShenoyFit(fChi2,4,0,13,2);
fitShenoy23 = ShenoyFit(fChi2,3,0,13,2);
fitShenoy37 = ShenoyFit(fChi2,7,0,13,3);
fitShenoy34 = ShenoyFit(fChi2,4,0,13,3);
fitShenoy33 = ShenoyFit(fChi2,3,0,13,3);
fitShenoy47 = ShenoyFit(fChi2,7,0,13,4);
fitShenoy44 = ShenoyFit(fChi2,4,0,13,4);
fitShenoy43 = ShenoyFit(fChi2,3,0,13,4);
fitShenoy42 = ShenoyFit(fChi2,2,0,13,4);



divShenoy17 = DivergenceKL2(fChi2,fitShenoy17[[1]],0,13);
divShenoy14 = DivergenceKL2(fChi2,fitShenoy14[[1]],0,13);
divShenoy13 = DivergenceKL2(fChi2,fitShenoy13[[1]],0,13);

sp2 = seq(0,13,length.out=3);
divShenoy27 = DivergenceKL2(fChi2,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fChi2,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fChi2,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fChi2,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fChi2,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fChi2,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,13,length.out=4);
divShenoy37 = DivergenceKL2(fChi2,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fChi2,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fChi2,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fChi2,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fChi2,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fChi2,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fChi2,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fChi2,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fChi2,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,13,length.out=5);
divShenoy47 = DivergenceKL2(fChi2,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi2,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi2,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi2,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fChi2,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi2,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi2,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi2,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fChi2,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi2,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi2,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi2,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fChi2,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi2,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi2,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi2,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ (CHI-2(2))
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
    xChi1 = rchisq(50,2);           xChi1 = xChi1[xChi1<=13];
    
    fitPB1 = univMoTBF(xChi1,POTENTIAL_TYPE="MOP",evalRange=c(0,13),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fChi2,fitPBpoly1,0,13);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi2 = rchisq(100,2);           xChi2 = xChi2[xChi2<=13];
    
    fitPB2 = univMoTBF(xChi2,POTENTIAL_TYPE="MOP",evalRange=c(0,13),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fChi2,fitPBpoly2,0,13);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi3 = rchisq(1000,2);          xChi3 = xChi3[xChi3<=13];
    
    fitPB3 = univMoTBF(xChi3,POTENTIAL_TYPE="MOP",evalRange=c(0,13),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fChi2,fitPBpoly3,0,13);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs (CHI-2(2))
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
    xChi1 = rchisq(50,2);           xChi1 = xChi1[xChi1<=13];
    
    fittmop1 = PolynomialFit(xChi1,maxDegree=grados[j],0,13,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fChi2,fittmop1,0,13);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi2 = rchisq(100,2);           xChi2 = xChi2[xChi2<=13];
    
    fittmop2 = PolynomialFit(xChi2,maxDegree=grados[j],0,13,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fChi2,fittmop2,0,13);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi3 = rchisq(1000,2);           xChi3 = xChi3[xChi3<=13];
    
    fittmop3 = PolynomialFit(xChi3,maxDegree=grados[j],0,13,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fChi2,fittmop3,0,13);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (CHI-2(2))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,3.25,length.out=100);              x2 = seq(3.25,6.5,length.out=100);
x3 = seq(6.5,9.75,length.out=100);            x4 = seq(9.75,13,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xChi3,POTENTIAL_TYPE="MOP",evalRange=c(0,13),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xChi3,maxDegree=7,0,13,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fChi2,c(0,13),col='#0072B2',type='l',xlim = c(0,13),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')











####################################################################

# MÉTODO SHENOY (CHI-2(4))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fChi4,7,0,17,1);
fitShenoy14 = ShenoyFit(fChi4,4,0,17,1);
fitShenoy13 = ShenoyFit(fChi4,3,0,17,1);
fitShenoy27 = ShenoyFit(fChi4,7,0,17,2);
fitShenoy24 = ShenoyFit(fChi4,4,0,17,2);
fitShenoy23 = ShenoyFit(fChi4,3,0,17,2);
fitShenoy37 = ShenoyFit(fChi4,7,0,17,3);
fitShenoy34 = ShenoyFit(fChi4,4,0,17,3);
fitShenoy33 = ShenoyFit(fChi4,3,0,17,3);
fitShenoy47 = ShenoyFit(fChi4,7,0,17,4);
fitShenoy44 = ShenoyFit(fChi4,4,0,17,4);
fitShenoy43 = ShenoyFit(fChi4,3,0,17,4);
fitShenoy42 = ShenoyFit(fChi4,2,0,17,4);



divShenoy17 = DivergenceKL2(fChi4,fitShenoy17[[1]],0,17);
divShenoy14 = DivergenceKL2(fChi4,fitShenoy14[[1]],0,17);
divShenoy13 = DivergenceKL2(fChi4,fitShenoy13[[1]],0,17);

sp2 = seq(0,17,length.out=3);
divShenoy27 = DivergenceKL2(fChi4,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fChi4,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fChi4,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fChi4,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fChi4,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fChi4,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,17,length.out=4);
divShenoy37 = DivergenceKL2(fChi4,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fChi4,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fChi4,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fChi4,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fChi4,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fChi4,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fChi4,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fChi4,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fChi4,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,17,length.out=5);
divShenoy47 = DivergenceKL2(fChi4,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi4,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi4,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi4,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fChi4,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi4,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi4,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi4,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fChi4,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi4,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi4,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi4,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fChi4,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi4,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi4,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi4,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ (CHI-2(4))
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
    xChi1 = rchisq(50,4);           xChi1 = xChi1[xChi1<=17];
    
    fitPB1 = univMoTBF(xChi1,POTENTIAL_TYPE="MOP",evalRange=c(0,17),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fChi4,fitPBpoly1,0,17);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi2 = rchisq(100,4);           xChi2 = xChi2[xChi2<=17];
    
    fitPB2 = univMoTBF(xChi2,POTENTIAL_TYPE="MOP",evalRange=c(0,17),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fChi4,fitPBpoly2,0,17);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi3 = rchisq(1000,4);          xChi3 = xChi3[xChi3<=17];
    
    fitPB3 = univMoTBF(xChi3,POTENTIAL_TYPE="MOP",evalRange=c(0,17),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fChi4,fitPBpoly3,0,17);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs (CHI-2(4))
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
    xChi1 = rchisq(50,4);           xChi1 = xChi1[xChi1<=17];
    
    fittmop1 = PolynomialFit(xChi1,maxDegree=grados[j],0,17,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fChi4,fittmop1,0,17);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi2 = rchisq(100,4);           xChi2 = xChi2[xChi2<=17];
    
    fittmop2 = PolynomialFit(xChi2,maxDegree=grados[j],0,17,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fChi4,fittmop2,0,17);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi3 = rchisq(1000,4);           xChi3 = xChi3[xChi3<=17];
    
    fittmop3 = PolynomialFit(xChi3,maxDegree=grados[j],0,17,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fChi4,fittmop3,0,17);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (CHI-2(4))
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
fitPB = univMoTBF(xChi3,POTENTIAL_TYPE="MOP",evalRange=c(0,17),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xChi3,maxDegree=7,0,17,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fChi4,c(0,17),col='#0072B2',type='l',xlim = c(0,17),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')








####################################################################

# MÉTODO SHENOY (CHI-2(6))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fChi6,7,0,22,1);
fitShenoy14 = ShenoyFit(fChi6,4,0,22,1);
fitShenoy13 = ShenoyFit(fChi6,3,0,22,1);
fitShenoy27 = ShenoyFit(fChi6,7,0,22,2);
fitShenoy24 = ShenoyFit(fChi6,4,0,22,2);
fitShenoy23 = ShenoyFit(fChi6,3,0,22,2);
fitShenoy37 = ShenoyFit(fChi6,7,0,22,3);
fitShenoy34 = ShenoyFit(fChi6,4,0,22,3);
fitShenoy33 = ShenoyFit(fChi6,3,0,22,3);
fitShenoy47 = ShenoyFit(fChi6,7,0,22,4);
fitShenoy44 = ShenoyFit(fChi6,4,0,22,4);
fitShenoy43 = ShenoyFit(fChi6,3,0,22,4);
fitShenoy42 = ShenoyFit(fChi6,2,0,22,4);



divShenoy17 = DivergenceKL2(fChi6,fitShenoy17[[1]],0,22);
divShenoy14 = DivergenceKL2(fChi6,fitShenoy14[[1]],0,22);
divShenoy13 = DivergenceKL2(fChi6,fitShenoy13[[1]],0,22);

sp2 = seq(0,22,length.out=3);
divShenoy27 = DivergenceKL2(fChi6,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fChi6,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fChi6,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fChi6,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fChi6,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fChi6,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,22,length.out=4);
divShenoy37 = DivergenceKL2(fChi6,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fChi6,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fChi6,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fChi6,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fChi6,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fChi6,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fChi6,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fChi6,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fChi6,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,22,length.out=5);
divShenoy47 = DivergenceKL2(fChi6,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi6,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi6,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi6,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fChi6,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi6,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi6,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi6,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fChi6,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi6,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi6,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi6,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fChi6,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi6,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi6,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi6,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ (CHI-2(6))
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
    xChi1 = rchisq(50,6);           xChi1 = xChi1[xChi1<=22];
    
    fitPB1 = univMoTBF(xChi1,POTENTIAL_TYPE="MOP",evalRange=c(0,22),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fChi6,fitPBpoly1,0,22);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi2 = rchisq(100,6);           xChi2 = xChi2[xChi2<=22];
    
    fitPB2 = univMoTBF(xChi2,POTENTIAL_TYPE="MOP",evalRange=c(0,22),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fChi6,fitPBpoly2,0,22);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi3 = rchisq(1000,6);          xChi3 = xChi3[xChi3<=22];
    
    fitPB3 = univMoTBF(xChi3,POTENTIAL_TYPE="MOP",evalRange=c(0,22),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fChi6,fitPBpoly3,0,22);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs (CHI-2(6))
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
    xChi1 = rchisq(50,6);           xChi1 = xChi1[xChi1<=22];
    
    fittmop1 = PolynomialFit(xChi1,maxDegree=grados[j],0,22,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fChi6,fittmop1,0,22);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi2 = rchisq(100,6);           xChi2 = xChi2[xChi2<=22];
    
    fittmop2 = PolynomialFit(xChi2,maxDegree=grados[j],0,22,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fChi6,fittmop2,0,22);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi3 = rchisq(1000,6);           xChi3 = xChi3[xChi3<=22];
    
    fittmop3 = PolynomialFit(xChi3,maxDegree=grados[j],0,22,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fChi6,fittmop3,0,22);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (CHI-2(6))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,5.5,length.out=100);               x2 = seq(5.5,11,length.out=100);
x3 = seq(11,16.5,length.out=100);              x4 = seq(16.5,22,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xChi3,POTENTIAL_TYPE="MOP",evalRange=c(0,22),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xChi3,maxDegree=7,0,22,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fChi6,c(0,22),col='#0072B2',type='l',xlim = c(0,22),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')








####################################################################

# MÉTODO SHENOY (CHI-2(8))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fChi8,7,0,25,1);
fitShenoy14 = ShenoyFit(fChi8,4,0,25,1);
fitShenoy13 = ShenoyFit(fChi8,3,0,25,1);
fitShenoy27 = ShenoyFit(fChi8,7,0,25,2);
fitShenoy24 = ShenoyFit(fChi8,4,0,25,2);
fitShenoy23 = ShenoyFit(fChi8,3,0,25,2);
fitShenoy37 = ShenoyFit(fChi8,7,0,25,3);
fitShenoy34 = ShenoyFit(fChi8,4,0,25,3);
fitShenoy33 = ShenoyFit(fChi8,3,0,25,3);
fitShenoy47 = ShenoyFit(fChi8,7,0,25,4);
fitShenoy44 = ShenoyFit(fChi8,4,0,25,4);
fitShenoy43 = ShenoyFit(fChi8,3,0,25,4);
fitShenoy42 = ShenoyFit(fChi8,2,0,25,4);



divShenoy17 = DivergenceKL2(fChi8,fitShenoy17[[1]],0,25);
divShenoy14 = DivergenceKL2(fChi8,fitShenoy14[[1]],0,25);
divShenoy13 = DivergenceKL2(fChi8,fitShenoy13[[1]],0,25);

sp2 = seq(0,25,length.out=3);
divShenoy27 = DivergenceKL2(fChi8,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fChi8,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fChi8,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fChi8,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fChi8,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fChi8,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,25,length.out=4);
divShenoy37 = DivergenceKL2(fChi8,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fChi8,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fChi8,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fChi8,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fChi8,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fChi8,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fChi8,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fChi8,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fChi8,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,25,length.out=5);
divShenoy47 = DivergenceKL2(fChi8,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi8,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi8,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi8,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fChi8,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi8,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi8,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi8,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fChi8,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi8,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi8,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi8,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fChi8,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi8,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi8,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi8,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ (CHI-2(8))
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
    xChi1 = rchisq(50,8);           xChi1 = xChi1[xChi1<=25];
    
    fitPB1 = univMoTBF(xChi1,POTENTIAL_TYPE="MOP",evalRange=c(0,25),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fChi8,fitPBpoly1,0,25);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi2 = rchisq(100,8);           xChi2 = xChi2[xChi2<=25];
    
    fitPB2 = univMoTBF(xChi2,POTENTIAL_TYPE="MOP",evalRange=c(0,25),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fChi8,fitPBpoly2,0,25);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi3 = rchisq(1000,8);          xChi3 = xChi3[xChi3<=25];
    
    fitPB3 = univMoTBF(xChi3,POTENTIAL_TYPE="MOP",evalRange=c(0,25),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fChi8,fitPBpoly3,0,25);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs (CHI-2(8))
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
    xChi1 = rchisq(50,8);           xChi1 = xChi1[xChi1<=25];
    
    fittmop1 = PolynomialFit(xChi1,maxDegree=grados[j],0,25,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fChi8,fittmop1,0,25);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi2 = rchisq(100,8);           xChi2 = xChi2[xChi2<=25];
    
    fittmop2 = PolynomialFit(xChi2,maxDegree=grados[j],0,25,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fChi8,fittmop2,0,25);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi3 = rchisq(1000,8);           xChi3 = xChi3[xChi3<=25];
    
    fittmop3 = PolynomialFit(xChi3,maxDegree=grados[j],0,25,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fChi8,fittmop3,0,25);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (CHI-2(8))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,6.25,length.out=100);               x2 = seq(6.25,12.5,length.out=100);
x3 = seq(12.5,18.75,length.out=100);           x4 = seq(18.75,25,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xChi3,POTENTIAL_TYPE="MOP",evalRange=c(0,25),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xChi3,maxDegree=7,0,25,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fChi8,c(0,25),col='#0072B2',type='l',xlim = c(0,25),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')







####################################################################

# MÉTODO SHENOY (CHI-2(10))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fChi10,7,0,27,1);
fitShenoy14 = ShenoyFit(fChi10,4,0,27,1);
fitShenoy13 = ShenoyFit(fChi10,3,0,27,1);
fitShenoy27 = ShenoyFit(fChi10,7,0,27,2);
fitShenoy24 = ShenoyFit(fChi10,4,0,27,2);
fitShenoy23 = ShenoyFit(fChi10,3,0,27,2);
fitShenoy37 = ShenoyFit(fChi10,7,0,27,3);
fitShenoy34 = ShenoyFit(fChi10,4,0,27,3);
fitShenoy33 = ShenoyFit(fChi10,3,0,27,3);
fitShenoy47 = ShenoyFit(fChi10,7,0,27,4);
fitShenoy44 = ShenoyFit(fChi10,4,0,27,4);
fitShenoy43 = ShenoyFit(fChi10,3,0,27,4);
fitShenoy42 = ShenoyFit(fChi10,2,0,27,4);



divShenoy17 = DivergenceKL2(fChi10,fitShenoy17[[1]],0,27);
divShenoy14 = DivergenceKL2(fChi10,fitShenoy14[[1]],0,27);
divShenoy13 = DivergenceKL2(fChi10,fitShenoy13[[1]],0,27);

sp2 = seq(0,27,length.out=3);
divShenoy27 = DivergenceKL2(fChi10,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fChi10,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fChi10,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fChi10,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fChi10,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fChi10,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,27,length.out=4);
divShenoy37 = DivergenceKL2(fChi10,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fChi10,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fChi10,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fChi10,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fChi10,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fChi10,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fChi10,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fChi10,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fChi10,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,27,length.out=5);
divShenoy47 = DivergenceKL2(fChi10,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi10,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi10,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi10,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fChi10,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi10,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi10,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi10,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fChi10,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi10,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi10,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi10,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fChi10,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fChi10,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fChi10,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fChi10,fitShenoy42[[4]],sp4[4],sp4[5]);




####################################################################

# MÉTODO PÉREZ-BERNABÉ (CHI-2(10))
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
    xChi1 = rchisq(50,10);           xChi1 = xChi1[xChi1<=27];
    
    fitPB1 = univMoTBF(xChi1,POTENTIAL_TYPE="MOP",evalRange=c(0,27),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fChi10,fitPBpoly1,0,27);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi2 = rchisq(100,10);           xChi2 = xChi2[xChi2<=27];
    
    fitPB2 = univMoTBF(xChi2,POTENTIAL_TYPE="MOP",evalRange=c(0,27),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fChi10,fitPBpoly2,0,27);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi3 = rchisq(1000,10);          xChi3 = xChi3[xChi3<=27];
    
    fitPB3 = univMoTBF(xChi3,POTENTIAL_TYPE="MOP",evalRange=c(0,27),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fChi10,fitPBpoly3,0,27);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs (CHI-2(10))
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
    xChi1 = rchisq(50,10);           xChi1 = xChi1[xChi1<=27];
    
    fittmop1 = PolynomialFit(xChi1,maxDegree=grados[j],0,27,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fChi10,fittmop1,0,27);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xChi2 = rchisq(100,10);           xChi2 = xChi2[xChi2<=27];
    
    fittmop2 = PolynomialFit(xChi2,maxDegree=grados[j],0,27,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fChi10,fittmop2,0,27);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    # HE CAMBIADO UN POCO LA SEMILLA PORQUE EN UNA DE ELLAS SALE COEFICIENTES
    # MUY GRANDES QUE DAN PROBLEMAS AL CALCULAR LA INTEGRAL DE LA DIVERGENCIA
    set.seed(i+j+4)
    xChi3 = rchisq(1000,10);           xChi3 = xChi3[xChi3<=27];
    
    fittmop3 = PolynomialFit(xChi3,maxDegree=grados[j],0,27,1,2);
    #cat('\n \n \n  ajuste = ',i);  print(fittmop3)
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fChi10,fittmop3,0,27);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divtmop1 = divtmop1 / 10;       divtmop2 = divtmop2 / 10;       divtmop3 = divtmop3 / 10;




# divergencias en la columna de resultados
resultados$`Div-KL`[125:155] =c(divShenoy17,divShenoy14,divShenoy13,divShenoy27,divShenoy24,
                              divShenoy23,divShenoy37,divShenoy34,divShenoy33,divShenoy47,
                              divShenoy44,divShenoy43,divShenoy42,divPB1,divPB2,divPB3,
                              divtmop1,divtmop2,divtmop3)



resultadosChi2 = resultados;

save(resultadosChi2,file='resultadosChi2.RData');





####################################################################

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (CHI-2(10))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,6.75,length.out=100);               x2 = seq(6.75,13.5,length.out=100);
x3 = seq(13.5,20.25,length.out=100);           x4 = seq(20.25,27,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xChi3,POTENTIAL_TYPE="MOP",evalRange=c(0,27),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xChi3,maxDegree=7,0,27,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];
ytmop[length(ytmop)] = ytmop[length(ytmop)-1];


## GUARDAR EN PDF
pdf("GraficatChi2-10.pdf",width=6,height=5)

# distribución original
plot(fChi10,c(0,27),col='#0072B2',type='l',xlim = c(0,27),lwd=4,lty=2,ylab=NA,xlab=NA)
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#a23d75',type='l',lwd=2)
legend(20,0.09, legend=c("Theoretical", "Shenoy","MoTBF","tMoP"),seg.len=2,lwd=2,
       col=c('#0072B2','#D55E00','#009E73','#a23d75'),cex=0.8,bty='n',lty=c(2,1,1,1))

dev.off()















