
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

# DISTRIBUCIÓN EXPONENCIAL de parámetro 1 en el intervalo [0,6]
fE1 = function(x)  exp(-x)
iE1 = pracma::integral(fE1,0,6);
fExp1 = function(x)  exp(-x) / iE1

# DISTRIBUCIÓN EXPONENCIAL de parámetro 3 en el intervalo [0,3]
fE3 = function(x)  3*exp(-3*x)
iE3 = pracma::integral(fE3,0,3);
fExp3 = function(x)  3*exp(-3*x) / iE3

# DISTRIBUCIÓN EXPONENCIAL de parámetro 1/3 en el intervalo [0,20]
fE03 = function(x)  1/3*exp(-x/3)
iE03 = pracma::integral(fE03,0,20);
fExp03 = function(x)  1/3*exp(-x/3) / iE03



####################################################################

# Shenoy
# - Número trozos: 1,2,3,4
# - Grado: 1 trozo (7,4,3), 2 trozos (7,4,3), 3 trozos (7,4,3), 4 trozos (7,4,3,2)

# tmop y Pérez-Bernabé
# - Grados: 7, 4 y 3

# Otros métodos:
# - Tamaño muestra: 50,100,1000





# data frame para los resultados
distribucion1 = rep('exp(1)',31);
intervalo1 = rep('[0,6]',31)
distribucion2 = rep('exp(3)',31);
intervalo2 = rep('[0,3]',31)
distribucion3 = rep('exp(1/3)',31);
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

# MÉTODO SHENOY (EXP(1))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fExp1,7,0,6,1);
fitShenoy14 = ShenoyFit(fExp1,4,0,6,1);
fitShenoy13 = ShenoyFit(fExp1,3,0,6,1);
fitShenoy27 = ShenoyFit(fExp1,7,0,6,2);
fitShenoy24 = ShenoyFit(fExp1,4,0,6,2);
fitShenoy23 = ShenoyFit(fExp1,3,0,6,2);
fitShenoy37 = ShenoyFit(fExp1,7,0,6,3);
fitShenoy34 = ShenoyFit(fExp1,4,0,6,3);
fitShenoy33 = ShenoyFit(fExp1,3,0,6,3);
fitShenoy47 = ShenoyFit(fExp1,7,0,6,4);
fitShenoy44 = ShenoyFit(fExp1,4,0,6,4);
fitShenoy43 = ShenoyFit(fExp1,3,0,6,4);
fitShenoy42 = ShenoyFit(fExp1,2,0,6,4);



divShenoy17 = DivergenceKL2(fExp1,fitShenoy17[[1]],0,6);
divShenoy14 = DivergenceKL2(fExp1,fitShenoy14[[1]],0,6);
divShenoy13 = DivergenceKL2(fExp1,fitShenoy13[[1]],0,6);

sp2 = seq(0,6,length.out=3);
divShenoy27 = DivergenceKL2(fExp1,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fExp1,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fExp1,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fExp1,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fExp1,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fExp1,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,6,length.out=4);
divShenoy37 = DivergenceKL2(fExp1,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fExp1,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fExp1,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fExp1,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fExp1,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fExp1,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fExp1,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fExp1,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fExp1,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,6,length.out=5);
divShenoy47 = DivergenceKL2(fExp1,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fExp1,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fExp1,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fExp1,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fExp1,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fExp1,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fExp1,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fExp1,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fExp1,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fExp1,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fExp1,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fExp1,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fExp1,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fExp1,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fExp1,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fExp1,fitShenoy42[[4]],sp4[4],sp4[5]);



####################################################################

# MÉTODO PÉREZ-BERNABÉ  (EXP(1))
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
    xExp1 = rexp(50,1);            xExp1 = xExp1[xExp1<=6];
    
    fitPB1 = univMoTBF(xExp1,POTENTIAL_TYPE="MOP",evalRange=c(0,6),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fExp1,fitPBpoly1,0,6);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xExp2 = rexp(100,1);            xExp2 = xExp2[xExp2<=6];
    
    fitPB2 = univMoTBF(xExp2,POTENTIAL_TYPE="MOP",evalRange=c(0,6),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fExp1,fitPBpoly2,0,6);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xExp3 = rexp(1000,1);            xExp3 = xExp3[xExp3<=6];
    
    fitPB3 = univMoTBF(xExp3,POTENTIAL_TYPE="MOP",evalRange=c(0,6),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fExp1,fitPBpoly3,0,6);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (EXP(1))
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
    xExp1 = rexp(50,1);            xExp1 = xExp1[xExp1<=6];
    
    fittmop1 = PolynomialFit(xExp1,maxDegree=grados[j],0,6,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fExp1,fittmop1,0,6);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xExp2 = rexp(100,1);            xExp2 = xExp2[xExp2<=6];
    
    fittmop2 = PolynomialFit(xExp2,maxDegree=grados[j],0,6,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fExp1,fittmop2,0,6);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xExp3 = rexp(1000,1);            xExp3 = xExp3[xExp3<=6];
    
    fittmop3 = PolynomialFit(xExp3,maxDegree=grados[j],0,6,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fExp1,fittmop3,0,6);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (EXP(1))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,1.5,length.out=100);          x2 = seq(1.5,3,length.out=100);
x3 = seq(3,4.5,length.out=100);          x4 = seq(4.5,6,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xExp3,POTENTIAL_TYPE="MOP",evalRange=c(0,6),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xExp3,maxDegree=7,0,6,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fExp1,c(0,6),col='#0072B2',type='l',xlim = c(0,6),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')









####################################################################

# MÉTODO SHENOY  (EXP(3))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fExp3,7,0,3,1);
fitShenoy14 = ShenoyFit(fExp3,4,0,3,1);
fitShenoy13 = ShenoyFit(fExp3,3,0,3,1);
fitShenoy27 = ShenoyFit(fExp3,7,0,3,2);
fitShenoy24 = ShenoyFit(fExp3,4,0,3,2);
fitShenoy23 = ShenoyFit(fExp3,3,0,3,2);
fitShenoy37 = ShenoyFit(fExp3,7,0,3,3);
fitShenoy34 = ShenoyFit(fExp3,4,0,3,3);
fitShenoy33 = ShenoyFit(fExp3,3,0,3,3);
fitShenoy47 = ShenoyFit(fExp3,7,0,3,4);
fitShenoy44 = ShenoyFit(fExp3,4,0,3,4);
fitShenoy43 = ShenoyFit(fExp3,3,0,3,4);
fitShenoy42 = ShenoyFit(fExp3,2,0,3,4);



divShenoy17 = DivergenceKL2(fExp3,fitShenoy17[[1]],0,3);
divShenoy14 = DivergenceKL2(fExp3,fitShenoy14[[1]],0,3);
divShenoy13 = DivergenceKL2(fExp3,fitShenoy13[[1]],0,3);

sp2 = seq(0,3,length.out=3);
divShenoy27 = DivergenceKL2(fExp3,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fExp3,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fExp3,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fExp3,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fExp3,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fExp3,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,3,length.out=4);
divShenoy37 = DivergenceKL2(fExp3,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fExp3,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fExp3,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fExp3,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fExp3,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fExp3,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fExp3,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fExp3,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fExp3,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,3,length.out=5);
divShenoy47 = DivergenceKL2(fExp3,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fExp3,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fExp3,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fExp3,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fExp3,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fExp3,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fExp3,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fExp3,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fExp3,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fExp3,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fExp3,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fExp3,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fExp3,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fExp3,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fExp3,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fExp3,fitShenoy42[[4]],sp4[4],sp4[5]);



####################################################################

# MÉTODO PÉREZ-BERNABÉ  (EXP(3))
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
    xExp1 = rexp(50,3);            xExp1 = xExp1[xExp1<=3];
    
    fitPB1 = univMoTBF(xExp1,POTENTIAL_TYPE="MOP",evalRange=c(0,3),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fExp3,fitPBpoly1,0,3);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xExp2 = rexp(100,3);            xExp2 = xExp2[xExp2<=3];
    
    fitPB2 = univMoTBF(xExp2,POTENTIAL_TYPE="MOP",evalRange=c(0,3),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fExp3,fitPBpoly2,0,3);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xExp3 = rexp(1000,3);            xExp3 = xExp3[xExp3<=3];
    
    fitPB3 = univMoTBF(xExp3,POTENTIAL_TYPE="MOP",evalRange=c(0,3),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fExp3,fitPBpoly3,0,3);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (EXP(3))
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
    xExp1 = rexp(50,3);            xExp1 = xExp1[xExp1<=3];
    
    fittmop1 = PolynomialFit(xExp1,maxDegree=grados[j],0,3,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fExp3,fittmop1,0,3);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xExp2 = rexp(100,3);            xExp2 = xExp2[xExp2<=3];
    
    fittmop2 = PolynomialFit(xExp2,maxDegree=grados[j],0,3,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fExp3,fittmop2,0,3);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xExp3 = rexp(1000,3);            xExp3 = xExp3[xExp3<=3];
    
    fittmop3 = PolynomialFit(xExp3,maxDegree=grados[j],0,3,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fExp3,fittmop3,0,3);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (EXP(3))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,0.75,length.out=100);          x2 = seq(0.75,1.5,length.out=100);
x3 = seq(1.5,2.25,length.out=100);        x4 = seq(2.25,3,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xExp3,POTENTIAL_TYPE="MOP",evalRange=c(0,3),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xExp3,maxDegree=7,0,3,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fExp3,c(0,3),col='#0072B2',type='l',xlim = c(0,3),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')








####################################################################

# MÉTODO SHENOY  (EXP(1/3))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fExp03,7,0,20,1);
fitShenoy14 = ShenoyFit(fExp03,4,0,20,1);
fitShenoy13 = ShenoyFit(fExp03,3,0,20,1);
fitShenoy27 = ShenoyFit(fExp03,7,0,20,2);
fitShenoy24 = ShenoyFit(fExp03,4,0,20,2);
fitShenoy23 = ShenoyFit(fExp03,3,0,20,2);
fitShenoy37 = ShenoyFit(fExp03,7,0,20,3);
fitShenoy34 = ShenoyFit(fExp03,4,0,20,3);
fitShenoy33 = ShenoyFit(fExp03,3,0,20,3);
fitShenoy47 = ShenoyFit(fExp03,7,0,20,4);
fitShenoy44 = ShenoyFit(fExp03,4,0,20,4);
fitShenoy43 = ShenoyFit(fExp03,3,0,20,4);
fitShenoy42 = ShenoyFit(fExp03,2,0,20,4);



divShenoy17 = DivergenceKL2(fExp03,fitShenoy17[[1]],0,20);
divShenoy14 = DivergenceKL2(fExp03,fitShenoy14[[1]],0,20);
divShenoy13 = DivergenceKL2(fExp03,fitShenoy13[[1]],0,20);

sp2 = seq(0,20,length.out=3);
divShenoy27 = DivergenceKL2(fExp03,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fExp03,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fExp03,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fExp03,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fExp03,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fExp03,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(0,20,length.out=4);
divShenoy37 = DivergenceKL2(fExp03,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fExp03,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fExp03,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fExp03,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fExp03,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fExp03,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fExp03,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fExp03,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fExp03,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(0,20,length.out=5);
divShenoy47 = DivergenceKL2(fExp03,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fExp03,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fExp03,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fExp03,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fExp03,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fExp03,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fExp03,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fExp03,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fExp03,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fExp03,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fExp03,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fExp03,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fExp03,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fExp03,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fExp03,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fExp03,fitShenoy42[[4]],sp4[4],sp4[5]);



####################################################################

# MÉTODO PÉREZ-BERNABÉ  (EXP(1/3))
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
    xExp1 = rexp(50,1/3);            xExp1 = xExp1[xExp1<=20];
    
    fitPB1 = univMoTBF(xExp1,POTENTIAL_TYPE="MOP",evalRange=c(0,20),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fExp03,fitPBpoly1,0,20);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xExp2 = rexp(100,1/3);            xExp2 = xExp2[xExp2<=20];
    
    fitPB2 = univMoTBF(xExp2,POTENTIAL_TYPE="MOP",evalRange=c(0,20),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fExp03,fitPBpoly2,0,20);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xExp3 = rexp(1000,1/3);            xExp3 = xExp3[xExp3<=20];
    
    fitPB3 = univMoTBF(xExp3,POTENTIAL_TYPE="MOP",evalRange=c(0,20),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fExp03,fitPBpoly3,0,20);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (EXP(1/3))
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
    xExp1 = rexp(50,1/3);            xExp1 = xExp1[xExp1<=20];
    
    fittmop1 = PolynomialFit(xExp1,maxDegree=grados[j],0,20,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fExp03,fittmop1,0,20);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xExp2 = rexp(100,1/3);            xExp2 = xExp2[xExp2<=20];
    
    fittmop2 = PolynomialFit(xExp2,maxDegree=grados[j],0,20,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fExp03,fittmop2,0,20);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xExp3 = rexp(1000,1/3);            xExp3 = xExp3[xExp3<=20];
    
    fittmop3 = PolynomialFit(xExp3,maxDegree=grados[j],0,20,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fExp03,fittmop3,0,20);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divtmop1 = divtmop1 / 10;       divtmop2 = divtmop2 / 10;       divtmop3 = divtmop3 / 10;




# divergencias en la columna de resultados
resultados$`Div-KL`[63:93] = c(divShenoy17,divShenoy14,divShenoy13,divShenoy27,divShenoy24,
                             divShenoy23,divShenoy37,divShenoy34,divShenoy33,divShenoy47,
                             divShenoy44,divShenoy43,divShenoy42,divPB1,divPB2,divPB3,
                             divtmop1,divtmop2,divtmop3)



resultadosExp = resultados;

save(resultadosExp,file='ResultadosExp.RData');





####################################################################

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (EXP(1/3))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(0,5,length.out=100);          x2 = seq(5,10,length.out=100);
x3 = seq(10,15,length.out=100);        x4 = seq(15,20,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xExp3,POTENTIAL_TYPE="MOP",evalRange=c(0,20),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xExp3,maxDegree=7,0,20,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fExp03,c(0,20),col='#0072B2',type='l',xlim = c(0,20),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')









