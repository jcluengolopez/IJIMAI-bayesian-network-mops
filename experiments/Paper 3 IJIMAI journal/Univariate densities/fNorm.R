
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

# DISTRIBUCIÓN NORMAL (0,1) en el intervalo [-3,3]
fN = function(x)  exp(-x^2/2)/sqrt(2*pi)
iN = pracma::integral(fN,-3,3);
fNorm = function(x)  (exp(-x^2/2)/sqrt(2*pi)) / iN


####################################################################

# Shenoy
# - Número trozos: 1,2,3,4
# - Grado: 1 trozo (7,4,3), 2 trozos (7,4,3), 3 trozos (7,4,3), 4 trozos (7,4,3,2)

# tmop y Pérez-Bernabé
# - Grados: 7, 4 y 3

# Otros métodos:
# - Tamaño muestra: 50,100,1000


# data frame para los resultados
distribucion = rep('N(0,1)',31);
intervalo = rep('[-3,3]',31)
metodo = c(rep('Shenoy',13),rep('Perez-Bernabe',9),rep('tMoPs',9));
nInter = c(rep(1,3),rep(2,3),rep(3,3),rep(4,4),rep(1,9),rep('1-3',9));
grado = c(rep(c(7,4,3),4),2,rep(c(7,4,3),6));
n = c(rep(NA,13),rep(c(50,50,50,100,100,100,1000,1000,1000),2));

resultados = data.frame(distribucion,intervalo,metodo,nInter,grado,n,rep(0,31))

colnames(resultados) = c('f','intervalo','Método','nInter','Grado','Muestra','Div-KL');




####################################################################

# MÉTODO SHENOY


# ponemos el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# carga la función del método
source('ShenoyFit.R')


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fNorm,7,-3,3,1);
fitShenoy14 = ShenoyFit(fNorm,4,-3,3,1);
fitShenoy13 = ShenoyFit(fNorm,3,-3,3,1);
fitShenoy27 = ShenoyFit(fNorm,7,-3,3,2);
fitShenoy24 = ShenoyFit(fNorm,4,-3,3,2);
fitShenoy23 = ShenoyFit(fNorm,3,-3,3,2);
fitShenoy37 = ShenoyFit(fNorm,7,-3,3,3);
fitShenoy34 = ShenoyFit(fNorm,4,-3,3,3);
fitShenoy33 = ShenoyFit(fNorm,3,-3,3,3);
fitShenoy47 = ShenoyFit(fNorm,7,-3,3,4);
fitShenoy44 = ShenoyFit(fNorm,4,-3,3,4);
fitShenoy43 = ShenoyFit(fNorm,3,-3,3,4);
fitShenoy42 = ShenoyFit(fNorm,2,-3,3,4);


divShenoy17 = DivergenceKL2(fNorm,fitShenoy17[[1]],-3,3);
divShenoy14 = DivergenceKL2(fNorm,fitShenoy14[[1]],-3,3);
divShenoy13 = DivergenceKL2(fNorm,fitShenoy13[[1]],-3,3);

sp2 = seq(-3,3,length.out=3);
divShenoy27 = DivergenceKL2(fNorm,fitShenoy27[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fNorm,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fNorm,fitShenoy24[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fNorm,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fNorm,fitShenoy23[[1]],sp2[1],sp2[2]) + 
              DivergenceKL2(fNorm,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(-3,3,length.out=4);
divShenoy37 = DivergenceKL2(fNorm,fitShenoy37[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fNorm,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fNorm,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fNorm,fitShenoy34[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fNorm,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fNorm,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fNorm,fitShenoy33[[1]],sp3[1],sp3[2]) + 
              DivergenceKL2(fNorm,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fNorm,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(-3,3,length.out=5);
divShenoy47 = DivergenceKL2(fNorm,fitShenoy47[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fNorm,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fNorm,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fNorm,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fNorm,fitShenoy44[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fNorm,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fNorm,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fNorm,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fNorm,fitShenoy43[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fNorm,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fNorm,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fNorm,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fNorm,fitShenoy42[[1]],sp4[1],sp4[2]) + 
              DivergenceKL2(fNorm,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fNorm,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fNorm,fitShenoy42[[4]],sp4[4],sp4[5]);






####################################################################

# MÉTODO PÉREZ-BERNABÉ (tiene un intervalo, grados 7,4,3 y tamaño de muestra: 50,100,1000)


# ajustes para distintos tamaños de muestra. Repetir cada uno con 10 muestras
divPB1 = rep(0,3);        divPB2 = rep(0,3);        divPB3 = rep(0,3);
parametros = c(8,5,4)

# 50
for (j in 1:length(parametros))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xNorm1 = rnorm(50,0,1);            xNorm1 = xNorm1[xNorm1>=-3 & xNorm1<=3];
    
    fitPB1 = univMoTBF(xNorm1,POTENTIAL_TYPE="MOP",evalRange=c(-3,3),maxParam=parametros[j]);
    fitPBpoly1 = coef(fitPB1);
    
    divPB1[j] = divPB1[j] + DivergenceKL(fNorm,fitPBpoly1,-3,3);
  }
}


# 100
for (j in 1:length(parametros))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xNorm2 = rnorm(100,0,1);            xNorm2 = xNorm2[xNorm2>=-3 & xNorm2<=3];
    
    fitPB2 = univMoTBF(xNorm2,POTENTIAL_TYPE="MOP",evalRange=c(-3,3),maxParam=parametros[j]);
    fitPBpoly2 = coef(fitPB2);
    
    divPB2[j] = divPB2[j] + DivergenceKL(fNorm,fitPBpoly2,-3,3);
  }
}


# 1000
for (j in 1:length(parametros))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xNorm3 = rnorm(1000,0,1);            xNorm3 = xNorm3[xNorm3>=-3 & xNorm3<=3];
    
    fitPB3 = univMoTBF(xNorm3,POTENTIAL_TYPE="MOP",evalRange=c(-3,3),maxParam=parametros[j]);
    fitPBpoly3 = coef(fitPB3);
    
    divPB3[j] = divPB3[j] + DivergenceKL(fNorm,fitPBpoly3,-3,3);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;








####################################################################

# MÉTODO tMoPs (grados 7,4,3 y tamaño de muestra: 50,100,1000)


# ajustes para distintos tamaños de muestra. Repetir cada uno con 10 muestras
divtmop1 = rep(0,3);        divtmop2 = rep(0,3);        divtmop3 = rep(0,3);
grados = c(7,4,3)

# 50
for (j in 1:length(parametros))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xNorm1 = rnorm(50,0,1);            xNorm1 = xNorm1[xNorm1>=-3 & xNorm1<=3];
    
    fittmop1 = PolynomialFit(xNorm1,maxDegree=grados[j],-3,3,1,2);
    
    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fNorm,fittmop1,-3,3);
  }
}


# 100
for (j in 1:length(parametros))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xNorm2 = rnorm(100,0,1);            xNorm2 = xNorm2[xNorm2>=-3 & xNorm2<=3];
    
    fittmop2 = PolynomialFit(xNorm2,maxDegree=grados[j],-3,3,1,2);
    
    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fNorm,fittmop2,-3,3);
  }
}



# 1000
for (j in 1:length(parametros))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xNorm3 = rnorm(1000,0,1);            xNorm3 = xNorm3[xNorm3>=-3 & xNorm3<=3];
    
    fittmop3 = PolynomialFit(xNorm3,maxDegree=grados[j],-3,3,1,2);
    
    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fNorm,fittmop3,-3,3);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divtmop1 = divtmop1 / 10;       divtmop2 = divtmop2 / 10;       divtmop3 = divtmop3 / 10;




# divergencias en la columna de resultados
resultados$`Div-KL` = c(divShenoy17,divShenoy14,divShenoy13,divShenoy27,divShenoy24,
                        divShenoy23,divShenoy37,divShenoy34,divShenoy33,divShenoy47,
                        divShenoy44,divShenoy43,divShenoy42,divPB1,divPB2,divPB3,
                        divtmop1,divtmop2,divtmop3)

resultadosNorm = resultados;

save(resultadosNorm,file='ResultadosNorm.RData');




####################################################################

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO
# Shenoy 4 intervalos y grado 7
# Pérez-Bernabé con tamaño 1000 y grado 7
# tMoP con tamaño 1000 y grado 7


x1 = seq(-3,-1.5,length.out=100);          x2 = seq(-1.5,0,length.out=100);
x3 = seq(0,1.5,length.out=100);            x4 = seq(1.5,3,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xNorm3,POTENTIAL_TYPE="MOP",evalRange=c(-3,3),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xNorm3,maxDegree=7,-3,3,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fNorm,c(-3,3),col='#0072B2',type='l',xlim = c(-3,3),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')












