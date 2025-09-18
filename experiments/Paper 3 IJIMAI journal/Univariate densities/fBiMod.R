
library(MoTBFs)
library(pracma)
library(polynom)
library(NlcOptim)
library(FamilyRank)  # simular la bimodal


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



# DISTRIBUCIÓN BIMODAL DE NORMAL (0,1) Y NORMAL (6,2) en el intervalo [-3,12]
fBM0162 = function(x)  0.25*exp(-x^2/2)/sqrt(2*pi) +
                          0.75*exp(-(x-6)^2/(2*2^2))/(2*sqrt(2*pi))
iBM0162 = pracma::integral(fBM0162,-3,12);
fBiMod0162 = function(x)  (0.25*exp(-x^2/2)/sqrt(2*pi) +
                           0.75*exp(-(x-6)^2/(2*2^2))/(2*sqrt(2*pi))) / iBM0162

# DISTRIBUCIÓN BIMODAL DE NORMAL (1,1) Y NORMAL (11,3) en el intervalo [-2,20]
fBM11113 = function(x)  0.4*exp(-(x-1)^2/2)/sqrt(2*pi) +
                           0.6*exp(-(x-11)^2/(2*3^2))/(3*sqrt(2*pi))
iBM11113 = pracma::integral(fBM11113,-2,20);
fBiMod11113 = function(x)  (0.4*exp(-(x-1)^2/2)/sqrt(2*pi) +
                            0.6*exp(-(x-11)^2/(2*3^2))/(3*sqrt(2*pi))) / iBM11113


# DISTRIBUCIÓN BIMODAL DE NORMAL (1,1) Y NORMAL (7,1) en el intervalo [-2,10]
fBM1171 = function(x)  0.5*exp(-(x-1)^2/2)/sqrt(2*pi) +
                          0.5*exp(-(x-7)^2/2)/sqrt(2*pi)
iBM1171 = pracma::integral(fBM1171,-2,10);
fBiMod1171 = function(x)  (0.5*exp(-(x-1)^2/2)/sqrt(2*pi) +
                           0.5*exp(-(x-7)^2/2)/sqrt(2*pi)) / iBM1171


####################################################################

# Shenoy
# - Número trozos: 1,2,3,4
# - Grado: 1 trozo (7,4,3), 2 trozos (7,4,3), 3 trozos (7,4,3), 4 trozos (7,4,3,2)

# tmop y Pérez-Bernabé
# - Grados: 7, 4 y 3

# Otros métodos:
# - Tamaño muestra: 50,100,1000



# data frame para los resultados
distribucion1 = rep('bimod(0,1)(6,2)',31);
intervalo1 = rep('[0,15]',31)
distribucion2 = rep('bimod(1,1)(11,3)',31);
intervalo2 = rep('[0,17]',31)
distribucion3 = rep('bimod(1,1)(7,1)',31);
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

# MÉTODO SHENOY  (BIMODAL(0,1)(6,2))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fBiMod0162,7,-3,12,1);
fitShenoy14 = ShenoyFit(fBiMod0162,4,-3,12,1);
fitShenoy13 = ShenoyFit(fBiMod0162,3,-3,12,1);
fitShenoy27 = ShenoyFit(fBiMod0162,7,-3,12,2);
fitShenoy24 = ShenoyFit(fBiMod0162,4,-3,12,2);
fitShenoy23 = ShenoyFit(fBiMod0162,3,-3,12,2);
fitShenoy37 = ShenoyFit(fBiMod0162,7,-3,12,3);
fitShenoy34 = ShenoyFit(fBiMod0162,4,-3,12,3);
fitShenoy33 = ShenoyFit(fBiMod0162,3,-3,12,3);

tshenoy1 = proc.time()
fitShenoy47 = ShenoyFit(fBiMod0162,7,-3,12,4);
tshenoy2 = proc.time()

fitShenoy44 = ShenoyFit(fBiMod0162,4,-3,12,4);
fitShenoy43 = ShenoyFit(fBiMod0162,3,-3,12,4);
fitShenoy42 = ShenoyFit(fBiMod0162,2,-3,12,4);



divShenoy17 = DivergenceKL2(fBiMod0162,fitShenoy17[[1]],-3,12);
divShenoy14 = DivergenceKL2(fBiMod0162,fitShenoy14[[1]],-3,12);
divShenoy13 = DivergenceKL2(fBiMod0162,fitShenoy13[[1]],-3,12);

sp2 = seq(-3,12,length.out=3);
divShenoy27 = DivergenceKL2(fBiMod0162,fitShenoy27[[1]],sp2[1],sp2[2]) +
              DivergenceKL2(fBiMod0162,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fBiMod0162,fitShenoy24[[1]],sp2[1],sp2[2]) +
              DivergenceKL2(fBiMod0162,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fBiMod0162,fitShenoy23[[1]],sp2[1],sp2[2]) +
              DivergenceKL2(fBiMod0162,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(-3,12,length.out=4);
divShenoy37 = DivergenceKL2(fBiMod0162,fitShenoy37[[1]],sp3[1],sp3[2]) +
              DivergenceKL2(fBiMod0162,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBiMod0162,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fBiMod0162,fitShenoy34[[1]],sp3[1],sp3[2]) +
              DivergenceKL2(fBiMod0162,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBiMod0162,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fBiMod0162,fitShenoy33[[1]],sp3[1],sp3[2]) +
              DivergenceKL2(fBiMod0162,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBiMod0162,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(-3,12,length.out=5);
divShenoy47 = DivergenceKL2(fBiMod0162,fitShenoy47[[1]],sp4[1],sp4[2]) +
              DivergenceKL2(fBiMod0162,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBiMod0162,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBiMod0162,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fBiMod0162,fitShenoy44[[1]],sp4[1],sp4[2]) +
              DivergenceKL2(fBiMod0162,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBiMod0162,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBiMod0162,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fBiMod0162,fitShenoy43[[1]],sp4[1],sp4[2]) +
              DivergenceKL2(fBiMod0162,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBiMod0162,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBiMod0162,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fBiMod0162,fitShenoy42[[1]],sp4[1],sp4[2]) +
              DivergenceKL2(fBiMod0162,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBiMod0162,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBiMod0162,fitShenoy42[[4]],sp4[4],sp4[5]);





####################################################################

# MÉTODO PÉREZ-BERNABÉ  (BIMODAL(0,1)(6,2))
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
    xBiMod1 = rbinorm(50,0,6,1,2,0.25);         xBiMod1 = xBiMod1[xBiMod1>=-3 & xBiMod1<=12];

    fitPB1 = univMoTBF(xBiMod1,POTENTIAL_TYPE="MOP",evalRange=c(-3,12),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);

    divPB1[j] = divPB1[j] + DivergenceKL(fBiMod0162,fitPBpoly1,-3,12);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBiMod2 = rbinorm(100,0,6,1,2,0.25);      xBiMod2 = xBiMod2[xBiMod2>=-3 & xBiMod2<=12];

    fitPB2 = univMoTBF(xBiMod2,POTENTIAL_TYPE="MOP",evalRange=c(-3,12),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);

    divPB2[j] = divPB2[j] + DivergenceKL(fBiMod0162,fitPBpoly2,-3,12);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBiMod3 = rbinorm(1000,0,6,1,2,0.25);       xBiMod3 = xBiMod3[xBiMod3>=-3 & xBiMod3<=12];

    fitPB3 = univMoTBF(xBiMod3,POTENTIAL_TYPE="MOP",evalRange=c(-3,12),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);

    divPB3[j] = divPB3[j] + DivergenceKL(fBiMod0162,fitPBpoly3,-3,12);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (BIMODAL(0,1)(6,2))
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
    xBiMod1 = rbinorm(50,0,6,1,2,0.25);         xBiMod1 = xBiMod1[xBiMod1>=-3 & xBiMod1<=12];

    fittmop1 = PolynomialFit(xBiMod1,maxDegree=grados[j],-3,12,1,2);

    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fBiMod0162,fittmop1,-3,12);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBiMod2 = rbinorm(100,0,6,1,2,0.25);       xBiMod2 = xBiMod2[xBiMod2>=-3 & xBiMod2<=12];

    fittmop2 = PolynomialFit(xBiMod2,maxDegree=grados[j],-3,12,1,2);

    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fBiMod0162,fittmop2,-3,12);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBiMod3 = rbinorm(1000,0,6,1,2,0.25);      xBiMod3 = xBiMod3[xBiMod3>=-3 & xBiMod3<=12];

    fittmop3 = PolynomialFit(xBiMod3,maxDegree=grados[j],-3,12,1,2);

    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fBiMod0162,fittmop3,-3,12);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (BIMODAL(0,1)(6,2))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(-3,0.75,length.out=100);            x2 = seq(0.75,4.5,length.out=100);
x3 = seq(4.5,8.25,length.out=100);           x4 = seq(8.25,12,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
tPB1 = proc.time()
fitPB = univMoTBF(xBiMod3,POTENTIAL_TYPE="MOP",evalRange=c(-3,12),maxParam=8);
tPB2 = proc.time()
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
ttmop1 = proc.time()
fittmop = PolynomialFit(xBiMod3,maxDegree=7,-3,12,1,2);
ttmop2 = proc.time()
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# comparación de tiempos
cat('\n tiempo Shenoy = ',tshenoy2-tshenoy1)
cat('\n tiempo Pérez-Bernabé = ',tPB2-tPB1)
cat('\n tiempo tmop = ',ttmop2-ttmop1)



## GUARDAR EN PDF
pdf("GraficaBiMod01-62.pdf",width=6,height=5)

# distribución original
plot(fBiMod0162,c(-3,12),col='#0072B2',type='l',xlim=c(-3,12),ylim=c(0,0.16),lwd=4,lty=2,
     ylab=NA,xlab=NA)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#a23d75',type='l',lwd=2)
legend(8,0.16, legend=c("Theoretical", "Shenoy","MoTBF","tMoP"),seg.len=2,lwd=2,
       col=c('#0072B2','#D55E00','#009E73','#a23d75'),cex=0.8,bty='n',lty=c(2,1,1,1))

dev.off()








####################################################################

# MÉTODO SHENOY  (BIMODAL(1,1)(11,3))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fBiMod11113,7,-2,20,1);
fitShenoy14 = ShenoyFit(fBiMod11113,4,-2,20,1);
fitShenoy13 = ShenoyFit(fBiMod11113,3,-2,20,1);
fitShenoy27 = ShenoyFit(fBiMod11113,7,-2,20,2);
fitShenoy24 = ShenoyFit(fBiMod11113,4,-2,20,2);
fitShenoy23 = ShenoyFit(fBiMod11113,3,-2,20,2);
fitShenoy37 = ShenoyFit(fBiMod11113,7,-2,20,3);
fitShenoy34 = ShenoyFit(fBiMod11113,4,-2,20,3);
fitShenoy33 = ShenoyFit(fBiMod11113,3,-2,20,3);
fitShenoy47 = ShenoyFit(fBiMod11113,7,-2,20,4);
fitShenoy44 = ShenoyFit(fBiMod11113,4,-2,20,4);
fitShenoy43 = ShenoyFit(fBiMod11113,3,-2,20,4);
fitShenoy42 = ShenoyFit(fBiMod11113,2,-2,20,4);



divShenoy17 = DivergenceKL2(fBiMod11113,fitShenoy17[[1]],-2,20);
divShenoy14 = DivergenceKL2(fBiMod11113,fitShenoy14[[1]],-2,20);
divShenoy13 = DivergenceKL2(fBiMod11113,fitShenoy13[[1]],-2,20);

sp2 = seq(-2,20,length.out=3);
divShenoy27 = DivergenceKL2(fBiMod11113,fitShenoy27[[1]],sp2[1],sp2[2]) +
              DivergenceKL2(fBiMod11113,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fBiMod11113,fitShenoy24[[1]],sp2[1],sp2[2]) +
              DivergenceKL2(fBiMod11113,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fBiMod11113,fitShenoy23[[1]],sp2[1],sp2[2]) +
              DivergenceKL2(fBiMod11113,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(-2,20,length.out=4);
divShenoy37 = DivergenceKL2(fBiMod11113,fitShenoy37[[1]],sp3[1],sp3[2]) +
              DivergenceKL2(fBiMod11113,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBiMod11113,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fBiMod11113,fitShenoy34[[1]],sp3[1],sp3[2]) +
              DivergenceKL2(fBiMod11113,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBiMod11113,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fBiMod11113,fitShenoy33[[1]],sp3[1],sp3[2]) +
              DivergenceKL2(fBiMod11113,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBiMod11113,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(-2,20,length.out=5);
divShenoy47 = DivergenceKL2(fBiMod11113,fitShenoy47[[1]],sp4[1],sp4[2]) +
              DivergenceKL2(fBiMod11113,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBiMod11113,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBiMod11113,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fBiMod11113,fitShenoy44[[1]],sp4[1],sp4[2]) +
              DivergenceKL2(fBiMod11113,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBiMod11113,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBiMod11113,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fBiMod11113,fitShenoy43[[1]],sp4[1],sp4[2]) +
              DivergenceKL2(fBiMod11113,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBiMod11113,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBiMod11113,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fBiMod11113,fitShenoy42[[1]],sp4[1],sp4[2]) +
              DivergenceKL2(fBiMod11113,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBiMod11113,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBiMod11113,fitShenoy42[[4]],sp4[4],sp4[5]);





####################################################################

# MÉTODO PÉREZ-BERNABÉ  (BIMODAL(1,1)(11,3))
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
    xBiMod1 = rbinorm(50,1,11,1,3,0.4);         xBiMod1 = xBiMod1[xBiMod1>=-2 & xBiMod1<=20];

    fitPB1 = univMoTBF(xBiMod1,POTENTIAL_TYPE="MOP",evalRange=c(-2,20),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);

    divPB1[j] = divPB1[j] + DivergenceKL(fBiMod11113,fitPBpoly1,-2,20);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBiMod2 = rbinorm(100,1,11,1,3,0.4);      xBiMod2 = xBiMod2[xBiMod2>=-2 & xBiMod2<=20];

    fitPB2 = univMoTBF(xBiMod2,POTENTIAL_TYPE="MOP",evalRange=c(-2,20),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);

    divPB2[j] = divPB2[j] + DivergenceKL(fBiMod11113,fitPBpoly2,-2,20);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBiMod3 = rbinorm(1000,1,11,1,3,0.4);       xBiMod3 = xBiMod3[xBiMod3>=-2 & xBiMod3<=20];

    fitPB3 = univMoTBF(xBiMod3,POTENTIAL_TYPE="MOP",evalRange=c(-2,20),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);

    divPB3[j] = divPB3[j] + DivergenceKL(fBiMod11113,fitPBpoly3,-2,20);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (BIMODAL(1,1)(11,3))
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
    xBiMod1 = rbinorm(50,1,11,1,3,0.4);         xBiMod1 = xBiMod1[xBiMod1>=-2 & xBiMod1<=20];

    fittmop1 = PolynomialFit(xBiMod1,maxDegree=grados[j],-2,20,1,2);

    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fBiMod11113,fittmop1,-2,20);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBiMod2 = rbinorm(100,1,11,1,3,0.4);       xBiMod2 = xBiMod2[xBiMod2>=-2 & xBiMod2<=20];

    fittmop2 = PolynomialFit(xBiMod2,maxDegree=grados[j],-2,20,1,2);

    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fBiMod11113,fittmop2,-2,20);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBiMod3 = rbinorm(1000,1,11,1,3,0.4);      xBiMod3 = xBiMod3[xBiMod3>=-2 & xBiMod3<=20];

    fittmop3 = PolynomialFit(xBiMod3,maxDegree=grados[j],-2,20,1,2);

    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fBiMod11113,fittmop3,-2,20);
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

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (BIMODAL(1,1)(11,3))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(-2,3.5,length.out=100);            x2 = seq(3.5,9,length.out=100);
x3 = seq(9,14.5,length.out=100);            x4 = seq(14.5,20,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xBiMod3,POTENTIAL_TYPE="MOP",evalRange=c(-2,20),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xBiMod3,maxDegree=7,-2,20,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fBiMod11113,c(-2,20),col='#0072B2',type='l',xlim = c(-2,20),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')







####################################################################

# MÉTODO SHENOY  (BIMODAL(1,1)(7,1))


# métodos con distintos parámetros
fitShenoy17 = ShenoyFit(fBiMod1171,7,-2,10,1);
fitShenoy14 = ShenoyFit(fBiMod1171,4,-2,10,1);
fitShenoy13 = ShenoyFit(fBiMod1171,3,-2,10,1);
fitShenoy27 = ShenoyFit(fBiMod1171,7,-2,10,2);
fitShenoy24 = ShenoyFit(fBiMod1171,4,-2,10,2);
fitShenoy23 = ShenoyFit(fBiMod1171,3,-2,10,2);
fitShenoy37 = ShenoyFit(fBiMod1171,7,-2,10,3);
fitShenoy34 = ShenoyFit(fBiMod1171,4,-2,10,3);
fitShenoy33 = ShenoyFit(fBiMod1171,3,-2,10,3);
fitShenoy47 = ShenoyFit(fBiMod1171,7,-2,10,4);
fitShenoy44 = ShenoyFit(fBiMod1171,4,-2,10,4);
fitShenoy43 = ShenoyFit(fBiMod1171,3,-2,10,4);
fitShenoy42 = ShenoyFit(fBiMod1171,2,-2,10,4);



divShenoy17 = DivergenceKL2(fBiMod1171,fitShenoy17[[1]],-2,10);
divShenoy14 = DivergenceKL2(fBiMod1171,fitShenoy14[[1]],-2,10);
divShenoy13 = DivergenceKL2(fBiMod1171,fitShenoy13[[1]],-2,10);

sp2 = seq(-2,10,length.out=3);
divShenoy27 = DivergenceKL2(fBiMod1171,fitShenoy27[[1]],sp2[1],sp2[2]) +
              DivergenceKL2(fBiMod1171,fitShenoy27[[2]],sp2[2],sp2[3]);
divShenoy24 = DivergenceKL2(fBiMod1171,fitShenoy24[[1]],sp2[1],sp2[2]) +
              DivergenceKL2(fBiMod1171,fitShenoy24[[2]],sp2[2],sp2[3]);
divShenoy23 = DivergenceKL2(fBiMod1171,fitShenoy23[[1]],sp2[1],sp2[2]) +
              DivergenceKL2(fBiMod1171,fitShenoy23[[2]],sp2[2],sp2[3]);

sp3 = seq(-2,10,length.out=4);
divShenoy37 = DivergenceKL2(fBiMod1171,fitShenoy37[[1]],sp3[1],sp3[2]) +
              DivergenceKL2(fBiMod1171,fitShenoy37[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBiMod1171,fitShenoy37[[3]],sp3[3],sp3[4]);
divShenoy34 = DivergenceKL2(fBiMod1171,fitShenoy34[[1]],sp3[1],sp3[2]) +
              DivergenceKL2(fBiMod1171,fitShenoy34[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBiMod1171,fitShenoy34[[3]],sp3[3],sp3[4]);
divShenoy33 = DivergenceKL2(fBiMod1171,fitShenoy33[[1]],sp3[1],sp3[2]) +
              DivergenceKL2(fBiMod1171,fitShenoy33[[2]],sp3[2],sp3[3]) +
              DivergenceKL2(fBiMod1171,fitShenoy33[[3]],sp3[3],sp3[4]);

sp4 = seq(-2,10,length.out=5);
divShenoy47 = DivergenceKL2(fBiMod1171,fitShenoy47[[1]],sp4[1],sp4[2]) +
              DivergenceKL2(fBiMod1171,fitShenoy47[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBiMod1171,fitShenoy47[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBiMod1171,fitShenoy47[[4]],sp4[4],sp4[5]);
divShenoy44 = DivergenceKL2(fBiMod1171,fitShenoy44[[1]],sp4[1],sp4[2]) +
              DivergenceKL2(fBiMod1171,fitShenoy44[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBiMod1171,fitShenoy44[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBiMod1171,fitShenoy44[[4]],sp4[4],sp4[5]);
divShenoy43 = DivergenceKL2(fBiMod1171,fitShenoy43[[1]],sp4[1],sp4[2]) +
              DivergenceKL2(fBiMod1171,fitShenoy43[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBiMod1171,fitShenoy43[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBiMod1171,fitShenoy43[[4]],sp4[4],sp4[5]);
divShenoy42 = DivergenceKL2(fBiMod1171,fitShenoy42[[1]],sp4[1],sp4[2]) +
              DivergenceKL2(fBiMod1171,fitShenoy42[[2]],sp4[2],sp4[3]) +
              DivergenceKL2(fBiMod1171,fitShenoy42[[3]],sp4[3],sp4[4]) +
              DivergenceKL2(fBiMod1171,fitShenoy42[[4]],sp4[4],sp4[5]);





####################################################################

# MÉTODO PÉREZ-BERNABÉ  (BIMODAL(1,1)(7,1))
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
    xBiMod1 = rbinorm(50,1,7,1,1,0.5);         xBiMod1 = xBiMod1[xBiMod1>=-2 & xBiMod1<=10];

    fitPB1 = univMoTBF(xBiMod1,POTENTIAL_TYPE="MOP",evalRange=c(-2,10),maxParam=p[j]);
    fitPBpoly1 = coef(fitPB1);

    divPB1[j] = divPB1[j] + DivergenceKL(fBiMod1171,fitPBpoly1,-2,10);
  }
}


# 100
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBiMod2 = rbinorm(100,1,7,1,1,0.5);      xBiMod2 = xBiMod2[xBiMod2>=-2 & xBiMod2<=10];

    fitPB2 = univMoTBF(xBiMod2,POTENTIAL_TYPE="MOP",evalRange=c(-2,10),maxParam=p[j]);
    fitPBpoly2 = coef(fitPB2);

    divPB2[j] = divPB2[j] + DivergenceKL(fBiMod1171,fitPBpoly2,-2,10);
  }
}


# 1000
for (j in 1:length(p))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBiMod3 = rbinorm(1000,1,7,1,1,0.5);       xBiMod3 = xBiMod3[xBiMod3>=-2 & xBiMod3<=10];

    fitPB3 = univMoTBF(xBiMod3,POTENTIAL_TYPE="MOP",evalRange=c(-2,10),maxParam=p[j]);
    fitPBpoly3 = coef(fitPB3);

    divPB3[j] = divPB3[j] + DivergenceKL(fBiMod1171,fitPBpoly3,-2,10);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divPB1 = divPB1 / 10;       divPB2 = divPB2 / 10;       divPB3 = divPB3 / 10;






####################################################################

# MÉTODO tMoPs  (BIMODAL(1,1)(7,1))
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
    xBiMod1 = rbinorm(50,1,7,1,1,0.5);         xBiMod1 = xBiMod1[xBiMod1>=-2 & xBiMod1<=10];

    fittmop1 = PolynomialFit(xBiMod1,maxDegree=grados[j],-2,10,1,2);

    divtmop1[j] = divtmop1[j] + DivergenceKLtMoP(fBiMod1171,fittmop1,-2,10);
  }
}


# 100
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBiMod2 = rbinorm(100,1,7,1,1,0.5);       xBiMod2 = xBiMod2[xBiMod2>=-2 & xBiMod2<=10];

    fittmop2 = PolynomialFit(xBiMod2,maxDegree=grados[j],-2,10,1,2);

    divtmop2[j] = divtmop2[j] + DivergenceKLtMoP(fBiMod1171,fittmop2,-2,10);
  }
}



# 1000
for (j in 1:length(grados))
{
  for (i in 1:10)
  {
    set.seed(i+j)
    xBiMod3 = rbinorm(1000,1,7,1,1,0.5);      xBiMod3 = xBiMod3[xBiMod3>=-2 & xBiMod3<=10];

    fittmop3 = PolynomialFit(xBiMod3,maxDegree=grados[j],-2,10,1,2);

    divtmop3[j] = divtmop3[j] + DivergenceKLtMoP(fBiMod1171,fittmop3,-2,10);
  }
}


# calcular la media de las divergencias dividiendo entre el número de elementos
divtmop1 = divtmop1 / 10;       divtmop2 = divtmop2 / 10;       divtmop3 = divtmop3 / 10;




# divergencias en la columna de resultados
resultados$`Div-KL`[63:93] = c(divShenoy17,divShenoy14,divShenoy13,divShenoy27,divShenoy24,
                              divShenoy23,divShenoy37,divShenoy34,divShenoy33,divShenoy47,
                              divShenoy44,divShenoy43,divShenoy42,divPB1,divPB2,divPB3,
                              divtmop1,divtmop2,divtmop3)



resultadosBiMod = resultados;

save(resultadosBiMod,file='resultadosBiMod.RData');






####################################################################

# GRÁFICA CON EL MEJOR MODELO DE CADA MÉTODO  (BIMODAL(1,1)(7,1))
# Shenoy 4 intervalos y grado 8 (modelo 7)
# Pérez-Bernabé con tamaño 1000 (modelo 3)
# tMoP con tamaño 1000 (modelo 3)


x1 = seq(-2,1,length.out=100);            x2 = seq(1,4,length.out=100);
x3 = seq(4,7,length.out=100);             x4 = seq(7,10,length.out=100);
x = c(x1,x2,x3,x4)

# Shenoy
yShenoy1 = fitShenoy47[[1]](x1);         yShenoy2 = fitShenoy47[[2]](x2);
yShenoy3 = fitShenoy47[[3]](x3);         yShenoy4 = fitShenoy47[[4]](x4);
yShenoy = c(yShenoy1,yShenoy2,yShenoy3,yShenoy4);

# Pérez Bernabé
fitPB = univMoTBF(xBiMod3,POTENTIAL_TYPE="MOP",evalRange=c(-2,10),maxParam=8);
yPB = predict(as.polynomial(coef(fitPB)),x);

# tMoP
fittmop = PolynomialFit(xBiMod3,maxDegree=7,-2,10,1,2);
ytmop = predicttMoP(fittmop,x);
ytmop[1] = ytmop[2];

# distribución original
plot(fBiMod1171,c(-2,10),col='#0072B2',type='l',xlim = c(-2,10),lwd=4,lty=2)     # original
points(x,yShenoy,col='#D55E00',type='l',lwd=2)
points(x,yPB,col='#009E73',type='l',lwd=2)
points(x,ytmop,col='#CC79A7',type='l',lwd=2)
legend(4,1.05, legend=c("Original", "Shenoy","Pérez-B","tMoP"),
       fill=c('#0072B2','#D55E00','#009E73','#CC79A7'),cex=0.5,bty='n')



























