
# CARGAR DEPENDENCIAS AL PRINCIPIO (polinomios)
library(polynom)


# ponemos el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# cambiamos el directorio para cargar las funciones
setwd("../../../core")
source("PolynomialFit.R");      source("BIC.R");       source("Prediction.R");
source("Simulation.R");

# poner el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# cambiamos el directorio para cargar las funciones
setwd("../Univariate densities")
source("DivergenceKL.R");

# poner el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))




####################################################################

source("Interpolation.R")
source("DivergenceKL.R")

load('tMoPstStu16.RData')

# cuidado al usar el comando integral con polinomios al cargar otros paquetes
# hay otros paquetes que hacen uso de esa función con prioridad sobre polynom



####################################################################

# dataframe para los resultados
# distribuciones con parámetro en el intervalo [1,6] cada 0.05 (ajuste en [-20,20])

distribuciones = c('tStu(1.18)', 'tStu(1.42)', 'tStu(1.69)', 'tStu(2.33)', 'tStu(2.77)',
                   'tStu(2.92)', 'tStu(3.24)', 'tStu(3.63)', 'tStu(3.71)', 'tStu(4.07)',
                   'tStu(4.43)', 'tStu(4.82)', 'tStu(5.12)', 'tStu(5.58)', 'tStu(5.94)');
modelo = c('interpolar','tmop')

# matriz con cada una de las combinaciones en cada fila
resultados = expand.grid(distribuciones,modelo,c(0))
colnames(resultados) = c('Distribucion','Modelo','DivKL')

# ordenar la matriz de combinaciones por distribución
resultados = resultados[order(resultados$Distribucion),]
row.names(resultados) = 1:nrow(resultados)




####################################################################

# Funciones de densidad

ftS1 = function(x)  gamma((1.18+1)/2)/(sqrt(1.18*pi)*gamma(1.18/2))*(1+x^2/1.18)^(-(1.18+1)/2)
itS1 = pracma::integral(ftS1,-20,20);
ftStudent1 = function(x)
  (gamma((1.18+1)/2) / (sqrt(1.18*pi)*gamma(1.18/2)) * (1+x^2/1.18) ^ (-(1.18+1)/2)) / itS1

ftS2 = function(x)  gamma((1.42+1)/2)/(sqrt(1.42*pi)*gamma(1.42/2))*(1+x^2/1.42)^(-(1.42+1)/2)
itS2 = pracma::integral(ftS2,-20,20);
ftStudent2 = function(x)
  (gamma((1.42+1)/2) / (sqrt(1.42*pi)*gamma(1.42/2)) * (1+x^2/1.42) ^ (-(1.42+1)/2)) / itS2

ftS3 = function(x)  gamma((1.69+1)/2)/(sqrt(1.69*pi)*gamma(1.69/2))*(1+x^2/1.69)^(-(1.69+1)/2)
itS3 = pracma::integral(ftS3,-20,20);
ftStudent3 = function(x)
  (gamma((1.69+1)/2) / (sqrt(1.69*pi)*gamma(1.69/2)) * (1+x^2/1.69) ^ (-(1.69+1)/2)) / itS3

ftStudent4 = function(x)
  gamma((2.33+1)/2) / (sqrt(2.33*pi) * gamma(2.33/2)) * (1 + x^2 / 2.33) ^ (-(2.33+1) / 2)
ftStudent5 = function(x)
  gamma((2.77+1)/2) / (sqrt(2.77*pi) * gamma(2.77/2)) * (1 + x^2 / 2.77) ^ (-(2.77+1) / 2)
ftStudent6 = function(x)
  gamma((2.92+1)/2) / (sqrt(2.92*pi) * gamma(2.92/2)) * (1 + x^2 / 2.92) ^ (-(2.92+1) / 2)
ftStudent7 = function(x)
  gamma((3.24+1)/2) / (sqrt(3.24*pi) * gamma(3.24/2)) * (1 + x^2 / 3.24) ^ (-(3.24+1) / 2)
ftStudent8 = function(x)
  gamma((3.63+1)/2) / (sqrt(3.63*pi) * gamma(3.63/2)) * (1 + x^2 / 3.63) ^ (-(3.63+1) / 2)
ftStudent9 = function(x)
  gamma((3.71+1)/2) / (sqrt(3.71*pi) * gamma(3.71/2)) * (1 + x^2 / 3.71) ^ (-(3.71+1) / 2)
ftStudent10 = function(x)
  gamma((4.07+1)/2) / (sqrt(4.07*pi) * gamma(4.07/2)) * (1 + x^2 / 4.07) ^ (-(4.07+1) / 2)
ftStudent11 = function(x)
  gamma((4.43+1)/2) / (sqrt(4.43*pi) * gamma(4.43/2)) * (1 + x^2 / 4.43) ^ (-(4.43+1) / 2)
ftStudent12 = function(x)
  gamma((4.82+1)/2) / (sqrt(4.82*pi) * gamma(4.82/2)) * (1 + x^2 / 4.82) ^ (-(4.82+1) / 2)
ftStudent13 = function(x)
  gamma((5.12+1)/2) / (sqrt(5.12*pi) * gamma(5.12/2)) * (1 + x^2 / 5.12) ^ (-(5.12+1) / 2)
ftStudent14 = function(x)
  gamma((5.58+1)/2) / (sqrt(5.58*pi) * gamma(5.58/2)) * (1 + x^2 / 5.58) ^ (-(5.58+1) / 2)
ftStudent15 = function(x)
  gamma((5.94+1)/2) / (sqrt(5.94*pi) * gamma(5.94/2)) * (1 + x^2 / 5.94) ^ (-(5.94+1) / 2)



####################################################################


## AJUSTE 1: tStu(1.18)

set.seed(5)
xtStu1 = rt(1000,1.18);            xtStu1 = xtStu1[xtStu1>=-20 & xtStu1<=20];

# ajustes
interpolar1 = Interpolate1param(1.18,tMoPstStu16,1,6,0.05)
fittMoP1 = PolynomialFit(xtStu1,maxDegree=7,-20,20,1,2);

# divergencia
div1 = DivergenceKLtMoP(ftStudent1,fittMoP1,-20,20);
divInter1 = DivergenceKLtMoP(ftStudent1,interpolar1,-20,20);



####################################################################


## AJUSTE 2: tStu(1.42)

set.seed(5)
xtStu2 = rt(1000,1.42);            xtStu2 = xtStu2[xtStu2>=-20 & xtStu2<=20];

# ajustes
interpolar2 = Interpolate1param(1.42,tMoPstStu16,1,6,0.05)
fittMoP2 = PolynomialFit(xtStu2,maxDegree=7,-20,20,1,2);

# divergencia
div2 = DivergenceKLtMoP(ftStudent2,fittMoP2,-20,20);
divInter2 = DivergenceKLtMoP(ftStudent2,interpolar2,-20,20);



####################################################################


## AJUSTE 3: tStu(1.69)

set.seed(5)
xtStu3 = rt(1000,1.69);            xtStu3 = xtStu3[xtStu3>=-20 & xtStu3<=20];

# ajustes
interpolar3 = Interpolate1param(1.69,tMoPstStu16,1,6,0.05)
fittMoP3 = PolynomialFit(xtStu3,maxDegree=7,-20,20,1,2);

# divergencia
div3 = DivergenceKLtMoP(ftStudent3,fittMoP3,-20,20);
divInter3 = DivergenceKLtMoP(ftStudent3,interpolar3,-20,20);



####################################################################


## AJUSTE 4: tStu(2.33)

set.seed(5)
xtStu4 = rt(1000,2.33);            xtStu4 = xtStu4[xtStu4>=-20 & xtStu4<=20];

# ajustes
interpolar4 = Interpolate1param(2.33,tMoPstStu16,1,6,0.05)
fittMoP4 = PolynomialFit(xtStu4,maxDegree=7,-20,20,1,2);

# divergencia
div4 = DivergenceKLtMoP(ftStudent4,fittMoP4,-20,20);
divInter4 = DivergenceKLtMoP(ftStudent4,interpolar4,-20,20);



####################################################################


## AJUSTE 5: tStu(2.77)

set.seed(5)
xtStu5 = rt(1000,2.77);            xtStu5 = xtStu5[xtStu5>=-20 & xtStu5<=20];

# ajustes
interpolar5 = Interpolate1param(2.77,tMoPstStu16,1,6,0.05)
fittMoP5 = PolynomialFit(xtStu5,maxDegree=7,-20,20,1,2);

# divergencia
div5 = DivergenceKLtMoP(ftStudent5,fittMoP5,-20,20);
divInter5 = DivergenceKLtMoP(ftStudent5,interpolar5,-20,20);



####################################################################


## AJUSTE 6: tStu(2.92)

set.seed(5)
xtStu6 = rt(1000,2.92);            xtStu6 = xtStu6[xtStu6>=-20 & xtStu6<=20];

# ajustes
interpolar6 = Interpolate1param(2.92,tMoPstStu16,1,6,0.05)
fittMoP6 = PolynomialFit(xtStu6,maxDegree=7,-20,20,1,2);

# divergencia
div6 = DivergenceKLtMoP(ftStudent6,fittMoP6,-20,20);
divInter6 = DivergenceKLtMoP(ftStudent6,interpolar6,-20,20);



####################################################################


## AJUSTE 7: tStu(3.24)

set.seed(5)
xtStu7 = rt(1000,3.24);            xtStu7 = xtStu7[xtStu7>=-20 & xtStu7<=20];

# ajustes
interpolar7 = Interpolate1param(3.24,tMoPstStu16,1,6,0.05)
fittMoP7 = PolynomialFit(xtStu7,maxDegree=7,-20,20,1,2);

# divergencia
div7 = DivergenceKLtMoP(ftStudent7,fittMoP7,-20,20);
divInter7 = DivergenceKLtMoP(ftStudent7,interpolar7,-20,20);



####################################################################


## AJUSTE 8: tStu(3.63)

set.seed(5)
xtStu8 = rt(1000,3.63);            xtStu8 = xtStu8[xtStu8>=-20 & xtStu8<=20];

# ajustes
interpolar8 = Interpolate1param(3.63,tMoPstStu16,1,6,0.05)
fittMoP8 = PolynomialFit(xtStu8,maxDegree=7,-20,20,1,2);

# divergencia
div8 = DivergenceKLtMoP(ftStudent8,fittMoP8,-20,20);
divInter8 = DivergenceKLtMoP(ftStudent8,interpolar8,-20,20);



####################################################################


## AJUSTE 9: tStu(3.71)

set.seed(5)
xtStu9 = rt(1000,3.71);            xtStu9 = xtStu9[xtStu9>=-20 & xtStu9<=20];

# ajustes
interpolar9 = Interpolate1param(3.71,tMoPstStu16,1,6,0.05)
fittMoP9 = PolynomialFit(xtStu9,maxDegree=7,-20,20,1,2);

# divergencia
div9 = DivergenceKLtMoP(ftStudent9,fittMoP9,-20,20);
divInter9 = DivergenceKLtMoP(ftStudent9,interpolar9,-20,20);



####################################################################


## AJUSTE 10: tStu(4.07)

set.seed(5)
xtStu10 = rt(1000,4.07);            xtStu10 = xtStu10[xtStu10>=-20 & xtStu10<=20];

# ajustes
interpolar10 = Interpolate1param(4.07,tMoPstStu16,1,6,0.05)
fittMoP10 = PolynomialFit(xtStu10,maxDegree=7,-20,20,1,2);

# divergencia
div10 = DivergenceKLtMoP(ftStudent10,fittMoP10,-20,20);
divInter10 = DivergenceKLtMoP(ftStudent10,interpolar10,-20,20);



####################################################################


## AJUSTE 11: tStu(4.43)

set.seed(5)
xtStu11 = rt(1000,4.43);            xtStu11 = xtStu11[xtStu11>=-20 & xtStu11<=20];

# ajustes
interpolar11 = Interpolate1param(4.43,tMoPstStu16,1,6,0.05)
fittMoP11 = PolynomialFit(xtStu11,maxDegree=7,-20,20,1,2);

# divergencia
div11 = DivergenceKLtMoP(ftStudent11,fittMoP11,-20,20);
divInter11 = DivergenceKLtMoP(ftStudent11,interpolar11,-20,20);



####################################################################


## AJUSTE 12: tStu(4.82)

set.seed(5)
xtStu12 = rt(1000,4.82);            xtStu12 = xtStu12[xtStu12>=-20 & xtStu12<=20];

# ajustes
interpolar12 = Interpolate1param(4.82,tMoPstStu16,1,6,0.05)
fittMoP12 = PolynomialFit(xtStu12,maxDegree=7,-20,20,1,2);

# divergencia
div12 = DivergenceKLtMoP(ftStudent12,fittMoP12,-20,20);
divInter12 = DivergenceKLtMoP(ftStudent12,interpolar12,-20,20);



####################################################################


## AJUSTE 13: tStu(5.12)

set.seed(5)
xtStu13 = rt(1000,5.12);            xtStu13 = xtStu13[xtStu13>=-20 & xtStu13<=20];

# ajustes
interpolar13 = Interpolate1param(5.12,tMoPstStu16,1,6,0.05)
fittMoP13 = PolynomialFit(xtStu13,maxDegree=7,-20,20,1,2);

# divergencia
div13 = DivergenceKLtMoP(ftStudent13,fittMoP13,-20,20);
divInter13 = DivergenceKLtMoP(ftStudent13,interpolar13,-20,20);



####################################################################


## AJUSTE 14: tStu(5.58)

set.seed(5)
xtStu14 = rt(1000,5.58);            xtStu14 = xtStu14[xtStu14>=-20 & xtStu14<=20];

# ajustes
interpolar14 = Interpolate1param(5.58,tMoPstStu16,1,6,0.05)
fittMoP14 = PolynomialFit(xtStu14,maxDegree=7,-20,20,1,2);

# divergencia
div14 = DivergenceKLtMoP(ftStudent14,fittMoP14,-20,20);
divInter14 = DivergenceKLtMoP(ftStudent14,interpolar14,-20,20);



####################################################################


## AJUSTE 15: tStu(5.94)

set.seed(5)
xtStu15 = rt(1000,5.94);            xtStu15 = xtStu15[xtStu15>=-20 & xtStu15<=20];

# ajustes
interpolar15 = Interpolate1param(5.94,tMoPstStu16,1,6,0.05)
fittMoP15 = PolynomialFit(xtStu15,maxDegree=7,-20,20,1,2);

# divergencia
div15 = DivergenceKLtMoP(ftStudent15,fittMoP15,-20,20);
divInter15 = DivergenceKLtMoP(ftStudent15,interpolar15,-20,20);




####################################################################


# divergencias en la columna de resultados
resultados$`DivKL` = c(divInter1,div1,divInter2,div2,divInter3,div3,divInter4,div4,
                       divInter5,div5,divInter6,div6,divInter7,div7,divInter8,div8,
                       divInter9,div9,divInter10,div10,divInter11,div11,divInter12,div12,
                       divInter13,div13,divInter14,div14,divInter15,div15);

resultadostStudent = resultados;

save(resultadostStudent,file='ResultadostStudent.RData');






####################################################################

# GRÁFICAS

####################################################################


## AJUSTE 8: tStu(3.63)

x = seq(-20,20,length.out=1000);
ytMoP8 = predicttMoP(fittMoP8,x);
ytMoP8[1] = ytMoP8[2];
ytMoPInter8 = predicttMoP(interpolar8,x);
ytMoPInter8[1] = ytMoPInter8[2];


## GUARDAR EN PDF
pdf("GraficaMetatStu363.pdf",width=7,height=5)

par(mar = c(3,3,3,3))
plot(ftStudent8,c(-20,20),col='#0072B2',type='l',xlim = c(-20,20),lwd=4,lty=2)     # original
points(x,ytMoPInter8,col='#D55E00',type='l',lwd=2)
points(x,ytMoP8,col='darkgreen',type='l',lwd=2)
legend(13,0.33, legend=c("Original", "Interpolation","tMoP"),
       fill=c('#0072B2','#D55E00','darkgreen'),cex=0.8,bty='n')

dev.off()




####################################################################


## AJUSTE 13: tStu(5.12)

x = seq(-20,20,length.out=1000);
ytMoP13 = predicttMoP(fittMoP13,x);
ytMoP13[1] = ytMoP13[2];
ytMoPInter13 = predicttMoP(interpolar13,x);
ytMoPInter13[1] = ytMoPInter13[2];


## GUARDAR EN PDF
pdf("GraficaMetatStu512.pdf",width=7,height=5)

par(mar = c(3,3,3,3))
plot(ftStudent13,c(-20,20),col='#0072B2',type='l',xlim = c(-20,20),lwd=4,lty=2)     # original
points(x,ytMoPInter13,col='#D55E00',type='l',lwd=2)
points(x,ytMoP13,col='darkgreen',type='l',lwd=2)
legend(13,0.33, legend=c("Original", "Interpolation","tMoP"),
       fill=c('#0072B2','#D55E00','darkgreen'),cex=0.8,bty='n')

dev.off()








