
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

load('tMoPsBeta0101.RData')
load('tMoPsBeta0112.RData')
load('tMoPsBeta1212.RData')

# cuidado al usar el comando integral con polinomios al cargar otros paquetes
# hay otros paquetes que hacen uso de esa función con prioridad sobre polynom



####################################################################

# dataframe para los resultados
# distribuciones en los intervalos [0,1]-[0,1]   [0,1]-[1,2]   [1,2]-[1,2] cada 0.1

distribuciones = c('Be(0.33,0.85)', 'Be(0.45,0.44)', 'Be(0.56,0.77)', 'Be(0.62,0.94)',
                   'Be(0.82,0.58)', 'Be(0.27,1.42)', 'Be(0.35,1.93)', 'Be(0.52,1.14)',
                   'Be(0.85,1.36)', 'Be(0.92,1.65)', 'Be(1.18,1.76)', 'Be(1.33,1.55)',
                   'Be(1.52,1.91)', 'Be(1.61,1.25)', 'Be(1.84,1.48)');
modelo = c('interpolar','tmop')

# matriz con cada una de las combinaciones en cada fila
resultados = expand.grid(distribuciones,modelo,c(0))
colnames(resultados) = c('Distribucion','Modelo','DivKL')

# ordenar la matriz de combinaciones por distribución
resultados = resultados[order(resultados$Distribucion),]
row.names(resultados) = 1:nrow(resultados)



####################################################################

# Funciones de densidad

fBeta1 = function(x)   x^(0.33-1)*(1-x)^(0.85-1)*gamma(0.33+0.85)/(gamma(0.33)*gamma(0.85));
fBeta2 = function(x)   x^(0.45-1)*(1-x)^(0.44-1)*gamma(0.45+0.44)/(gamma(0.45)*gamma(0.44));
fBeta3 = function(x)   x^(0.56-1)*(1-x)^(0.77-1)*gamma(0.56+0.77)/(gamma(0.56)*gamma(0.77));
fBeta4 = function(x)   x^(0.62-1)*(1-x)^(0.94-1)*gamma(0.62+0.94)/(gamma(0.62)*gamma(0.94));
fBeta5 = function(x)   x^(0.82-1)*(1-x)^(0.58-1)*gamma(0.82+0.58)/(gamma(0.82)*gamma(0.58));
fBeta6 = function(x)   x^(0.27-1)*(1-x)^(1.42-1)*gamma(0.27+1.42)/(gamma(0.27)*gamma(1.42));
fBeta7 = function(x)   x^(0.35-1)*(1-x)^(1.93-1)*gamma(0.35+1.93)/(gamma(0.35)*gamma(1.93));
fBeta8 = function(x)   x^(0.52-1)*(1-x)^(1.14-1)*gamma(0.52+1.14)/(gamma(0.52)*gamma(1.14));
fBeta9 = function(x)   x^(0.85-1)*(1-x)^(1.36-1)*gamma(0.85+1.36)/(gamma(0.85)*gamma(1.36));
fBeta10 = function(x)   x^(0.92-1)*(1-x)^(1.65-1)*gamma(0.92+1.65)/(gamma(0.92)*gamma(1.65));
fBeta11 = function(x)   x^(1.18-1)*(1-x)^(1.76-1)*gamma(1.18+1.76)/(gamma(1.18)*gamma(1.76));
fBeta12 = function(x)   x^(1.33-1)*(1-x)^(1.55-1)*gamma(1.33+1.55)/(gamma(1.33)*gamma(1.55));
fBeta13 = function(x)   x^(1.52-1)*(1-x)^(1.91-1)*gamma(1.52+1.91)/(gamma(1.52)*gamma(1.91));
fBeta14 = function(x)   x^(1.61-1)*(1-x)^(1.25-1)*gamma(1.61+1.25)/(gamma(1.61)*gamma(1.25));
fBeta15 = function(x)   x^(1.84-1)*(1-x)^(1.48-1)*gamma(1.84+1.48)/(gamma(1.84)*gamma(1.48));





####################################################################


## AJUSTE 1: Beta(0.33, 0.85)

set.seed(5)
xBeta1 = rbeta(1000,0.33,0.85);

# ajustes
interpolar1 = Interpolate2param(0.33,0.85,tMoPsBeta0101,0,1,0,1,0.1)
fittMoP1 = PolynomialFit(xBeta1,maxDegree=7,0,1,1,2);

# divergencia
div1 = DivergenceKLtMoP(fBeta1,fittMoP1,0,1);
divInter1 = DivergenceKLtMoP(fBeta1,interpolar1,0,1);



####################################################################


## AJUSTE 2: Beta(0.45, 0.44)

set.seed(5)
xBeta2 = rbeta(1000,0.45,0.44);

# ajustes
interpolar2 = Interpolate2param(0.45,0.44,tMoPsBeta0101,0,1,0,1,0.1)
fittMoP2 = PolynomialFit(xBeta2,maxDegree=7,0,1,1,2);

# divergencia
div2 = DivergenceKLtMoP(fBeta2,fittMoP2,0,1);
divInter2 = DivergenceKLtMoP(fBeta2,interpolar2,0,1);



####################################################################


## AJUSTE 3: Beta(0.56, 0.77)

set.seed(5)
xBeta3 = rbeta(1000,0.56,0.77);

# ajustes
interpolar3 = Interpolate2param(0.56,0.77,tMoPsBeta0101,0,1,0,1,0.1)
fittMoP3 = PolynomialFit(xBeta3,maxDegree=7,0,1,1,2);

# divergencia
div3 = DivergenceKLtMoP(fBeta3,fittMoP3,0,1);
divInter3 = DivergenceKLtMoP(fBeta3,interpolar3,0,1);



####################################################################


## AJUSTE 4: Beta(0.62, 0.94)

set.seed(5)
xBeta4 = rbeta(1000,0.62,0.94);

# ajustes
interpolar4 = Interpolate2param(0.62,0.94,tMoPsBeta0101,0,1,0,1,0.1)
fittMoP4 = PolynomialFit(xBeta4,maxDegree=7,0,1,1,2);

# divergencia
div4 = DivergenceKLtMoP(fBeta4,fittMoP4,0,1);
divInter4 = DivergenceKLtMoP(fBeta4,interpolar4,0,1);



####################################################################


## AJUSTE 5: Beta(0.82, 0.58)

set.seed(5)
xBeta5 = rbeta(1000,0.82,0.58);

# ajustes
interpolar5 = Interpolate2param(0.82,0.58,tMoPsBeta0101,0,1,0,1,0.1)
fittMoP5 = PolynomialFit(xBeta5,maxDegree=7,0,1,1,2);

# divergencia
div5 = DivergenceKLtMoP(fBeta5,fittMoP5,0,1);
divInter5 = DivergenceKLtMoP(fBeta5,interpolar5,0,1);



####################################################################


## AJUSTE 6: Beta(0.27, 1.42)

set.seed(5)
xBeta6 = rbeta(1000,0.27,1.42);

# ajustes
interpolar6 = Interpolate2param(0.27,1.42,tMoPsBeta0112,0,1,1,2,0.1)
fittMoP6 = PolynomialFit(xBeta6,maxDegree=7,0,1,1,2);

#divergencia
div6 = DivergenceKLtMoP(fBeta6,fittMoP6,0,1);
divInter6 = DivergenceKLtMoP(fBeta6,interpolar6,0,1);



####################################################################


## AJUSTE 7: Beta(0.35, 1.93)

set.seed(5)
xBeta7 = rbeta(1000,0.35,1.93);

# ajustes
interpolar7 = Interpolate2param(0.35,1.93,tMoPsBeta0112,0,1,1,2,0.1)
fittMoP7 = PolynomialFit(xBeta7,maxDegree=7,0,1,1,2);

#divergencia
div7 = DivergenceKLtMoP(fBeta7,fittMoP7,0,1);
divInter7 = DivergenceKLtMoP(fBeta7,interpolar7,0,1);



####################################################################


## AJUSTE 8: Beta(0.52, 1.14)

set.seed(5)
xBeta8 = rbeta(1000,0.52,1.14);

# ajustes
interpolar8 = Interpolate2param(0.52,1.14,tMoPsBeta0112,0,1,1,2,0.1)
fittMoP8 = PolynomialFit(xBeta8,maxDegree=7,0,1,1,2);

#divergencia
div8 = DivergenceKLtMoP(fBeta8,fittMoP8,0,1);
divInter8 = DivergenceKLtMoP(fBeta8,interpolar8,0,1);



####################################################################


## AJUSTE 9: Beta(0.85, 1.36)

set.seed(5)
xBeta9 = rbeta(1000,0.85,1.36);

# ajustes
interpolar9 = Interpolate2param(0.85,1.36,tMoPsBeta0112,0,1,1,2,0.1)
fittMoP9 = PolynomialFit(xBeta9,maxDegree=7,0,1,1,2);

#divergencia
div9 = DivergenceKLtMoP(fBeta9,fittMoP9,0,1);
divInter9 = DivergenceKLtMoP(fBeta9,interpolar9,0,1);



####################################################################


## AJUSTE 10: Beta(0.92, 1.65)

set.seed(5)
xBeta10 = rbeta(1000,0.92,1.65);

# ajustes
interpolar10 = Interpolate2param(0.92,1.65,tMoPsBeta0112,0,1,1,2,0.1)
fittMoP10 = PolynomialFit(xBeta10,maxDegree=7,0,1,1,2);

#divergencia
div10 = DivergenceKLtMoP(fBeta10,fittMoP10,0,1);
divInter10 = DivergenceKLtMoP(fBeta10,interpolar10,0,1);



####################################################################


## AJUSTE 11: Beta(1.18, 1.76)

set.seed(5)
xBeta11 = rbeta(1000,1.18,1.76);

# ajustes
interpolar11 = Interpolate2param(1.18,1.76,tMoPsBeta0112,1,2,1,2,0.1)
fittMoP11 = PolynomialFit(xBeta11,maxDegree=7,0,1,1,2);

#divergencia
div11 = DivergenceKLtMoP(fBeta11,fittMoP11,0,1);
divInter11 = DivergenceKLtMoP(fBeta11,interpolar11,0,1);



####################################################################


## AJUSTE 12: Beta(1.33, 1.55)

set.seed(5)
xBeta12 = rbeta(1000,1.33,1.55);

# ajustes
interpolar12 = Interpolate2param(1.33,1.55,tMoPsBeta0112,1,2,1,2,0.1)
fittMoP12 = PolynomialFit(xBeta12,maxDegree=7,0,1,1,2);

#divergencia
div12 = DivergenceKLtMoP(fBeta12,fittMoP12,0,1);
divInter12 = DivergenceKLtMoP(fBeta12,interpolar12,0,1);



####################################################################


## AJUSTE 13: Beta(1.52, 1.91)

set.seed(5)
xBeta13 = rbeta(1000,1.52,1.91);

# ajustes
interpolar13 = Interpolate2param(1.52,1.91,tMoPsBeta0112,1,2,1,2,0.1)
fittMoP13 = PolynomialFit(xBeta13,maxDegree=7,0,1,1,2);

#divergencia
div13 = DivergenceKLtMoP(fBeta13,fittMoP13,0,1);
divInter13 = DivergenceKLtMoP(fBeta13,interpolar13,0,1);



####################################################################


## AJUSTE 14: Beta(1.61, 1.25)

set.seed(5)
xBeta14 = rbeta(1000,1.61,1.25);

# ajustes
interpolar14 = Interpolate2param(1.61,1.25,tMoPsBeta0112,1,2,1,2,0.1)
fittMoP14 = PolynomialFit(xBeta14,maxDegree=7,0,1,1,2);

#divergencia
div14 = DivergenceKLtMoP(fBeta14,fittMoP14,0,1);
divInter14 = DivergenceKLtMoP(fBeta14,interpolar14,0,1);



####################################################################


## AJUSTE 15: Beta(1.84, 1.48)

set.seed(5)
xBeta15 = rbeta(1000,1.84,1.48);

# ajustes
interpolar15 = Interpolate2param(1.84,1.48,tMoPsBeta0112,1,2,1,2,0.1)
fittMoP15 = PolynomialFit(xBeta15,maxDegree=7,0,1,1,2);

#divergencia
div15 = DivergenceKLtMoP(fBeta15,fittMoP15,0,1);
divInter15 = DivergenceKLtMoP(fBeta15,interpolar15,0,1);






####################################################################


# divergencias en la columna de resultados
resultados$`DivKL` = c(divInter1,div1,divInter2,div2,divInter3,div3,divInter4,div4,
                       divInter5,div5,divInter6,div6,divInter7,div7,divInter8,div8,
                       divInter9,div9,divInter10,div10,divInter11,div11,divInter12,div12,
                       divInter13,div13,divInter14,div14,divInter15,div15);

resultadosBeta = resultados;

save(resultadosBeta,file='ResultadosBeta.RData');





####################################################################

## GRÁFICAS

####################################################################


## AJUSTE 3: Beta(0.56, 0.77)

x = seq(0,1,length.out=1000);
ytMoP3 = predicttMoP(fittMoP3,x);
ytMoP3[1] = ytMoP3[2];
ytMoP3[1000] = ytMoP3[999];
ytMoPInter3 = predicttMoP(interpolar3,x);
ytMoPInter3[1] = ytMoPInter3[2];
ytMoPInter3[1000] = ytMoPInter3[999];

plot(fBeta3,c(0,1),col='#0072B2',type='l',xlim=c(0,1),lwd=4,lty=2,xlab='',ylab='') # original
points(x,ytMoPInter3,col='#D55E00',type='l',lwd=2)
points(x,ytMoP3,col='darkgreen',type='l',lwd=2)
legend(0.75,3.6, legend=c("Original", "Interpolar","tMoP"),
       fill=c('#0072B2','#D55E00','darkgreen'),cex=0.5,bty='n')


####################################################################


## AJUSTE 12: Beta(1.33, 1.55)

x = seq(0,1,length.out=1000);
ytMoP12 = predicttMoP(fittMoP12,x);
ytMoP12[1] = ytMoP12[2];
ytMoP12[1000] = ytMoP12[999];
ytMoPInter12 = predicttMoP(interpolar12,x);
ytMoPInter12[1] = ytMoPInter12[2];

plot(fBeta12,c(0,1),col='#0072B2',type='l',xlim = c(0,1),lwd=4,lty=2,xlab='',ylab='')
points(x,ytMoPInter12,col='#D55E00',type='l',lwd=2)
points(x,ytMoP12,col='darkgreen',type='l',lwd=2)
legend(0.75,1.58, legend=c("Original", "Interpolar","tMoP"),
       fill=c('#0072B2','#D55E00','darkgreen'),cex=0.5,bty='n')



####################################################################


## AJUSTE 9: Beta(0.85, 1.36)

x = seq(0,1,length.out=1000);
ytMoP9 = predicttMoP(fittMoP9,x);
ytMoP9[1] = ytMoP9[2];
ytMoP9[1000] = ytMoP9[999];
ytMoPInter9 = predicttMoP(interpolar9,x);
ytMoPInter9[1] = ytMoPInter9[2];
ytMoPInter9[1000] = ytMoPInter9[999];


## GUARDAR EN PDF
pdf("GraficaMetaBeta085-136.pdf",width=7,height=5)

par(mar = c(3,3,3,3))
plot(fBeta9,c(0,1),col='#0072B2',type='l',xlim = c(0,1),lwd=4,lty=2,xlab='',ylab='')
points(x,ytMoPInter9,col='#D55E00',type='l',lwd=2)
points(x,ytMoP9,col='darkgreen',type='l',lwd=2)
legend(0.75,2.1, legend=c("Original", "Interpolation","tMoP fit"),
       fill=c('#0072B2','#D55E00','darkgreen'),cex=0.8,bty='n')

dev.off()








