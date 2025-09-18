
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

load('tMoPsLogNorm01212.RData')

# cuidado al usar el comando integral con polinomios al cargar otros paquetes
# hay otros paquetes que hacen uso de esa función con prioridad sobre polynom


####################################################################

# dataframe para los resultados
# distribuciones en los intervalos [0,1]-[0.2,1.2] cada 0.1

distribuciones = c('LN(0.08,1.14)', 'LN(0.11,0.76)', 'LN(0.16,0.25)', 'LN(0.27,0.69)',
                   'LN(0.35,0.87)', 'LN(0.39,1.04)', 'LN(0.44,0.62)', 'LN(0.53,1.12)',
                   'LN(0.57,0.46)', 'LN(0.64,0.88)', 'LN(0.69,1.07)', 'LN(0.77,0.57)',
                   'LN(0.81,1.15)', 'LN(0.91,0.93)', 'LN(0.94,0.36)');
modelo = c('interpolar','tmop')

# matriz con cada una de las combinaciones en cada fila
resultados = expand.grid(distribuciones,modelo,c(0))
colnames(resultados) = c('Distribucion','Modelo','DivKL')

# ordenar la matriz de combinaciones por distribución
resultados = resultados[order(resultados$Distribucion),]
row.names(resultados) = 1:nrow(resultados)



####################################################################

# Funciones de densidad
fLogNorm1 = function(x)   1/(x*sqrt(2*pi*1.14)) * exp(-(log(x)-0.08)^2/(2*1.14))
fLogNorm2 = function(x)   1/(x*sqrt(2*pi*0.76)) * exp(-(log(x)-0.11)^2/(2*0.76))
fLogNorm3 = function(x)   1/(x*sqrt(2*pi*0.25)) * exp(-(log(x)-0.16)^2/(2*0.25))
fLogNorm4 = function(x)   1/(x*sqrt(2*pi*0.69)) * exp(-(log(x)-0.27)^2/(2*0.69))
fLogNorm5 = function(x)   1/(x*sqrt(2*pi*0.87)) * exp(-(log(x)-0.35)^2/(2*0.87))
fLogNorm6 = function(x)   1/(x*sqrt(2*pi*1.04)) * exp(-(log(x)-0.39)^2/(2*1.04))
fLogNorm7 = function(x)   1/(x*sqrt(2*pi*0.62)) * exp(-(log(x)-0.44)^2/(2*0.62))
fLogNorm8 = function(x)   1/(x*sqrt(2*pi*1.12)) * exp(-(log(x)-0.53)^2/(2*1.12))
fLogNorm9 = function(x)   1/(x*sqrt(2*pi*0.46)) * exp(-(log(x)-0.57)^2/(2*0.46))
fLogNorm10 = function(x)   1/(x*sqrt(2*pi*0.88)) * exp(-(log(x)-0.64)^2/(2*0.88))
fLogNorm11 = function(x)   1/(x*sqrt(2*pi*1.07)) * exp(-(log(x)-0.69)^2/(2*1.07))
fLogNorm12 = function(x)   1/(x*sqrt(2*pi*0.57)) * exp(-(log(x)-0.77)^2/(2*0.57))
fLogNorm13 = function(x)   1/(x*sqrt(2*pi*1.15)) * exp(-(log(x)-0.81)^2/(2*1.15))
fLogNorm14 = function(x)   1/(x*sqrt(2*pi*0.93)) * exp(-(log(x)-0.91)^2/(2*0.93))
fLogNorm15 = function(x)   1/(x*sqrt(2*pi*0.36)) * exp(-(log(x)-0.94)^2/(2*0.36))





####################################################################


## AJUSTE 1: LN(0.08,1.14)

set.seed(5)
xLogNorm1 = rlnorm(1000,0.08,1.14);

# ajustes
interpolar1 = Interpolate2param(0.08,1.14,tMoPsLogNorm01212,0,1,0.2,1.2,0.1)
fittMoP1 = PolynomialFit(xLogNorm1,maxDegree=7,0,15,1,2);

# divergencia
div1 = DivergenceKLtMoP(fLogNorm1,fittMoP1,0,15);
divInter1 = DivergenceKLtMoP(fLogNorm1,interpolar1,0,15);



####################################################################


## AJUSTE 2: LN(0.11,0.76)

set.seed(5)
xLogNorm2 = rlnorm(1000,0.11,0.76);

# ajustes
interpolar2 = Interpolate2param(0.11,0.76,tMoPsLogNorm01212,0,1,0.2,1.2,0.1)
fittMoP2 = PolynomialFit(xLogNorm2,maxDegree=7,0,15,1,2);

# divergencia
div2 = DivergenceKLtMoP(fLogNorm2,fittMoP2,0,15);
divInter2 = DivergenceKLtMoP(fLogNorm2,interpolar2,0,15);



####################################################################


## AJUSTE 3: LN(0.16,0.25)

set.seed(5)
xLogNorm3 = rlnorm(1000,0.16,0.25);

# ajustes
interpolar3 = Interpolate2param(0.16,0.25,tMoPsLogNorm01212,0,1,0.2,1.2,0.1)
fittMoP3 = PolynomialFit(xLogNorm3,maxDegree=7,0,15,1,2);

# divergencia
div3 = DivergenceKLtMoP(fLogNorm3,fittMoP3,0,15);
divInter3 = DivergenceKLtMoP(fLogNorm3,interpolar3,0,15);



####################################################################


## AJUSTE 4: LN(0.27,0.69)

set.seed(5)
xLogNorm4 = rlnorm(1000,0.27,0.69);

# ajustes
interpolar4 = Interpolate2param(0.27,0.69,tMoPsLogNorm01212,0,1,0.2,1.2,0.1)
fittMoP4 = PolynomialFit(xLogNorm4,maxDegree=7,0,15,1,2);

# divergencia
div4 = DivergenceKLtMoP(fLogNorm4,fittMoP4,0,15);
divInter4 = DivergenceKLtMoP(fLogNorm4,interpolar4,0,15);



####################################################################


## AJUSTE 5: LN(0.35,0.87)

set.seed(5)
xLogNorm5 = rlnorm(1000,0.35,0.87);

# ajustes
interpolar5 = Interpolate2param(0.35,0.87,tMoPsLogNorm01212,0,1,0.2,1.2,0.1)
fittMoP5 = PolynomialFit(xLogNorm5,maxDegree=7,0,15,1,2);

# divergencia
div5 = DivergenceKLtMoP(fLogNorm5,fittMoP5,0,15);
divInter5 = DivergenceKLtMoP(fLogNorm5,interpolar5,0,15);



####################################################################


## AJUSTE 6: LN(0.39,1.04)

set.seed(5)
xLogNorm6 = rlnorm(1000,0.39,1.04);

# ajustes
interpolar6 = Interpolate2param(0.39,1.04,tMoPsLogNorm01212,0,1,0.2,1.2,0.1)
fittMoP6 = PolynomialFit(xLogNorm6,maxDegree=7,0,15,1,2);

# divergencia
div6 = DivergenceKLtMoP(fLogNorm6,fittMoP6,0,15);
divInter6 = DivergenceKLtMoP(fLogNorm6,interpolar6,0,15);



####################################################################


## AJUSTE 7: LN(0.44,0.62)

set.seed(5)
xLogNorm7 = rlnorm(1000,0.44,0.62);

# ajustes
interpolar7 = Interpolate2param(0.44,0.62,tMoPsLogNorm01212,0,1,0.2,1.2,0.1)
fittMoP7 = PolynomialFit(xLogNorm7,maxDegree=7,0,15,1,2);

# divergencia
div7 = DivergenceKLtMoP(fLogNorm7,fittMoP7,0,15);
divInter7 = DivergenceKLtMoP(fLogNorm7,interpolar7,0,15);



####################################################################


## AJUSTE 8: LN(0.53,1.12)

set.seed(5)
xLogNorm8 = rlnorm(1000,0.53,1.12);

# ajustes
interpolar8 = Interpolate2param(0.53,1.12,tMoPsLogNorm01212,0,1,0.2,1.2,0.1)
fittMoP8 = PolynomialFit(xLogNorm8,maxDegree=7,0,15,1,2);

# divergencia
div8 = DivergenceKLtMoP(fLogNorm8,fittMoP8,0,15);
divInter8 = DivergenceKLtMoP(fLogNorm8,interpolar8,0,15);



####################################################################


## AJUSTE 9: LN(0.57,0.46)

set.seed(5)
xLogNorm9 = rlnorm(1000,0.57,0.46);

# ajustes
interpolar9 = Interpolate2param(0.57,0.46,tMoPsLogNorm01212,0,1,0.2,1.2,0.1)
fittMoP9 = PolynomialFit(xLogNorm9,maxDegree=7,0,15,1,2);

# divergencia
div9 = DivergenceKLtMoP(fLogNorm9,fittMoP9,0,15);
divInter9 = DivergenceKLtMoP(fLogNorm9,interpolar9,0,15);



####################################################################


## AJUSTE 10: LN(0.64,0.88)

set.seed(5)
xLogNorm10 = rlnorm(1000,0.64,0.88);

# ajustes
interpolar10 = Interpolate2param(0.64,0.88,tMoPsLogNorm01212,0,1,0.2,1.2,0.1)
fittMoP10 = PolynomialFit(xLogNorm10,maxDegree=7,0,15,1,2);

# divergencia
div10 = DivergenceKLtMoP(fLogNorm10,fittMoP10,0,15);
divInter10 = DivergenceKLtMoP(fLogNorm10,interpolar10,0,15);



####################################################################


## AJUSTE 11: LN(0.69,1.07)

set.seed(5)
xLogNorm11 = rlnorm(1000,0.69,1.07);

# ajustes
interpolar11 = Interpolate2param(0.69,1.07,tMoPsLogNorm01212,0,1,0.2,1.2,0.1)
fittMoP11 = PolynomialFit(xLogNorm11,maxDegree=7,0,15,1,2);

# divergencia
div11 = DivergenceKLtMoP(fLogNorm11,fittMoP11,0,15);
divInter11 = DivergenceKLtMoP(fLogNorm11,interpolar11,0,15);



####################################################################


## AJUSTE 12: LN(0.77,0.57)

set.seed(5)
xLogNorm12 = rlnorm(1000,0.77,0.57);

# ajustes
interpolar12 = Interpolate2param(0.77,0.57,tMoPsLogNorm01212,0,1,0.2,1.2,0.1)
fittMoP12 = PolynomialFit(xLogNorm12,maxDegree=7,0,15,1,2);

# divergencia
div12 = DivergenceKLtMoP(fLogNorm12,fittMoP12,0,15);
divInter12 = DivergenceKLtMoP(fLogNorm12,interpolar12,0,15);



####################################################################


## AJUSTE 13: LN(0.81,1.15)

set.seed(5)
xLogNorm13 = rlnorm(1000,0.81,1.15);

# ajustes
interpolar13 = Interpolate2param(0.81,1.15,tMoPsLogNorm01212,0,1,0.2,1.2,0.1)
fittMoP13 = PolynomialFit(xLogNorm13,maxDegree=7,0,15,1,2);

# divergencia
div13 = DivergenceKLtMoP(fLogNorm13,fittMoP13,0,15);
divInter13 = DivergenceKLtMoP(fLogNorm13,interpolar13,0,15);



####################################################################


## AJUSTE 14: LN(0.91,0.93)

set.seed(5)
xLogNorm14 = rlnorm(1000,0.91,0.93);

# ajustes
interpolar14 = Interpolate2param(0.91,0.93,tMoPsLogNorm01212,0,1,0.2,1.2,0.1)
fittMoP14 = PolynomialFit(xLogNorm14,maxDegree=7,0,15,1,2);

# divergencia
div14 = DivergenceKLtMoP(fLogNorm14,fittMoP14,0,15);
divInter14 = DivergenceKLtMoP(fLogNorm14,interpolar14,0,15);



####################################################################


## AJUSTE 15: LN(0.94,0.36)

set.seed(5)
xLogNorm15 = rlnorm(1000,0.94,0.36);

# ajustes
interpolar15 = Interpolate2param(0.94,0.36,tMoPsLogNorm01212,0,1,0.2,1.2,0.1)
fittMoP15 = PolynomialFit(xLogNorm15,maxDegree=7,0,15,1,2);

# divergencia
div15 = DivergenceKLtMoP(fLogNorm15,fittMoP15,0,15);
divInter15 = DivergenceKLtMoP(fLogNorm15,interpolar15,0,15);






####################################################################


# divergencias en la columna de resultados
resultados$`DivKL` = c(divInter1,div1,divInter2,div2,divInter3,div3,divInter4,div4,
                       divInter5,div5,divInter6,div6,divInter7,div7,divInter8,div8,
                       divInter9,div9,divInter10,div10,divInter11,div11,divInter12,div12,
                       divInter13,div13,divInter14,div14,divInter15,div15);

resultadosLogNorm = resultados;

save(resultadosLogNorm,file='ResultadosLogNorm.RData');





####################################################################

## GRÁFICAS

####################################################################


## AJUSTE 4: LN(0.27,0.69)

x = seq(0,15,length.out=1000);
ytMoP4 = predicttMoP(fittMoP4,x);
ytMoP4[1] = ytMoP4[2];
ytMoP4[1000] = ytMoP4[999];
ytMoPInter4 = predicttMoP(interpolar4,x);
ytMoPInter4[1] = ytMoPInter4[2];
ytMoPInter4[1000] = ytMoPInter4[999];

plot(fLogNorm4,c(0,15),col='#0072B2',type='l',xlim=c(0,15),lwd=4,lty=2,xlab='',ylab='')
points(x,ytMoPInter4,col='#D55E00',type='l',lwd=2)
points(x,ytMoP4,col='darkgreen',type='l',lwd=2)
legend(11,0.5, legend=c("Theoretical", "Interpolation","Direct fit"),seg.len=2,lwd=2,
       fill=c('#0072B2','#D55E00','darkgreen'),cex=0.5,bty='n',lty=c(2,1,1))



####################################################################


## AJUSTE 12: LN(0.77,0.57)

x = seq(0,15,length.out=1000);
ytMoP12 = predicttMoP(fittMoP12,x);
ytMoP12[1] = ytMoP12[2];
ytMoP12[1000] = ytMoP12[999];
ytMoPInter12 = predicttMoP(interpolar12,x);
ytMoPInter12[1] = ytMoPInter12[2];
ytMoPInter12[1000] = ytMoPInter12[999];


## GUARDAR EN PDF
pdf("GraficaMetaLogNorm077-057.pdf",width=7,height=5)

#par(mar = c(2, 2, 0.5, 0.5))
plot(fLogNorm12,c(0,15),col='#0072B2',type='l',xlim=c(0,15),lwd=4,lty=2,xlab='',ylab='')
points(x,ytMoPInter12,col='#D55E00',type='l',lwd=2)
points(x,ytMoP12,col='darkgreen',type='l',lwd=2)
legend(11,0.3, legend=c("Theoretical", "Interpolation","Direct tMoP fit"),seg.len=2,lwd=2,
       col=c('#0072B2','#D55E00','darkgreen'),cex=0.8,bty='n',lty=c(2,1,1))

dev.off()






