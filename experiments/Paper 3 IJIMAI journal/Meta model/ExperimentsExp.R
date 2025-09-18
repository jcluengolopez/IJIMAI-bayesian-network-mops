
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

load('tMoPsExp16.RData')

# cuidado al usar el comando integral con polinomios al cargar otros paquetes
# hay otros paquetes que hacen uso de esa función con prioridad sobre polynom



####################################################################

# dataframe para los resultados
# distribuciones con parámetro en el intervalo [1,6] cada 0.05 (ajuste en [0,5])

distribuciones = c('Exp(1.33)', 'Exp(1.41)', 'Exp(1.67)', 'Exp(2.12)', 'Exp(2.44)',
                   'Exp(2.93)', 'Exp(3.08)', 'Exp(3.51)', 'Exp(3.72)', 'Exp(4.11)',
                   'Exp(4.52)', 'Exp(4.96)', 'Exp(5.28)', 'Exp(5.47)', 'Exp(5.81)');
modelo = c('interpolar','tmop')

# matriz con cada una de las combinaciones en cada fila
resultados = expand.grid(distribuciones,modelo,c(0))
colnames(resultados) = c('Distribucion','Modelo','DivKL')

# ordenar la matriz de combinaciones por distribución
resultados = resultados[order(resultados$Distribucion),]
row.names(resultados) = 1:nrow(resultados)



####################################################################

# Funciones de densidad

fExp1 = function(x)  1.33*exp(-1.33*x)
fExp2 = function(x)  1.41*exp(-1.41*x)
fExp3 = function(x)  1.67*exp(-1.67*x)
fExp4 = function(x)  2.12*exp(-2.12*x)
fExp5 = function(x)  2.44*exp(-2.44*x)
fExp6 = function(x)  2.93*exp(-2.93*x)
fExp7 = function(x)  3.08*exp(-3.08*x)
fExp8 = function(x)  3.51*exp(-3.51*x)
fExp9 = function(x)  3.72*exp(-3.72*x)
fExp10 = function(x)  4.11*exp(-4.11*x)
fExp11 = function(x)  4.52*exp(-4.52*x)
fExp12 = function(x)  4.96*exp(-4.96*x)
fExp13 = function(x)  5.28*exp(-5.28*x)
fExp14 = function(x)  5.47*exp(-5.47*x)
fExp15 = function(x)  5.81*exp(-5.81*x)


####################################################################


## AJUSTE 1: Exp(1.33)

set.seed(5)
xExp1 = rexp(1000,1.33);

# ajustes
interpolar1 = Interpolate1param(1.33,tMoPsExp16,1,6,0.05)
fittMoP1 = PolynomialFit(xExp1,maxDegree=7,0,5,1,2);

# divergencia
div1 = DivergenceKLtMoP(fExp1,fittMoP1,0,5);
divInter1 = DivergenceKLtMoP(fExp1,interpolar1,0,5);



####################################################################


## AJUSTE 2: Exp(1.41)

set.seed(5)
xExp2 = rexp(1000,1.41);

# ajustes
interpolar2 = Interpolate1param(1.41,tMoPsExp16,1,6,0.05)
fittMoP2 = PolynomialFit(xExp2,maxDegree=7,0,5,1,2);

# divergencia
div2 = DivergenceKLtMoP(fExp2,fittMoP2,0,5);
divInter2 = DivergenceKLtMoP(fExp2,interpolar2,0,5);



####################################################################


## AJUSTE 3: Exp(1.67)

set.seed(5)
xExp3 = rexp(1000,1.67);

# ajustes
interpolar3 = Interpolate1param(1.67,tMoPsExp16,1,6,0.05)
fittMoP3 = PolynomialFit(xExp3,maxDegree=7,0,5,1,2);

# divergencia
div3 = DivergenceKLtMoP(fExp3,fittMoP3,0,5);
divInter3 = DivergenceKLtMoP(fExp3,interpolar3,0,5);



####################################################################


## AJUSTE 4: Exp(2.12)

set.seed(5)
xExp4 = rexp(1000,2.12);

# ajustes
interpolar4 = Interpolate1param(2.12,tMoPsExp16,1,6,0.05)
fittMoP4 = PolynomialFit(xExp4,maxDegree=7,0,5,1,2);

# divergencia
div4 = DivergenceKLtMoP(fExp4,fittMoP4,0,5);
divInter4 = DivergenceKLtMoP(fExp4,interpolar4,0,5);



####################################################################


## AJUSTE 5: Exp(2.44)

set.seed(5)
xExp5 = rexp(1000,2.44);

# ajustes
interpolar5 = Interpolate1param(2.44,tMoPsExp16,1,6,0.05)
fittMoP5 = PolynomialFit(xExp5,maxDegree=7,0,5,1,2);

# divergencia
div5 = DivergenceKLtMoP(fExp5,fittMoP5,0,5);
divInter5 = DivergenceKLtMoP(fExp5,interpolar5,0,5);



####################################################################


## AJUSTE 6: Exp(2.93)

set.seed(5)
xExp6 = rexp(1000,2.93);

# ajustes
interpolar6 = Interpolate1param(2.93,tMoPsExp16,1,6,0.05)
fittMoP6 = PolynomialFit(xExp6,maxDegree=7,0,5,1,2);

# divergencia
div6 = DivergenceKLtMoP(fExp6,fittMoP6,0,5);
divInter6 = DivergenceKLtMoP(fExp6,interpolar6,0,5);



####################################################################


## AJUSTE 7: Exp(3.08)

set.seed(5)
xExp7 = rexp(1000,3.08);

# ajustes
interpolar7 = Interpolate1param(3.08,tMoPsExp16,1,6,0.05)
fittMoP7 = PolynomialFit(xExp7,maxDegree=7,0,5,1,2);

# divergencia
div7 = DivergenceKLtMoP(fExp7,fittMoP7,0,5);
divInter7 = DivergenceKLtMoP(fExp7,interpolar7,0,5);



####################################################################


## AJUSTE 8: Exp(3.51)

set.seed(5)
xExp8 = rexp(1000,3.51);

# ajustes
interpolar8 = Interpolate1param(3.51,tMoPsExp16,1,6,0.05)
fittMoP8 = PolynomialFit(xExp8,maxDegree=7,0,5,1,2);

# divergencia
div8 = DivergenceKLtMoP(fExp8,fittMoP8,0,5);
divInter8 = DivergenceKLtMoP(fExp8,interpolar8,0,5);



####################################################################


## AJUSTE 9: Exp(3.72)

set.seed(5)
xExp9 = rexp(1000,3.72);

# ajustes
interpolar9 = Interpolate1param(3.72,tMoPsExp16,1,6,0.05)
fittMoP9 = PolynomialFit(xExp9,maxDegree=7,0,5,1,2);

# divergencia
div9 = DivergenceKLtMoP(fExp9,fittMoP9,0,5);
divInter9 = DivergenceKLtMoP(fExp9,interpolar9,0,5);



####################################################################


## AJUSTE 10: Exp(4.11)

set.seed(5)
xExp10 = rexp(1000,4.11);

# ajustes
interpolar10 = Interpolate1param(4.11,tMoPsExp16,1,6,0.05)
fittMoP10 = PolynomialFit(xExp10,maxDegree=7,0,5,1,2);

# divergencia
div10 = DivergenceKLtMoP(fExp10,fittMoP10,0,5);
divInter10 = DivergenceKLtMoP(fExp10,interpolar10,0,5);



####################################################################


## AJUSTE 11: Exp(4.52)

set.seed(5)
xExp11 = rexp(1000,4.52);

# ajustes
interpolar11 = Interpolate1param(4.52,tMoPsExp16,1,6,0.05)
fittMoP11 = PolynomialFit(xExp11,maxDegree=7,0,5,1,2);

# divergencia
div11 = DivergenceKLtMoP(fExp11,fittMoP11,0,5);
divInter11 = DivergenceKLtMoP(fExp11,interpolar11,0,5);



####################################################################


## AJUSTE 12: Exp(4.96)

set.seed(5)
xExp12 = rexp(1000,4.96);

# ajustes
interpolar12 = Interpolate1param(4.96,tMoPsExp16,1,6,0.05)
fittMoP12 = PolynomialFit(xExp12,maxDegree=7,0,5,1,2);

# divergencia
div12 = DivergenceKLtMoP(fExp12,fittMoP12,0,5);
divInter12 = DivergenceKLtMoP(fExp12,interpolar12,0,5);



####################################################################


## AJUSTE 13: Exp(5.28)

set.seed(5)
xExp13 = rexp(1000,5.28);

# ajustes
interpolar13 = Interpolate1param(5.28,tMoPsExp16,1,6,0.05)
fittMoP13 = PolynomialFit(xExp13,maxDegree=7,0,5,1,2);

# divergencia
div13 = DivergenceKLtMoP(fExp13,fittMoP13,0,5);
divInter13 = DivergenceKLtMoP(fExp13,interpolar13,0,5);



####################################################################


## AJUSTE 14: Exp(5.47)

set.seed(5)
xExp14 = rexp(1000,5.47);

# ajustes
interpolar14 = Interpolate1param(5.47,tMoPsExp16,1,6,0.05)
fittMoP14 = PolynomialFit(xExp14,maxDegree=7,0,5,1,2);

# divergencia
div14 = DivergenceKLtMoP(fExp14,fittMoP14,0,5);
divInter14 = DivergenceKLtMoP(fExp14,interpolar14,0,5);



####################################################################


## AJUSTE 15: Exp(5.81)

set.seed(5)
xExp15 = rexp(1000,5.81);

# ajustes
interpolar15 = Interpolate1param(5.81,tMoPsExp16,1,6,0.05)
fittMoP15 = PolynomialFit(xExp15,maxDegree=7,0,5,1,2);

# divergencia
div15 = DivergenceKLtMoP(fExp15,fittMoP15,0,5);
divInter15 = DivergenceKLtMoP(fExp15,interpolar15,0,5);




####################################################################


# divergencias en la columna de resultados
resultados$`DivKL` = c(divInter1,div1,divInter2,div2,divInter3,div3,divInter4,div4,
                       divInter5,div5,divInter6,div6,divInter7,div7,divInter8,div8,
                       divInter9,div9,divInter10,div10,divInter11,div11,divInter12,div12,
                       divInter13,div13,divInter14,div14,divInter15,div15);

resultadosExp = resultados;

save(resultadosExp,file='ResultadosExp.RData');





####################################################################

## GRÁFICAS

####################################################################


## AJUSTE 4: Exp(2.12)

x = seq(0,5,length.out=1000);
ytMoP4 = predicttMoP(fittMoP4,x);
ytMoP4[1] = ytMoP4[2];
ytMoPInter4 = predicttMoP(interpolar4,x);
ytMoPInter4[1] = ytMoPInter4[2];


## GUARDAR EN PDF
pdf("GraficaMetaExp212.pdf",width=7,height=5)

par(mar = c(3,3,3,3))
plot(fExp4,c(0,5),col='#0072B2',type='l',xlim = c(0,5),lwd=4,lty=2)     # original
points(x,ytMoPInter4,col='#D55E00',type='l',lwd=2)
points(x,ytMoP4,col='darkgreen',type='l',lwd=2)
legend(4,2, legend=c("Original", "Interpolar","tMoP"),
       fill=c('#0072B2','#D55E00','darkgreen'),cex=0.8,bty='n')

dev.off()






