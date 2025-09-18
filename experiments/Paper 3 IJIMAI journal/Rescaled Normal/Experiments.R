
# CARGAR DEPENDENCIAS AL PRINCIPIO (polinomios)
library(polynom)


# poner el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# cargar función de divergencia y tipificación
source("DivergenceKL.R");         source("RescalefNorm.R");


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
fNorm01 = function(x)  exp(-x^2/2)/sqrt(2*pi)


set.seed(2)
xNorm01 = rnorm(2000,0,1);           xNorm01 = xNorm01[xNorm01>=-3 & xNorm01<=3];
xNorm01 = xNorm01[1:1000];


# tmop
fittMoP01 = PolynomialFit(xNorm01,maxDegree=7,-3,3,1,2);



####################################################################


# DISTRIBUCIÓN NORMAL (2,4) en el intervalo [-10,14]
fNorm24 = function(x)  exp(-(x-2)^2/(2*4^2))/(4*sqrt(2*pi))

set.seed(2)
xNorm24 = rnorm(2000,2,4);           xNorm01 = xNorm01[xNorm01>=-10 & xNorm01<=14];
xNorm24 = xNorm24[1:1000];

# ajuste tmop y tipificado
fittMoP24tip = RescalefNorm(fittMoP01,2,4,-3,3,-10,14)
fittMoP24 = PolynomialFit(xNorm24,maxDegree=7,-10,14,1,2);

# divergencia
DivergenceKLtMoP(fNorm24,fittMoP24tip,-10,14);
DivergenceKLtMoP(fNorm24,fittMoP24,-10,14);





####################################################################

## GRÁFICAS

####################################################################


## AJUSTE 1: N(2,4)

x = seq(-10,14,length.out=10000);
ytMoP24 = predicttMoP(fittMoP24,x);
ytMoP24[1] = ytMoP24[2];
#ytMoP24[1000] = ytMoP24[999];
ytMoPInter24 = predicttMoP(fittMoP24tip,x);
ytMoPInter24[1] = ytMoPInter24[2];
#ytMoPInter24[1000] = ytMoPInter24[999];


## GUARDAR EN PDF
pdf("GraficaMetaNorm2-4.pdf",width=7,height=5)

par(mar = c(3,3,3,3))
plot(fNorm24,c(-10,14),col='#0072B2',type='l',xlim = c(-10,14),lwd=4,lty=2)     # original
points(x,ytMoPInter24,col='#D55E00',type='l',lwd=3)
points(x,ytMoP24,col='darkgreen',type='l',lwd=2)
legend(8.5,0.09, legend=c("Theoretical", "Standarization", "Direct tMoP fit"),seg.len=2,lwd=2,
       col=c('#0072B2','#D55E00','darkgreen'),cex=0.8,bty='n',lty=c(2,1,1))

dev.off()











p1 = c(3,2,-1)
p2 = c(3,2/3,-1/9)
p3 = c(3,2/2,-1/4)

plot(as.polynomial(p1),xlim=c(-2,4))
plot(as.polynomial(p2),xlim=c(-2,4))
solve(as.polynomial(p1))
solve(as.polynomial(p2))
solve(as.polynomial(p3))


q1 = c(6,-5,-2,1)
q2 = c(6,-5/3,-2/9,1/27)

plot(as.polynomial(q1),xlim=c(-3,4))
plot(as.polynomial(q2),xlim=c(-3,4))
solve(as.polynomial(q1))
solve(as.polynomial(q2))
solve(as.polynomial(q3))


coef(change.origin(as.polynomial(p1),2))





