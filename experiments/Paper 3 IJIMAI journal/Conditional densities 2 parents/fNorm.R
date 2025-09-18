
# PRUEBA SOBRE APRENDIZAJEDE VARIABLES CONTINUAS CON PADRES CONTINUOS

library(squash);
library(plotly);
library(pracma);
library(ggplot2);
library(plot3D);
library(classInt);
library(MASS);


# ponemos el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# cambiamos el directorio para cargar las funciones
setwd("../../../core")
source("Adjust1Parent.R");          source("Adjust2Parents.R");      source("BIC.R");
source("PolynomialFit.R");          source("ShowModel.R");           source("Prediction.R");
source("AuxiliarFunctions.R");      source("Adjust2ParentsCases.R");

# ponemos el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("DivergenceKL.R");





########################################################################

# DISTRIBUCIONES NORMALES TRIVARIANTES

########################################################################


# X-N(0,1) definida en el intervalo [-3,3]
# Y-N(0,1)  definida en el intervalo [-3,3]
# Z-N(0,1)  definida en el intervalo [-4,4]
# correlación_XY = 0,   correlación_XZ= 1/3,   correlación_YZ = 2/3

meanX = 0;         stdX = 1;
meanY = 0;         stdY = 1;
meanZ = 0;         stdZ = 1;
roXY = 0;     roXZ = 1/3;     roYZ = 2/3;

# generar valores aleatorios usando el paquete MASS
mu = c(0,0,0)
sigma = matrix(c(1,0,1/3,0,1,2/3,1/3,2/3,1),nrow=3,ncol=3)

set.seed(2)
data = mvrnorm(2000,mu,sigma);

# restringir el dominio a [-3,3] en X e Y y a [-4,4] en Z
data = data[data[,1]>=-3& data[,1]<=3& data[,2]>=-3& data[,2]<=3& data[,3]>=-4 &data[,3]<=4,]

# seleccionar 1000 datos y poner la variable hija Z en la primera columna
data = data[1:1000,c(3,1,2)];
colnames(data) = c('Z','X','Y')



# ajuste de variable continua con padre continuo (X e Y son padres de Z)
fittMoP = Adjust2ParentsCCC(data[,1],-4,4,data[,2],-3,3,data[,3],-3,3,n=1,7)

values = fittMoP[[1]]
poly = fittMoP[[2]]

integral(as.polynomial(poly[4,1:7]),values[,10])
plot(as.polynomial(poly[3,1:7]),xlim=c(values[1,7],values[2,7]))



# discretización con equal width (3 intervalos)
a3 = DiscretizationEqWid(data,3);
dataDisc3 = a3[[1]]
breaks3 = rbind(seq(-4,4,length.out=4),seq(-3,3,2),seq(-3,3,2))

# discretización con equal width (5 intervalos)
a5 = DiscretizationEqWid(data,5);
dataDisc5 = a5[[1]]
breaks5 = rbind(seq(-4,4,1.6),seq(-3,3,1.2),seq(-3,3,1.2))


# ajuste discreto (1ª dim = hijo X, 2ª dim = padre X, 3ª dim = padre Y)
fitDisc3 = Adjust2ParentsDDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3,dataDisc3[,3],1:3)
fitDisc5 = Adjust2ParentsDDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5,dataDisc5[,3],1:5)

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1];
area5 = diff(breaks5[1,])[1];
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;


# DIVERGENCIA KL

# X-N(0,1),  Y-N(0,1)     
fX = function(x)        exp(-(x-0)^2/(2*1^2)) / (1*sqrt(2*pi)) / 0.9973002
fY = function(y)        exp(-(y-0)^2/(2*1^2)) / (1*sqrt(2*pi)) / 0.9973002


fZXY = function(z,x,y)
{
  exp(-x^2/8 - y^2/2 - 9*z^2/8 - x*y/2 + 3*x*z/4 + 3*y*z/2) / (2/3 * sqrt(2*pi))
}



pracma::integral(fX,-3,3)
pracma::integral(fY,-3,3)

s = c(seq(-3,-2,0.2),-1.5,-1,0,1,1.5,seq(2,3,0.2))
s = seq(-3,3,0.5)
for (i in 1:length(s))
{
  for (j in 1:length(s))
  {
    f3 = function(z)   fZXY(z,s[i],s[j])
    cat('\n sX =',s[i],'sY =',s[j],'   integral = ',pracma::integral(f3,-4,4))
  }
}


DivKLtmopXYZ(fX,fY,fZXY,fittMoP,-3,3,-3,3,-4,4)
DivKLDiscretXYZ(fX,fY,fZXY,prob3,breaks3)
DivKLDiscretXYZ(fX,fY,fZXY,prob5,breaks5)









########################################################################


# X-N(0,1) definida en el intervalo [-3,3]
# Y-N(0,3)  definida en el intervalo [-6,6]
# Z-N(0,2)  definida en el intervalo [-7,7]
# correlación_XY = 0,   correlación_XZ= 1/2,   correlación_YZ = 1/2

meanX = 0;         stdX = 1;
meanY = 0;         stdY = 3;
meanZ = 0;         stdZ = 2;
roXY = 0;     roXZ = 1/2;     roYZ = 1/2;

# generar valores aleatorios usando el paquete MASS
mu = c(0,0,0)
sigma = matrix(c(1,0,1,0,9,3,1,3,4),nrow=3,ncol=3)


set.seed(2)
data = mvrnorm(2000,mu,sigma);

# restringir el dominio a [-3,3] en X e Y y a [-4,4] en Z
data = data[data[,1]>=-3& data[,1]<=3& data[,2]>=-6& data[,2]<=6& data[,3]>=-7 &data[,3]<=7,]

# seleccionar 1000 datos y poner la variable hija Z en la primera columna
data = data[1:1000,c(3,1,2)];
colnames(data) = c('Z','X','Y')

# ajuste de variable continua con padre continuo (X e Y son padres de Z)
fittMoP = Adjust2ParentsCCC(data[,1],-7,7,data[,2],-3,3,data[,3],-6,6,n=1,7)



# discretización con equal width (3 intervalos)
a3 = DiscretizationEqWid(data,3);
dataDisc3 = a3[[1]]
breaks3 = rbind(seq(-7,7,length.out=4),seq(-3,3,2),seq(-6,6,4))

# discretización con equal width (5 intervalos)
a5 = DiscretizationEqWid(data,5);
dataDisc5 = a5[[1]]
breaks5 = rbind(seq(-7,7,2.8),seq(-3,3,1.2),seq(-6,6,2.4))


# ajuste discreto (1ª dim = hijo X, 2ª dim = padre X, 3ª dim = padre Y)
fitDisc3 = Adjust2ParentsDDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3,dataDisc3[,3],1:3)
fitDisc5 = Adjust2ParentsDDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5,dataDisc5[,3],1:5)

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1];
area5 = diff(breaks5[1,])[1];
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;





# DIVERGENCIA KL

# X-N(0,1),  Y-N(0,3)     
fX = function(x)        exp(-(x-0)^2/(2*1^2)) / (1*sqrt(2*pi)) / 0.9973002
fY = function(y)        exp(-(y-0)^2/(2*3^2)) / (3*sqrt(2*pi)) / 0.9544997


fZXY = function(z,x,y)
{
  exp(-x^2/4 - y^2/36 - z^2/4 - x*y/6 + x*z/2 + y*z/6) / (sqrt(2*pi*2))
}



pracma::integral(fX,-3,3)
pracma::integral(fY,-6,6)

s1 = seq(-3,3,0.5)
s2 = seq(-6,6,1)
for (i in 1:length(s1))
{
  for (j in 1:length(s2))
  {
    f3 = function(z)   fZXY(z,s1[i],s2[j])
    cat('\n sX =',s1[i],'sY =',s2[j],'   integral = ',pracma::integral(f3,-7,7))
  }
}




DivKLtmopXYZ(fX,fY,fZXY,fittMoP,-3,3,-6,6,-7,7)
DivKLDiscretXYZ(fX,fY,fZXY,prob3,breaks3)
DivKLDiscretXYZ(fX,fY,fZXY,prob5,breaks5)






########################################################################


# X-N(0,2) definida en el intervalo [-5,5]
# Y-N(0,1)  definida en el intervalo [-3,3]
# Z-N(0,2)  definida en el intervalo [-6,6]
# correlación_XY = 1/2,   correlación_XZ= 1/4,   correlación_YZ = 1/2

meanX = 0;         stdX = 2;
meanY = 0;         stdY = 1;
meanZ = 0;         stdZ = 2;
roXY = 1/2;     roXZ = 1/4;     roYZ = 1/2;

# generar valores aleatorios usando el paquete MASS
mu = c(0,0,0)
sigma = matrix(c(4,1,1,1,1,1,1,1,4),nrow=3,ncol=3)
inv(sigma)
det(sigma)
-c(1,2,1)%*%inv(sigma)%*%c(1,2,1) / 2


sigma2 = matrix(c(4,1,1,1),nrow=2)


set.seed(2)
data = mvrnorm(2000,mu,sigma);

# restringir el dominio a [-3,3] en X e Y y a [-4,4] en Z
data = data[data[,1]>=-5& data[,1]<=5& data[,2]>=-3& data[,2]<=3& data[,3]>=-6 &data[,3]<=6,]

# seleccionar 1000 datos y poner la variable hija Z en la primera columna
data = data[1:1000,c(3,1,2)];
colnames(data) = c('Z','X','Y')



# ajuste de variable continua con padre continuo (X e Y son padres de Z)
fittMoP = Adjust2ParentsCCC(data[,1],-6,6,data[,2],-5,5,data[,3],-3,3,n=1,7)



# discretización con equal width (3 intervalos)
a3 = DiscretizationEqWid(data,3);
dataDisc3 = a3[[1]]
breaks3 = rbind(seq(-6,6,4),seq(-5,5,length.out=4),seq(-3,3,2))

# discretización con equal width (5 intervalos)
a5 = DiscretizationEqWid(data,5);
dataDisc5 = a5[[1]]
breaks5 = rbind(seq(-6,6,2.4),seq(-5,5,2),seq(-3,3,1.2))


# ajuste discreto (1ª dim = hijo X, 2ª dim = padre X, 3ª dim = padre Y)
fitDisc3 = Adjust2ParentsDDD(dataDisc3[,1],1:3,dataDisc3[,2],1:3,dataDisc3[,3],1:3)
fitDisc5 = Adjust2ParentsDDD(dataDisc5[,1],1:5,dataDisc5[,2],1:5,dataDisc5[,3],1:5)

# cáclulo de las alturas del histograma dividiendo entre el área de los rectángulos
area3 = diff(breaks3[1,])[1];
area5 = diff(breaks5[1,])[1];
prob3 = fitDisc3[[1]] / area3;
prob5 = fitDisc5[[1]] / area5;






# DIVERGENCIA KL

# X-N(0,2),  Y|X-N(x/4,sqrt(3/4))     
fX = function(x)        exp(-(x-0)^2/(2*2^2)) / (2*sqrt(2*pi)) / 0.9875807
fYX = function(y,x)        exp(-(y-x/4)^2 / (2*3/4)) / (sqrt(3/4) * sqrt(2*pi))
fZXY = function(z,x,y)     exp(- y^2/6 - z^2/6 + y*z/3) / (sqrt(2*pi*3))




pracma::integral(fX,-5,5)

s = seq(-5,5,0.5)
for (i in 1:length(s))
{
  f2 = function(y)   fYX(y,s[i])
  cat('\n s =',s[i],'   integral = ',pracma::integral(f2,-3,3))
}


s1 = seq(-6,6,0.5)
s2 = seq(-3,3,1)
for (i in 1:length(s1))
{
  for (j in 1:length(s2))
  {
    f3 = function(z)   fZXY(z,s1[i],s2[j])
    cat('\n sX =',s1[i],'sY =',s2[j],'   integral = ',pracma::integral(f3,-6,6))
  }
}


DivKLtmopXYZ2(fX,fYX,fZXY,fittMoP,-5,5,-3,3,-6,6)
DivKLDiscretXYZ2(fX,fYX,fZXY,prob3,breaks3)
DivKLDiscretXYZ2(fX,fYX,fZXY,prob5,breaks5)











