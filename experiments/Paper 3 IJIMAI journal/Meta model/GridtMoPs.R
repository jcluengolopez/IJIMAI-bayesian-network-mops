
# CARGAR DEPENDENCIAS AL PRINCIPIO (polinomios)
library(polynom)


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

# poner el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))




####################################################################

# DISTRIBUCIÓN BETA CON PARÁMETROS ALPHA Y BETA ENTRE 0 Y 1 CADA 0,1


# la distancia entre los puntos de la malla queda como un parámetro
d = 0.1;
puntos = seq(0,1,d);

# crear lista para almacenar todos los tMoPs
tMoPs = vector('list',length(puntos)^2);

# para cada uno de los puntos de la malla
for (i in 1:length(puntos))
{
  for (j in 1:length(puntos))
  {
    # muestra aleatoria (fijando la semilla)
    set.seed(i+j)
    xBeta = rbeta(1000,puntos[i],puntos[j]);
    
    # crear el tMoP
    n = (i-1) * length(puntos) + j     # fila i, columna j de la malla
    tMoPs[[n]] = PolynomialFit(xBeta,maxDegree=7,0,1,1,2);
  }
}

tMoPsBeta0101 = tMoPs;

save(tMoPsBeta0101,file='tMoPsBeta0101.RData');

####################################################################





####################################################################

# DISTRIBUCIÓN BETA CON PARÁMETROS ALPHA ENTRE 0 Y 1 Y BETA ENTRE 1 Y 2 CADA 0,1


# la distancia entre los puntos de la malla queda como un parámetro
d = 0.1;
puntos1 = seq(0,1,d);
puntos2 = seq(1,2,d);

# crear lista para almacenar todos los tMoPs
tMoPs = vector('list',length(puntos1)*length(puntos2));

# para cada uno de los puntos de la malla
for (i in 1:length(puntos1))
{
  for (j in 1:length(puntos2))
  {
    # muestra aleatoria (fijando la semilla)
    set.seed(i+j)
    xBeta = rbeta(1000,puntos1[i],puntos2[j]);
    
    # crear el tMoP
    n = (i-1) * length(puntos1) + j     # fila i, columna j de la malla
    tMoPs[[n]] = PolynomialFit(xBeta,maxDegree=7,0,1,1,2);
  }
}

tMoPsBeta0112 = tMoPs;

save(tMoPsBeta0112,file='tMoPsBeta0112.RData');

####################################################################




####################################################################

# DISTRIBUCIÓN BETA CON PARÁMETROS ALPHA Y BETA ENTRE 1 Y 2 CADA 0,1


# la distancia entre los puntos de la malla queda como un parámetro
d = 0.1;
puntos = seq(1,2,d);

# crear lista para almacenar todos los tMoPs
tMoPs = vector('list',length(puntos)^2);

# para cada uno de los puntos de la malla
for (i in 1:length(puntos))
{
  for (j in 1:length(puntos))
  {
    # muestra aleatoria (fijando la semilla)
    set.seed(i+j)
    xBeta = rbeta(1000,puntos[i],puntos[j]);
    
    # crear el tMoP
    n = (i-1) * length(puntos) + j     # fila i, columna j de la malla
    tMoPs[[n]] = PolynomialFit(xBeta,maxDegree=7,0,1,1,2);
  }
}

tMoPsBeta1212 = tMoPs;

save(tMoPsBeta1212,file='tMoPsBeta1212.RData');

####################################################################




####################################################################

# DISTRIBUCIÓN EXPONENCIAL CON PARÁMETRO k ENTRE 1 Y 6 CADA 0.05 (ajuste en [0,5])


# la distancia entre los puntos de la malla queda como un parámetro
d = 0.05;
puntos = seq(1,6,d);

# crear lista para almacenar todos los tMoPs
tMoPs = vector('list',length(puntos));

# para cada uno de los puntos del vector
for (i in 1:length(puntos))
{
  # muestra aleatoria (fijando la semilla)
  set.seed(i)
  xExp = rexp(1000,puntos[i]);
  
  # crear el tMoP
  tMoPs[[i]] = PolynomialFit(xExp,maxDegree=7,0,5,1,2);
}

tMoPsExp16 = tMoPs;

save(tMoPsExp16,file='tMoPsExp16.RData');

####################################################################




####################################################################

# DISTRIBUCIÓN T-STUDENT CON PARÁMETRO k ENTRE 1 Y 6 CADA 0.05 (intervalo -20,20)


# la distancia entre los puntos de la malla queda como un parámetro
d = 0.05;
puntos = seq(1,6,d);

# crear lista para almacenar todos los tMoPs
tMoPs = vector('list',length(puntos));

# para cada uno de los puntos del vector
for (i in 1:length(puntos))
{
  # muestra aleatoria (fijando la semilla)
  set.seed(i)
  xtStu = rt(1000,puntos[i]);
  
  # crear el tMoP
  tMoPs[[i]] = PolynomialFit(xtStu,maxDegree=7,-20,20,1,2);
}

tMoPstStu16 = tMoPs;

save(tMoPstStu16,file='tMoPstStu16.RData');

####################################################################






####################################################################

# DISTRIBUCIÓN LOGNORMAL CON PARÁMETROS mu ENTRE 0 Y 1 Y sigma ENTRE 0.2 Y 1.2 CADA 0,1
# intervalo [0,15]


# la distancia entre los puntos de la malla queda como un parámetro
d = 0.1;
puntos1 = seq(0,1,d);
puntos2 = seq(0.2,1.2,d);

# gráficas para probar parámetros
#f = function(x)   1/(x*sqrt(2*pi*1)) * exp(-(log(x)-2)^2/(2*1))
#plot(f,xlim=c(0,5))
#pracma::integral(f,0,15)

# crear lista para almacenar todos los tMoPs
tMoPs = vector('list',length(puntos1)*length(puntos2));

# para cada uno de los puntos de la malla
for (i in 1:length(puntos1))
{
  for (j in 1:length(puntos2))
  {
    # muestra aleatoria (fijando la semilla)
    set.seed(i+j)
    xLogNorm = rlnorm(1000,puntos1[i],puntos2[j]);
    
    # crear el tMoP
    n = (i-1) * length(puntos1) + j     # fila i, columna j de la malla
    tMoPs[[n]] = PolynomialFit(xLogNorm,maxDegree=7,0,15,1,2);
  }
}

tMoPsLogNorm01212 = tMoPs;

save(tMoPsLogNorm01212,file='tMoPsLogNorm01212.RData');

####################################################################





