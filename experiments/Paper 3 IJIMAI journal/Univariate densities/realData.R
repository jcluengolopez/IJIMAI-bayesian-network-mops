
# PRUEBA SOBRE APRENDIZAJEDE VARIABLES CONTINUAS CON PADRES CONTINUOS

library(MoTBFs)
library(polynom)


# ponemos el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# cambiamos el directorio para cargar las funciones
setwd("../../../core")
source("Adjust1Parent.R");          source("Adjust2Parents.R");      source("BIC.R");
source("PolynomialFit.R");          source("ShowModel.R");           source("Prediction.R");
source("AuxiliarFunctions.R");      source("Adjust2ParentsCases.R"); source("Simulation.R");


# ponemos el directorio de este script y lo cambiamos para cargar los datos
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../../../data/UCI")
load('appendicitis.RData');     load('cloud.RData');          load('concrete.RData');
load('housing.RData');          load('hungarian.RData');

#load('istanbul.RData');         load('adult.RData');         load('forestfires.RData');
#load('machine.RData');          load('strikes.RData');

# Volvemos al directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))






####################################################################

# tmop y Pérez-Bernabé
# - Grados: 7, 4 y 3


# data frame para los resultados
var = c(rep('Appendicitis-V2',6),rep('Appendicitis-V3',6),rep('Appendicitis-V4',6),
        rep('Cloud-V4',6),rep('Cloud-V6',6),rep('Cloud-V8',6),
        rep('Concrete-V1',6),rep('Concrete-V3',6),rep('Concrete-V6',6),
        rep('Housing-V5',6),rep('Housing-V7',6),rep('Housing-V14',6),
        rep('Hungarian-V1',6),rep('Hungarian-V5',6),rep('Hungarian-V8',6));

#distribucion = c(distribucion1,distribucion2,distribucion3,distribucion4,distribucion5)
metodo = rep(c(rep('Perez-Bernabe',3),rep('tMoPs',3)),15);
grado = rep(c(7,4,3),30);

resultados = data.frame(var,metodo,grado,rep(0,90))
colnames(resultados) = c('Variable','Método','Grado','BIC');





########################################################################

# DATAAPPENDICITIS (V2, V3, V4)

########################################################################


varTypeAp = DiscOrCont(dataAppendicitis,1);
which(varTypeAp==1)
hist(dataAppendicitis[,4],breaks=30)

dataMinAp = as.numeric(apply(dataAppendicitis,2,min));
dataMaxAp = as.numeric(apply(dataAppendicitis,2,max));

# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrainAp = sample(1:nrow(dataAppendicitis),floor(nrow(dataAppendicitis)*0.8));
dataTrainAp = dataAppendicitis[selectTrainAp,];
dataTestAp = dataAppendicitis[-selectTrainAp,];



########################################################################

# variable V2

min1 = dataMinAp[2];                   max1 = dataMaxAp[2];
dataTrain1 = dataTrainAp[,2];          dataTest1 = dataTestAp[,2];

grado = c(7,4,3);
bicPB1 = rep(0,3);
bictmop1 = rep(0,3);

for (i in 1:3)
{
  # ajuste Pérez-Bernabé
  fitPB1 = univMoTBF(dataTrain1,"MOP",c(min1,max1),maxParam=grado[i]+1);
  fitPBpoly1 = coef(fitPB1);
  dataValPB1 = cbind(c(min1,max1),c(min1,max1),c(0,0))

  # ajuste tmop
  fittmop1 = PolynomialFit(dataTrain1,maxDegree=grado[i],min1,max1,1,2);
  fittmopPoly1 = fittmop1[[2]];
  dataValtmop1 = cbind(fittmop1[[1]],c(min1,max1),fittmop1[[3]])

  # BIC
  bicPB1[i] = BIC0Parent(fitPBpoly1,dataTest1,dataValPB1,1)
  bictmop1[i] = BIC0Parent(fittmopPoly1,dataTest1,dataValtmop1,1)
}




########################################################################

# variable V3

min2 = dataMinAp[3];                   max2 = dataMaxAp[3];
dataTrain2 = dataTrainAp[,3];          dataTest2 = dataTestAp[,3];

grado = c(7,4,3);
bicPB2 = rep(0,3);
bictmop2 = rep(0,3);

for (i in 1:3)
{
  # ajuste Pérez-Bernabé
  fitPB2 = univMoTBF(dataTrain2,"MOP",c(min2,max2),maxParam=grado[i]+1);
  fitPBpoly2 = coef(fitPB2);
  dataValPB2 = cbind(c(min2,max2),c(min2,max2),c(0,0))

  # ajuste tmop
  fittmop2 = PolynomialFit(dataTrain2,maxDegree=grado[i],min2,max2,1,2);
  fittmopPoly2 = fittmop2[[2]];
  dataValtmop2 = cbind(fittmop2[[1]],c(min2,max2),fittmop2[[3]])

  # BIC
  bicPB2[i] = BIC0Parent(fitPBpoly2,dataTest2,dataValPB2,1)
  bictmop2[i] = BIC0Parent(fittmopPoly2,dataTest2,dataValtmop2,1)
}




########################################################################

# variable V4

min3 = dataMinAp[4];                   max3 = dataMaxAp[4];
dataTrain3 = dataTrainAp[,4];          dataTest3 = dataTestAp[,4];

grado = c(7,4,3);
bicPB3 = rep(0,3);
bictmop3 = rep(0,3);

for (i in 1:3)
{
  # ajuste Pérez-Bernabé
  fitPB3 = univMoTBF(dataTrain3,"MOP",c(min3,max3),maxParam=grado[i]+1);
  fitPBpoly3 = coef(fitPB3);
  dataValPB3 = cbind(c(min3,max3),c(min3,max3),c(0,0))

  # ajuste tmop
  fittmop3 = PolynomialFit(dataTrain3,maxDegree=grado[i],min3,max3,1,2);
  fittmopPoly3 = fittmop3[[2]];
  dataValtmop3 = cbind(fittmop3[[1]],c(min3,max3),fittmop3[[3]])

  # BIC
  bicPB3[i] = BIC0Parent(fitPBpoly3,dataTest3,dataValPB3,1)
  bictmop3[i] = BIC0Parent(fittmopPoly3,dataTest3,dataValtmop3,1)
}





########################################################################

# DATACLOUD (V4, V6, V8)

########################################################################


varTypeCl = DiscOrCont(dataCloud,1);
which(varTypeCl==1)
hist(dataCloud[,4],breaks=30)

dataMinCl = as.numeric(apply(dataCloud,2,min));
dataMaxCl = as.numeric(apply(dataCloud,2,max));

# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrainCl = sample(1:nrow(dataCloud),floor(nrow(dataCloud)*0.8));
dataTrainCl = dataCloud[selectTrainCl,];
dataTestCl = dataCloud[-selectTrainCl,];



########################################################################

# variable V4

min4 = dataMinCl[4];                   max4 = dataMaxCl[4];
dataTrain4 = dataTrainCl[,4];          dataTest4 = dataTestCl[,4];

grado = c(7,4,3);
bicPB4 = rep(0,3);
bictmop4 = rep(0,3);

for (i in 1:3)
{
  # ajuste Pérez-Bernabé
  fitPB4 = univMoTBF(dataTrain4,"MOP",c(min4,max4),maxParam=grado[i]+1);
  fitPBpoly4 = coef(fitPB4);
  dataValPB4 = cbind(c(min4,max4),c(min4,max4),c(0,0))

  # ajuste tmop
  fittmop4 = PolynomialFit(dataTrain4,maxDegree=grado[i],min4,max4,1,2);
  fittmopPoly4 = fittmop4[[2]];
  dataValtmop4 = cbind(fittmop4[[1]],c(min4,max4),fittmop4[[3]])

  # BIC
  bicPB4[i] = BIC0Parent(fitPBpoly4,dataTest4,dataValPB4,1)
  bictmop4[i] = BIC0Parent(fittmopPoly4,dataTest4,dataValtmop4,1)
}



########################################################################

# variable V6

min5 = dataMinCl[6];                   max5 = dataMaxCl[6];
dataTrain5 = dataTrainCl[,6];          dataTest5 = dataTestCl[,6];

grado = c(7,4,3);
bicPB5 = rep(0,3);
bictmop5 = rep(0,3);

for (i in 1:3)
{
  # ajuste Pérez-Bernabé
  fitPB5 = univMoTBF(dataTrain5,"MOP",c(min5,max5),maxParam=grado[i]+1);
  fitPBpoly5 = coef(fitPB5);
  dataValPB5 = cbind(c(min5,max5),c(min5,max5),c(0,0))

  # ajuste tmop
  fittmop5 = PolynomialFit(dataTrain5,maxDegree=grado[i],min5,max5,1,2);
  fittmopPoly5 = fittmop5[[2]];
  dataValtmop5 = cbind(fittmop5[[1]],c(min5,max5),fittmop5[[3]])

  # BIC
  bicPB5[i] = BIC0Parent(fitPBpoly5,dataTest5,dataValPB5,1)
  bictmop5[i] = BIC0Parent(fittmopPoly5,dataTest5,dataValtmop5,1)
}





########################################################################

# variable V8

min6 = dataMinCl[8];                   max6 = dataMaxCl[8];
dataTrain6 = dataTrainCl[,8];          dataTest6 = dataTestCl[,8];

grado = c(7,4,3);
bicPB6 = rep(0,3);
bictmop6 = rep(0,3);

for (i in 1:3)
{
  # ajuste Pérez-Bernabé
  fitPB6 = univMoTBF(dataTrain6,"MOP",c(min6,max6),maxParam=grado[i]+1);
  fitPBpoly6 = coef(fitPB6);
  dataValPB6 = cbind(c(min6,max6),c(min6,max6),c(0,0))

  # ajuste tmop
  fittmop6 = PolynomialFit(dataTrain6,maxDegree=grado[i],min6,max6,1,2);
  fittmopPoly6 = fittmop6[[2]];
  dataValtmop6 = cbind(fittmop6[[1]],c(min6,max6),fittmop6[[3]])

  # BIC
  bicPB6[i] = BIC0Parent(fitPBpoly6,dataTest6,dataValPB6,1)
  bictmop6[i] = BIC0Parent(fittmopPoly6,dataTest6,dataValtmop6,1)
}









########################################################################

# DATACONCRETE (V1, V3, V6)

########################################################################


varTypeCo = DiscOrCont(dataConcrete,1);
which(varTypeCo==1)
hist(dataConcrete[,2],breaks=30)

dataMinCo = as.numeric(apply(dataConcrete,2,min));
dataMaxCo = as.numeric(apply(dataConcrete,2,max));

# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrainCo = sample(1:nrow(dataConcrete),floor(nrow(dataConcrete)*0.8));
dataTrainCo = dataConcrete[selectTrainCo,];
dataTestCo = dataConcrete[-selectTrainCo,];




########################################################################

# variable V1

min7 = dataMinCo[1];                   max7 = dataMaxCo[1];
dataTrain7 = dataTrainCo[,1];          dataTest7 = dataTestCo[,1];

grado = c(7,4,3);
bicPB7 = rep(0,3);
bictmop7 = rep(0,3);

for (i in 1:3)
{
  # ajuste Pérez-Bernabé
  fitPB7 = univMoTBF(dataTrain7,"MOP",c(min7,max7),maxParam=grado[i]+1);
  fitPBpoly7 = coef(fitPB7);
  dataValPB7 = cbind(c(min7,max7),c(min7,max7),c(0,0))

  # ajuste tmop
  fittmop7 = PolynomialFit(dataTrain7,maxDegree=grado[i],min7,max7,1,2);
  fittmopPoly7 = fittmop7[[2]];
  dataValtmop7 = cbind(fittmop7[[1]],c(min7,max7),fittmop7[[3]])

  # BIC
  bicPB7[i] = BIC0Parent(fitPBpoly7,dataTest7,dataValPB7,1)
  bictmop7[i] = BIC0Parent(fittmopPoly7,dataTest7,dataValtmop7,1)
}




########################################################################

# variable V3

min8 = dataMinCo[3];                   max8 = dataMaxCo[3];
dataTrain8 = dataTrainCo[,3];          dataTest8 = dataTestCo[,3];

grado = c(7,4,3);
bicPB8 = rep(0,3);
bictmop8 = rep(0,3);

for (i in 1:3)
{
  # ajuste Pérez-Bernabé
  fitPB8 = univMoTBF(dataTrain8,"MOP",c(min8,max8),maxParam=grado[i]+1);
  fitPBpoly8 = coef(fitPB8);
  dataValPB8 = cbind(c(min8,max8),c(min8,max8),c(0,0))

  # ajuste tmop
  fittmop8 = PolynomialFit(dataTrain8,maxDegree=grado[i],min8,max8,1,2);
  fittmopPoly8 = fittmop8[[2]];
  dataValtmop8 = cbind(fittmop8[[1]],c(min8,max8),fittmop8[[3]])

  # BIC
  bicPB8[i] = BIC0Parent(fitPBpoly8,dataTest8,dataValPB8,1)
  bictmop8[i] = BIC0Parent(fittmopPoly8,dataTest8,dataValtmop8,1)
}






########################################################################

# variable V6

min9 = dataMinCo[6];                   max9 = dataMaxCo[6];
dataTrain9 = dataTrainCo[,6];          dataTest9 = dataTestCo[,6];

grado = c(7,4,3);
bicPB9 = rep(0,3);
bictmop9 = rep(0,3);

for (i in 1:3)
{
  # ajuste Pérez-Bernabé
  fitPB9 = univMoTBF(dataTrain9,"MOP",c(min9,max9),maxParam=grado[i]+1);
  fitPBpoly9 = coef(fitPB9);
  dataValPB9 = cbind(c(min9,max9),c(min9,max9),c(0,0))

  # ajuste tmop
  fittmop9 = PolynomialFit(dataTrain9,maxDegree=grado[i],min9,max9,1,2);
  fittmopPoly9 = fittmop9[[2]];
  dataValtmop9 = cbind(fittmop9[[1]],c(min9,max9),fittmop9[[3]])

  # BIC
  bicPB9[i] = BIC0Parent(fitPBpoly9,dataTest9,dataValPB9,1)
  bictmop9[i] = BIC0Parent(fittmopPoly9,dataTest9,dataValtmop9,1)
}









########################################################################

# DATAHOUSING (V5, V7, V14)

########################################################################


varTypeHo = DiscOrCont(dataHousing,1);
which(varTypeHo==1)
hist(dataHousing[,14],breaks=30)

dataMinHo = as.numeric(apply(dataHousing,2,min));
dataMaxHo = as.numeric(apply(dataHousing,2,max));

# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrainHo = sample(1:nrow(dataHousing),floor(nrow(dataHousing)*0.8));
dataTrainHo = dataHousing[selectTrainHo,];
dataTestHo = dataHousing[-selectTrainHo,];




########################################################################

# variable V5

min10 = dataMinHo[5];                   max10 = dataMaxHo[5];
dataTrain10 = dataTrainHo[,5];          dataTest10 = dataTestHo[,5];

grado = c(7,4,3);
bicPB10 = rep(0,3);
bictmop10 = rep(0,3);

for (i in 1:3)
{
  # ajuste Pérez-Bernabé
  fitPB10 = univMoTBF(dataTrain10,"MOP",c(min10,max10),maxParam=grado[i]+1);
  fitPBpoly10 = coef(fitPB10);
  dataValPB10 = cbind(c(min10,max10),c(min10,max10),c(0,0))

  # ajuste tmop
  fittmop10 = PolynomialFit(dataTrain10,maxDegree=grado[i],min10,max10,1,2);
  fittmopPoly10 = fittmop10[[2]];
  dataValtmop10 = cbind(fittmop10[[1]],c(min10,max10),fittmop10[[3]])

  # BIC
  bicPB10[i] = BIC0Parent(fitPBpoly10,dataTest10,dataValPB10,1)
  bictmop10[i] = BIC0Parent(fittmopPoly10,dataTest10,dataValtmop10,1)
}




########################################################################

# variable V7

min11 = dataMinHo[7];                   max11 = dataMaxHo[7];
dataTrain11 = dataTrainHo[,7];          dataTest11 = dataTestHo[,7];

grado = c(7,4,3);
bicPB11 = rep(0,3);
bictmop11 = rep(0,3);

for (i in 1:3)
{
  # ajuste Pérez-Bernabé
  fitPB11 = univMoTBF(dataTrain11,"MOP",c(min11,max11),maxParam=grado[i]+1);
  fitPBpoly11 = coef(fitPB11);
  dataValPB11 = cbind(c(min11,max11),c(min11,max11),c(0,0))

  # ajuste tmop
  fittmop11 = PolynomialFit(dataTrain11,maxDegree=grado[i],min11,max11,1,2);
  fittmopPoly11 = fittmop11[[2]];
  dataValtmop11 = cbind(fittmop11[[1]],c(min11,max11),fittmop11[[3]])

  # BIC
  bicPB11[i] = BIC0Parent(fitPBpoly11,dataTest11,dataValPB11,1)
  bictmop11[i] = BIC0Parent(fittmopPoly11,dataTest11,dataValtmop11,1)
}




########################################################################

# variable V14

min12 = dataMinHo[14];                   max12 = dataMaxHo[14];
dataTrain12 = dataTrainHo[,14];          dataTest12 = dataTestHo[,14];

grado = c(7,4,3);
bicPB12 = rep(0,3);
bictmop12 = rep(0,3);

for (i in 1:3)
{
  # ajuste Pérez-Bernabé
  fitPB12 = univMoTBF(dataTrain12,"MOP",c(min12,max12),maxParam=grado[i]+1);
  fitPBpoly12 = coef(fitPB12);
  dataValPB12 = cbind(c(min12,max12),c(min12,max12),c(0,0))

  # ajuste tmop
  fittmop12 = PolynomialFit(dataTrain12,maxDegree=grado[i],min12,max12,1,2);
  fittmopPoly12 = fittmop12[[2]];
  dataValtmop12 = cbind(fittmop12[[1]],c(min12,max12),fittmop12[[3]])

  # BIC
  bicPB12[i] = BIC0Parent(fitPBpoly12,dataTest12,dataValPB12,1)
  bictmop12[i] = BIC0Parent(fittmopPoly12,dataTest12,dataValtmop12,1)
}






########################################################################

# DATAHUNGARIAN

########################################################################


# Variables V1, V5 y V8 de dataHungarian
varTypeHu = DiscOrCont(dataHungarian,1);
which(varTypeHu==1)
hist(dataHungarian[,8],breaks=20)

dataMinHu = as.numeric(apply(dataHungarian,2,min));
dataMaxHu = as.numeric(apply(dataHungarian,2,max));

# dataTrain y dataTest (80% - 20%)
set.seed(2)
selectTrainHu = sample(1:nrow(dataHungarian),floor(nrow(dataHungarian)*0.8));
dataTrainHu = dataHungarian[selectTrainHu,];
dataTestHu = dataHungarian[-selectTrainHu,];




########################################################################

# variable V1

min13 = dataMinHu[1];                   max13 = dataMaxHu[1];
dataTrain13 = dataTrainHu[,1];          dataTest13 = dataTestHu[,1];

grado = c(7,4,3);
bicPB13 = rep(0,3);
bictmop13 = rep(0,3);

for (i in 1:3)
{
  # ajuste Pérez-Bernabé
  fitPB13 = univMoTBF(dataTrain13,"MOP",c(min13,max13),maxParam=grado[i]+1);
  fitPBpoly13 = coef(fitPB13);
  dataValPB13 = cbind(c(min13,max13),c(min13,max13),c(0,0))

  # ajuste tmop
  fittmop13 = PolynomialFit(dataTrain13,maxDegree=grado[i],min13,max13,1,2);
  fittmopPoly13 = fittmop13[[2]];
  dataValtmop13 = cbind(fittmop13[[1]],c(min13,max13),fittmop13[[3]])

  # BIC
  bicPB13[i] = BIC0Parent(fitPBpoly13,dataTest13,dataValPB13,1)
  bictmop13[i] = BIC0Parent(fittmopPoly13,dataTest13,dataValtmop13,1)
}



########################################################################

# variable V5

min14 = dataMinHu[5];                   max14 = dataMaxHu[5];
dataTrain14 = dataTrainHu[,5];          dataTest14 = dataTestHu[,5];

grado = c(7,4,3);
bicPB14 = rep(0,3);
bictmop14 = rep(0,3);

for (i in 1:3)
{
  # ajuste Pérez-Bernabé
  fitPB14 = univMoTBF(dataTrain14,"MOP",c(min14,max14),maxParam=grado[i]+1);
  fitPBpoly14 = coef(fitPB14);
  dataValPB14 = cbind(c(min14,max14),c(min14,max14),c(0,0))

  # ajuste tmop
  fittmop14 = PolynomialFit(dataTrain14,maxDegree=grado[i],min14,max14,1,2);
  fittmopPoly14 = fittmop14[[2]];
  dataValtmop14 = cbind(fittmop14[[1]],c(min14,max14),fittmop14[[3]])

  # BIC
  bicPB14[i] = BIC0Parent(fitPBpoly14,dataTest14,dataValPB14,1)
  bictmop14[i] = BIC0Parent(fittmopPoly14,dataTest14,dataValtmop14,1)
}




########################################################################

# variable V8

min15 = dataMinHu[8];                   max15 = dataMaxHu[8];
dataTrain15 = dataTrainHu[,8];          dataTest15 = dataTestHu[,8];

grado = c(7,4,3);
bicPB15 = rep(0,3);
bictmop15 = rep(0,3);

for (i in 1:3)
{
  # ajuste Pérez-Bernabé
  fitPB15 = univMoTBF(dataTrain15,"MOP",c(min15,max15),maxParam=grado[i]+1);
  fitPBpoly15 = coef(fitPB15);
  dataValPB15 = cbind(c(min15,max15),c(min15,max15),c(0,0))

  # ajuste tmop
  fittmop15 = PolynomialFit(dataTrain15,maxDegree=grado[i],min15,max15,1,2);
  fittmopPoly15 = fittmop15[[2]];
  dataValtmop15 = cbind(fittmop15[[1]],c(min15,max15),fittmop15[[3]])

  # BIC
  bicPB15[i] = BIC0Parent(fitPBpoly15,dataTest15,dataValPB15,1)
  bictmop15[i] = BIC0Parent(fittmopPoly15,dataTest15,dataValtmop15,1)
}





resultados$BIC = c(bicPB1,bictmop1,bicPB2,bictmop2,bicPB3,bictmop3,bicPB4,bictmop4,
                   bicPB5,bictmop5,bicPB6,bictmop6,bicPB7,bictmop7,bicPB8,bictmop8,
                   bicPB9,bictmop9,bicPB10,bictmop10,bicPB11,bictmop11,bicPB12,bictmop12,
                   bicPB13,bictmop13,bicPB14,bictmop14,bicPB15,bictmop15)


resultadosRealData = resultados;

save(resultadosRealData,file='ResultadosRealData.RData');






########################################################################

# GRÁFICAS

########################################################################

# V1 de Concrete

x = seq(min4,max4,length.out=1000)

# Pérez Bernabé
fitPB4 = univMoTBF(dataTrain4,"MOP",c(min4,max4),maxParam=5);
yPB = predict(as.polynomial(coef(fitPB4)),x);

# tmop
fittmop4 = PolynomialFit(dataTrain4,maxDegree=4,min4,max4,1,2);
ytmop = predicttMoP(fittmop4,x);
ytmop[1] = ytmop[2];
ytmop[1000] = ytmop[999];

## GUARDAR EN PDF
pdf("GraficaUnivHousingV12.pdf",width=8,height=5)

par(mar = c(4,4,2,2))
plot(x,ytmop,col='#0072B2',type='l',xlim = c(min4,max4),lwd=2);
legend(x=40,y=0.245,legend=c('MoTBF','tMoP'),col=c('#CC79A7','#0072B2'),lwd=2.5,
       box.lty=0,cex=0.9,lty=c(0,1,1),pch=20,x.intersp=2,y.intersp=1.35,seg.len=0.7)
points(x,yPB,type='l',lwd=2,col='#CC79A7')
hist(dataTrain4,probability=TRUE,add=TRUE,breaks=30,col=rgb(255,215,181,100,
                                                               maxColorValue=255))

dev.off()







########################################################################

# V14 de Housing

x = seq(min12,max12,length.out=1000)

# Pérez Bernabé
fitPB12 = univMoTBF(dataTrain12,"MOP",c(min12,max12),maxParam=5);
yPB = predict(as.polynomial(coef(fitPB12)),x);

# tmop
fittmop12 = PolynomialFit(dataTrain12,maxDegree=4,min12,max12,1,2);
ytmop = predicttMoP(fittmop12,x);
ytmop[1000] = ytmop[999];


## GUARDAR EN PDF
pdf("GraficaUnivHousingV14.pdf",width=7,height=5)

par(mar = c(4,4,2,2))
par(mar = c(2, 2, 0.5, 0.5))
plot(x,ytmop,col='#0072B2',type='l',xlim = c(min12,max12),ylim=c(0,0.07),lwd=2,
     ylab=NA,xlab=NA,cex.axis=1.1);
points(x,yPB,type='l',lwd=2,col='#a23d75')
hist(dataTrain12,probability=TRUE,add=TRUE,breaks=30,col=rgb(255,215,181,100,
                                                               maxColorValue=255))
legend(40,0.06, legend=c("tMoP","MoTBF"),seg.len=2,lwd=2,
       col=c('#0072B2','#a23d75'),cex=1.1,bty='n',lty=c(1,1))

dev.off()

















