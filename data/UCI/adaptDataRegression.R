

# ADAPT REGRESSION DATASETS TO USE WITH THE TAN PROGRAM


# ponemos el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))






# DATASET adult
dataAdult = read.csv('adult.csv',header=FALSE,sep=',');
summary(dataAdult)
dim(dataAdult)
# eliminar los datos faltantes
B = dataAdult[apply(dataAdult, 1, function(x) !any(x==' ?')), , drop=F]
dataAdult = B;
# cambiar las variables formadas por cadenas a character, incluida la clase
dataAdult[,2] = as.character(dataAdult[,2]);
dataAdult[,4] = as.character(dataAdult[,4]);
dataAdult[,6] = as.character(dataAdult[,6]);
dataAdult[,7] = as.character(dataAdult[,7]);
dataAdult[,8] = as.character(dataAdult[,8]);
dataAdult[,9] = as.character(dataAdult[,9]);
dataAdult[,10] = as.character(dataAdult[,10]);
dataAdult[,14] = as.character(dataAdult[,14]);
# la clase debe ser un car?cter
dataAdult[,15] = as.character(dataAdult[,15]);
# guardar datos en fichero .RData
save(dataAdult,file='adult.RData');


# DATASET airfoil
dataAirfoil = read.csv('airfoil.csv',header=FALSE,sep=';');
summary(dataAirfoil)
dim(dataAirfoil)
# guardar datos en fichero .RData
save(dataAirfoil,file='airfoil.csv');




# DATASET auto
dataAuto = read.csv('auto.csv',header=FALSE,sep=';');
summary(dataAuto)
dim(dataAuto)
# la clase es la primera variable y hay que ponerla al final
a = dataAuto[,1];
dataAuto[,1] = dataAuto[,8];
dataAuto[,8] = a;
# guardar datos en fichero .RData
save(dataAuto,file='auto.RData');



# DATASET bikeSharingDay (se quitan las dos primeras porque son ?ndices y las dos
# pen?ltimas porque al sumarlas se obtiene la clase)
dataBikeSharingDay = read.csv('bikeSharingDay.csv',header=TRUE,sep=',');
summary(dataBikeSharingDay)
dim(dataBikeSharingDay)
dataBikeSharingDay = dataBikeSharingDay[,c(3:13,16)];
# guardar datos en fichero .RData
save(dataBikeSharingDay,file='bikeSharingDay.RData');




# DATASET bikeSharingHour (se quitan las dos primeras porque son ?ndices y las dos
# pen?ltimas porque al sumarlas se obtiene la clase)
dataBikeSharingHour = read.csv('bikeSharingHour.csv',header=TRUE,sep=',');
summary(dataBikeSharingHour)
dim(dataBikeSharingHour)
dataBikeSharingHour = dataBikeSharingHour[,c(3:14,17)];
# guardar datos en fichero .RData
save(dataBikeSharingHour,file='bikeSharingHour.RData');




# DATASET bodyfat
dataBodyfat = read.csv('bodyfat.csv',header=FALSE);
summary(dataBodyfat)
dim(dataBodyfat)
# guardar datos en fichero .RData
save(dataBodyfat,file='bodyfat.RData');



# DATASET chickenpox
dataChickenpox = read.csv('chickenpox.csv',header=TRUE,sep=',');
summary(dataChickenpox)
dim(dataChickenpox)
# eliminar la primera variable
dataChickenpox = dataChickenpox[,-1];
# guardar datos en fichero .RData
save(dataChickenpox,file='chickenpox.RData');





# DATASET cloud
dataCloud = read.csv('cloud.csv',header=FALSE);
summary(dataCloud)
dim(dataCloud)
# cambiar las variables formadas por cadenas a character
dataCloud[,1] = as.double(dataCloud[,1]);
dataCloud[,2] = as.character(dataCloud[,2]);
dataCloud[,3] = as.character(dataCloud[,3]);
# guardar datos en fichero .RData
save(dataCloud,file='cloud.RData');




# DATASET concrete
dataConcrete = read.table('concrete.txt',header=FALSE);
summary(dataConcrete)
dim(dataConcrete)
# pasar a num?ricas las variables que son factores
for (i in 1:9)   {   dataConcrete[,i] = as.numeric(dataConcrete[,i]);   }
# guardar datos en fichero .RData
save(dataConcrete,file='concrete.RData');




# DATASET election
dataElection = read.csv('election.csv',header=TRUE);
summary(dataElection)
dim(dataElection)
dataElection = dataElection[,-2];
# guardar datos en fichero .RData
save(dataElection,file='election.RData');






# DATASET forecasting
dataForecasting = read.csv('forecasting.csv',header=TRUE,sep=';');
summary(dataForecasting)
dim(dataForecasting)
# guardar datos en fichero .RData
save(dataForecasting,file='forecasting.RData');




# DATASET forestfires
dataForestfires = read.csv('forestfires.csv',header=FALSE,sep=',');
summary(dataForestfires)
dim(dataForestfires)
# cambiar las variables formadas por cadenas a character
dataForestfires[,3] = as.character(dataForestfires[,3]);
dataForestfires[,4] = as.character(dataForestfires[,4]);
# guardar datos en fichero .RData
save(dataForestfires,file='forestfires.RData');




# DATASET housing
dataHousing = read.csv('housing.csv',header=FALSE,sep=',');
summary(dataHousing)
dim(dataHousing)
# guardar datos en fichero .RData
save(dataHousing,file='housing.RData');



# DATASET istanbul
dataIstanbul = read.csv('istanbul.csv',header=TRUE,sep=';');
summary(dataIstanbul)
dim(dataIstanbul)
# la clase es la primera variable y hay que ponerla al final
a = dataIstanbul[,1];     b = colnames(dataIstanbul)[1]
dataIstanbul[,1] = dataIstanbul[,9];
colnames(dataIstanbul)[1] = colnames(dataIstanbul)[9];
dataIstanbul[,9] = a;
colnames(dataIstanbul)[9] = b;
# guardar datos en fichero .RData
save(dataIstanbul,file='istanbul.RData');




# DATASET localization
dataLocalization = read.csv('localization.csv',header=TRUE,sep=',');
summary(dataLocalization)
dim(dataLocalization)
# la ?litma columna es la desviaci?n estandar de ALE, que es la variable a predecir
dataLocalization = dataLocalization[,-6];
# guardar datos en fichero .RData
save(dataLocalization,file='localization.RData');




# DATASET machine
dataMachine = read.csv('machine.csv',header=FALSE,sep=',');
summary(dataMachine)
dim(dataMachine)
# cambiar las variables formadas por factores a caracteres
dataMachine[,1] = as.character(dataMachine[,1]);
# quitar la variable 2 porque son nombres ?nicos
dataMachine = dataMachine[,-2];
# le pongo bien los nombres de las variables
names(dataMachine) = c('V1','V2','V3','V4','V5','V6','V7','V8','V9');
# guardar datos en fichero .RData
save(dataMachine,file='machine.RData');




# DATASET metro
dataMetro = read.csv('metro.csv',header=TRUE,sep=',');
summary(dataMetro)
dim(dataMetro)
# eliminar la columna de fechas y horas
dataMetro = dataMetro[,-8];
# guardar datos en fichero .RData
save(dataMetro,file='metro.RData');




# DATASET pollution
dataPollution = read.table('pollution.txt',header=FALSE);
summary(dataPollution)
dim(dataPollution)
# guardar datos en fichero .RData
save(dataPollution,file='pollution.RData');




# DATASET powerplant
dataPowerplant = read.csv('powerplant.csv',header=TRUE,sep=';');
summary(dataPowerplant)
dim(dataPowerplant)
# guardar datos en fichero .RData
save(dataPowerplant,file='powerplant.RData');




# DATASET servo
dataServo = read.csv('servo.csv',header=FALSE,sep=',');
summary(dataServo)
dim(dataServo)
# cambiar las variables formadas por factores a caracteres
dataServo[,1] = as.character(dataServo[,1]);
dataServo[,2] = as.character(dataServo[,2]);
# guardar datos en fichero .RData
save(dataServo,file='servo.RData');




# DATASET strikes
dataStrikes = read.table('strikes.txt',header=FALSE);
summary(dataStrikes)
dim(dataStrikes)
# la clase al final. Es la variable 3 (volumeen de jornadas de huelga)
aux = dataStrikes[,3];
dataStrikes[,3] = dataStrikes[,7];
dataStrikes[,7] = aux;
# guardar datos en fichero .RData
save(dataStrikes,file='strikes.RData');




# DATASET synchronous
dataSynchronous = read.csv('synchronous.csv',header=TRUE,sep=';');
summary(dataSynchronous)
dim(dataSynchronous)
# guardar datos en fichero .RData
save(dataSynchronous,file='synchronous.RData');





# DATASET toxicity
dataToxicity = read.table('toxicity.csv',header=FALSE,sep=';');
summary(dataToxicity)
dim(dataToxicity)
# guardar datos en fichero .RData
save(dataToxicity,file='toxicity.RData');




# DATASET trafficSaoPaulo
dataTrafficSaoPaulo = read.table('trafficSaoPaulo.csv',header=TRUE,sep=';');
summary(dataTrafficSaoPaulo)
dim(dataTrafficSaoPaulo)
# guardar datos en fichero .RData
save(dataTrafficSaoPaulo,file='trafficSaoPaulo.RData');



# DATASET veteran
dataVeteran = read.table('veteran.txt',header=FALSE);
summary(dataVeteran)
dim(dataVeteran)
# la clase al final. Es la variable 3 (d?as de supervivencia)
aux = dataVeteran[,3];
dataVeteran[,3] = dataVeteran[,8];
dataVeteran[,8] = aux;
# guardar datos en fichero .RData
save(dataVeteran,file='veteran.RData');









