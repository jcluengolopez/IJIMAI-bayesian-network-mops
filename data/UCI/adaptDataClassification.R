

# ADAPT CLASSFIFICATION DATASETS TO USE WITH THE TAN PROGRAM


# ponemos el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))





# DATASET abalone
dataAbalone = read.csv('abalone.csv',header=FALSE,sep=',');
summary(dataAbalone)
dim(dataAbalone)
# la clase es la primera variable y hay que ponerla al final
a = dataAbalone[,1];
dataAbalone[,1] = dataAbalone[,9];
dataAbalone[,9] = a;
# la clase debe ser un carácter
dataAbalone[,9] = as.character(dataAbalone[,9]);
# guardar datos en fichero .RData
save(dataAbalone,file='abalone.RData');




# DATASET appendicitis
dataAppendicitis = read.csv('appendicitis.csv',header=FALSE,sep=',');
summary(dataAppendicitis)
# la clase debe ser un carácter
dataAppendicitis[,8] = as.character(dataAppendicitis[,8]);
# guardar datos en fichero .RData
save(dataAppendicitis,file='appendicitis.RData');




# DATASET australian
dataAustralian = read.csv('australian.csv',header=FALSE,sep=',');
summary(dataAustralian)
dim(dataAustralian)
# la clase debe ser un carácter
dataAustralian[,15] = as.character(dataAustralian[,15]);
# guardar datos en fichero .RData
save(dataAustralian,file='australian.RData');




# DATASET auto
dataAuto = read.csv('auto.csv',header=FALSE,sep=',');
summary(dataAuto)
dim(dataAuto)
# eliminar datos faltantes
B = dataAuto[apply(dataAuto, 1, function(x) !any(x=='?')), , drop=F];
length(B[B=='?']);
dataAuto = B;
# cambiar las variables formadas por cadenas a character, incluida la clase
dataAuto[,3] = as.character(dataAuto[,3]);
dataAuto[,4] = as.character(dataAuto[,4]);
dataAuto[,5] = as.character(dataAuto[,5]);
dataAuto[,6] = as.character(dataAuto[,6]);
dataAuto[,7] = as.character(dataAuto[,7]);
dataAuto[,8] = as.character(dataAuto[,8]);
dataAuto[,9] = as.character(dataAuto[,9]);
dataAuto[,15] = as.character(dataAuto[,15]);
dataAuto[,16] = as.character(dataAuto[,16]);
dataAuto[,18] = as.character(dataAuto[,18]);
# cambiar a numérico las variables que aparecen como factores
dataAuto[,2] = as.double(dataAuto[,2]);
dataAuto[,19] = as.double(dataAuto[,19]);
dataAuto[,20] = as.double(dataAuto[,20]);
dataAuto[,22] = as.double(dataAuto[,22]);
dataAuto[,23] = as.double(dataAuto[,23]);
dataAuto[,26] = as.double(dataAuto[,26]);
# características de los datos
summary(dataAuto)
# guardar datos en fichero .RData
save(dataAuto,file='auto.RData');




# DATASET band
dataBand = read.csv('band.csv',header=FALSE,sep=',');
summary(dataBand)
dim(dataBand)
# eliminar datos faltantes
B = dataBand[apply(dataBand, 1, function(x) !any(x=='?')), , drop=F];
length(B[B=='?']);
dataBand = B;
# cambiar las variables formadas por cadenas a character, incluida la clase
dataBand[,2] = as.character(dataBand[,2]);
dataBand[,3] = as.character(dataBand[,3]);
for (i in 5:15)   {   dataBand[,i] = as.character(dataBand[,i]);   }
dataBand[,18] = as.character(dataBand[,18]);
dataBand[,19] = as.character(dataBand[,19]);
dataBand[,40] = as.character(dataBand[,40]);
# cambiar a numérico las variables que aparecen como factores
dataBand[,1] = as.double(dataBand[,1]);
dataBand[,4] = as.double(dataBand[,4]);
dataBand[,16] = as.double(dataBand[,16]);
for(i in 20:39)   {   dataBand[,i] = as.double(dataBand[,i]);   }
# características de los datos
summary(dataBand)
# guardar datos en fichero .RData
save(dataBand,file='band.RData');




# DATASET cleveland
dataCleveland = read.csv('cleveland.csv',header=FALSE,sep=',');
summary(dataCleveland)
dim(dataCleveland)
# eliminar datos faltantes
B = dataCleveland[apply(dataCleveland, 1, function(x) !any(x=='?')), , drop=F];
length(B[B=='?']);
dataCleveland = B;
# cambiar las variables formadas por cadenas a character, incluida la clase
dataCleveland[,14] = as.character(dataCleveland[,14]);
# cambiar a numérico las variables que aparecen como factores
dataCleveland[,12] = as.double(dataCleveland[,12]);
dataCleveland[,13] = as.double(dataCleveland[,13]);
# características de los datos
summary(dataCleveland)
# guardar datos en fichero .RData
save(dataCleveland,file='cleveland.RData');




# DATASET credit
dataCredit = read.csv('credit.csv',header=FALSE,sep=',');
summary(dataCredit)
dim(dataCredit)
# eliminar datos faltantes
B = dataCredit[apply(dataCredit, 1, function(x) !any(x=='?')), , drop=F];
length(B[B=='?']);
dataCredit = B;
# cambiar las variables formadas por cadenas a character, incluida la clase
dataCredit[,1] = as.character(dataCredit[,1]);
dataCredit[,4] = as.character(dataCredit[,4]);
dataCredit[,5] = as.character(dataCredit[,5]);
dataCredit[,6] = as.character(dataCredit[,6]);
dataCredit[,7] = as.character(dataCredit[,7]);
dataCredit[,9] = as.character(dataCredit[,9]);
dataCredit[,10] = as.character(dataCredit[,10]);
dataCredit[,12] = as.character(dataCredit[,12]);
dataCredit[,13] = as.character(dataCredit[,13]);
dataCredit[,16] = as.character(dataCredit[,16]);
# cambiar a numérico las variables que aparecen como factores
dataCredit[,2] = as.double(dataCredit[,2]);
dataCredit[,14] = as.double(dataCredit[,14]);
# características de los datos
summary(dataCredit)
# guardar datos en fichero .RData
save(dataCredit,file='credit.RData');




# DATASET echocardiogram
dataEchocardiogram = read.csv('echocardiogram.csv',header=FALSE,sep=',');
summary(dataEchocardiogram)
dim(dataEchocardiogram)
# eliminar datos faltantes
B = dataEchocardiogram[apply(dataEchocardiogram, 1, function(x) !any(x=='?')), , drop=F];
length(B[B=='?']);
dataEchocardiogram = B;
# cambiar las variables formadas por cadenas a character, incluida la clase
dataEchocardiogram[,11] = as.character(dataEchocardiogram[,11]);
dataEchocardiogram[,13] = as.character(dataEchocardiogram[,13]);
# cambiar a numérico las variables que aparecen como factores
for (i in 1:10)   {   dataEchocardiogram[,i] = as.double(dataEchocardiogram[,i]);   }
dataEchocardiogram[,12] = as.double(dataEchocardiogram[,12]);
# características de los datos
summary(dataEchocardiogram)
# guardar datos en fichero .RData
save(dataEchocardiogram,file='echocardiogram.RData');




# DATASET fourclass (revisar!!!!)
dataFourclass = read.csv('fourclass.csv',header=FALSE,sep=' ');
summary(dataFourclass)
dim(dataFourclass)
# al ser el espacio el separador, toma una última variable lógica que hay que borrar
dataFourclass = dataFourclass[,1:3];
# la clase es la primera variable y hay que ponerla al final
a = dataFourclass[,1];
dataFourclass[,1] = dataFourclass[,3];
dataFourclass[,3] = a;
# la clase debe ser un carácter
dataFourclass[,3] = as.factor(dataFourclass[,3]);
dataFourclass[,3] = as.character(dataFourclass[,3]);
# guardar datos en fichero .RData
save(dataFourclass,file='fourclass.RData');




# DATASET german
dataGerman = read.csv('german.csv',header=FALSE,sep=',');
summary(dataGerman)
dim(dataGerman)
# se elimina la primera variable porque aparece una en blanco
dataGerman = dataGerman[,2:26];
# la clase debe ser un carácter
dataGerman[,25] = as.character(dataGerman[,25]);
# guardar datos en fichero .RData
save(dataGerman,file='german.RData');




# DATASET glass
dataGlass = read.csv('glass.csv',header=FALSE,sep=',');
summary(dataGlass)
dim(dataGlass)
# la primera columna es id, no es necesaria
dataGlass = dataGlass[,2:11]
# la clase debe ser un carácter
dataGlass[,10] = as.character(dataGlass[,10]);
# guardar datos en fichero .RData
save(dataGlass,file='glass.RData');




# DATASET haberman
dataHaberman = read.csv('haberman.csv',header=FALSE,sep=',');
summary(dataHaberman)
dim(dataHaberman)
# la clase debe ser un carácter
dataHaberman[,4] = as.character(dataHaberman[,4]);
# guardar datos en fichero .RData
save(dataHaberman,file='haberman.RData');




# DATASET hepatitis
dataHepatitis = read.csv('hepatitis.csv',header=FALSE,sep=',');
summary(dataHepatitis)
dim(dataHepatitis)
# eliminar datos faltantes
B = dataHepatitis[apply(dataHepatitis, 1, function(x) !any(x=='?')), , drop=F];
length(B[B=='?']);
dataHepatitis = B;
# cambiar las variables formadas por cadenas a character, incluida la clase
dataHepatitis[,20] = as.character(dataHepatitis[,20]);
# cambiar a numérico las variables que aparecen como factores
dataHepatitis[,4] = as.double(dataHepatitis[,4]);
for(i in 1:6:19)   {   dataHepatitis[,i] = as.double(dataHepatitis[,i]);   }
# características de los datos
summary(dataHepatitis)
# guardar datos en fichero .RData
save(dataHepatitis,file='hepatitis.RData');




# DATASET hungarian 
dataHungarian = read.csv('hungarian.csv',header=FALSE,sep=',');
summary(dataHungarian)
dim(dataHungarian)
# le quito tres variables por ser datos desconocidos
dataHungarian = dataHungarian[,c(1,2,3,4,5,6,7,8,9,10,14)]
# eliminar datos faltantes
B = dataHungarian[apply(dataHungarian, 1, function(x) !any(x=='?')), , drop=F];
length(B[B=='?']);
dataHungarian = B;
# cambiar las variables formadas por cadenas a character, incluida la clase
dataHungarian[,11] = as.character(dataHungarian[,11]);
# cambiar a numérico las variables que aparecen como factores
for (i in 4:9)   {   dataHungarian[,i] = as.double(dataHungarian[,i]);   }
# características de los datos
summary(dataHungarian)
# guardar datos en fichero .RData
save(dataHungarian,file='hungarian.RData');




# DATASET ionosphere
dataIonosphere = read.csv('ionosphere.csv',header=FALSE,sep=',');
summary(dataIonosphere)
dim(dataIonosphere)
# la clase debe ser un carácter
dataIonosphere[,35] = as.character(dataIonosphere[,35]);
# guardar datos en fichero .RData
save(dataIonosphere,file='ionosphere.RData');




# DATASET iris
dataIris = read.csv('iris.csv',header=FALSE,sep=',');
summary(dataIris)
dim(dataIris)
# la clase debe ser un carácter
dataIris[,5] = as.character(dataIris[,5]);
# guardar datos en fichero .RData
save(dataIris,file='iris.RData');




# DATASET liver
dataLiver = read.csv('liver.csv',header=FALSE,sep=',');
summary(dataLiver)
dim(dataLiver)
# la clase debe ser un carácter
dataLiver[,7] = as.character(dataLiver[,7]);
# guardar datos en fichero .RData
save(dataLiver,file='liver.RData');




# DATASET newthyroid
dataNewthyroid = read.csv('newthyroid.csv',header=FALSE,sep=',');
summary(dataNewthyroid)
dim(dataNewthyroid)
# la clase es la primera variable y hay que ponerla al final
a = dataNewthyroid[,1];
dataNewthyroid[,1] = dataNewthyroid[,6];
dataNewthyroid[,6] = a;
# la clase debe ser un carácter
dataNewthyroid[,6] = as.character(dataNewthyroid[,6]);
# guardar datos en fichero .RData
save(dataNewthyroid,file='newthyroid.RData');




# DATASET phoneme
dataPhoneme = read.csv('phoneme.csv',header=FALSE,sep=',');
summary(dataPhoneme)
dim(dataPhoneme)
# la clase debe ser un carácter
dataPhoneme[,6] = as.character(dataPhoneme[,6]);
# guardar datos en fichero .RData
save(dataPhoneme,file='phoneme.RData');




# DATASET switzerland
dataSwitzerland = read.csv('switzerland.csv',header=FALSE,sep=',');
summary(dataSwitzerland)
dim(dataSwitzerland)
# le quito tres variables por ser datos desconocidos
dataSwitzerland = dataSwitzerland[,c(1,2,3,4,5,7,8,9,10,11,14)]
# eliminar datos faltantes
B = dataSwitzerland[apply(dataSwitzerland, 1, function(x) !any(x=='?')), , drop=F];
length(B[B=='?']);
dataSwitzerland = B;
# cambiar las variables formadas por cadenas a character, incluida la clase
dataSwitzerland[,11] = as.character(dataSwitzerland[,11]);
# cambiar a numérico las variables que aparecen como factores
dataSwitzerland[,4] = as.double(dataSwitzerland[,4]);
for (i in 6:10)   {   dataSwitzerland[,i] = as.double(dataSwitzerland[,i]);   }
# características de los datos
summary(dataSwitzerland)
# guardar datos en fichero .RData
save(dataSwitzerland,file='switzerland.RData');




# DATASET vehicle
dataVehicle = read.csv('vehicle.csv',header=FALSE,sep=',');
summary(dataVehicle)
dim(dataVehicle)
# cambiar las variables formadas por cadenas a character, incluida la clase
dataVehicle[,19] = as.character(dataVehicle[,19]);
# guardar datos en fichero .RData
save(dataVehicle,file='vehicle.RData');





# DATASET waveform
dataWaveform = read.csv('waveform.csv',header=FALSE,sep=',');
summary(dataWaveform)
dim(dataWaveform)
# la clase debe ser un carácter
dataWaveform[,22] = as.character(dataWaveform[,22]);
# guardar datos en fichero .RData
save(dataWaveform,file='waveform.RData');




# DATASET wdbc
dataWdbc = read.csv('wdbc.csv',header=FALSE,sep=',');
summary(dataWdbc)
dim(dataWdbc)
# la primera columna es id, no es necesaria
dataWdbc = dataWdbc[,2:32];
# la clase es la primera variable y hay que ponerla al final
a = dataWdbc[,1];
dataWdbc[,1] = dataWdbc[,31];
dataWdbc[,31] = a;
# la clase debe ser un carácter
dataWdbc[,31] = as.character(dataWdbc[,31]);
# guardar datos en fichero .RData
save(dataWdbc,file='wdbc.RData');




# DATASET wine
dataWine = read.csv('wine.csv',header=FALSE,sep=',');
summary(dataWine)
dim(dataWine)
# la clase es la primera variable y hay que ponerla al final
a = dataWine[,1];
dataWine[,1] = dataWine[,14];
dataWine[,14] = a;
# la clase debe ser un carácter
dataWine[,14] = as.character(dataWine[,14]);
# guardar datos en fichero .RData
save(dataWine,file='wine.RData');








