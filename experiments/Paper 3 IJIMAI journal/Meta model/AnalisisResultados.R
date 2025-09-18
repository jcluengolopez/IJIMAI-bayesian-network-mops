
# analisis de resultados usando el paquete exreport
library(exreport)
library(ggplot2)


# poner el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("Resultados");

# cargar los dataframses con los resultados
load('ResultadosExp.RData');             
load('ResultadosBeta.RData');
load('ResultadostStudent.RData');
load('ResultadosLogNorm.RData')

rBeta = resultadosBeta;
rExp = resultadosExp;
rtStu = resultadostStudent;
rLogNorm = resultadosLogNorm;




##########################################################################

# CONTRASTES DISTRIBUCIÓN BETA

##########################################################################


# TEST 1: TODOS LOS MODELOS

# crear el objeto exreport con % de aciertos
exp1 = expCreate(rBeta,methods='Modelo',problems='Distribucion',
                  name='Interpolation vs tMoPs (Beta)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m1 = aggregate(rBeta$DivKL,list(rBeta$Distribucion),FUN=min);
m11 = rep(m1$x,each=2);
r1Grafica = rBeta;
r1Grafica$DivKL = rBeta$DivKL / m11;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp1Grafica = expCreate(r1Grafica,methods='Modelo',problems='Distribucion',
                         name='Interpolation vs tMoPs (Beta)');

# contraste sobre el error
test1 = testPaired(exp1,output='DivKL',rankOrder='min');

# graficas
plot11 = plotExpSummary(exp1Grafica,'DivKL');
plot12 = plotExpSummary(exp1Grafica,'DivKL',columns=5);


# titulo del report
report1 = exreport('Report for interpolation vs tMoPs (Beta)');

# añadir las diferentes partes del report
report1 = exreportAdd(report1, exp1);
report1 = exreportAdd(report1, test1);
report1 = exreportAdd(report1, list(plot11,plot12));

# añade un resumen al report
table1 = tabularExpSummary(exp1,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report1 = exreportAdd(report1, table1);

# report en html
rr1 = exreportRender(report1,target='html',visualize=T);



##########################################################################


# TEST 2: MODELOS EN [0,1] - [0,2]

rBeta2 = rBeta[1:20,];

# crear el objeto exreport con % de aciertos
exp2 = expCreate(rBeta2,methods='Modelo',problems='Distribucion',
                 name='Interpolation vs tMoPs (Beta in [0,1] - [0,2])');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m2 = aggregate(rBeta2$DivKL,list(rBeta2$Distribucion),FUN=min);
m21 = rep(m2$x,each=2);
r2Grafica = rBeta2;
r2Grafica$DivKL = rBeta2$DivKL / m21;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp2Grafica = expCreate(r2Grafica,methods='Modelo',problems='Distribucion',
                        name='Interpolation vs tMoPs (Beta in [0,1] - [0,2])');

# contraste sobre el error
test2 = testPaired(exp2,output='DivKL',rankOrder='min');

# graficas
plot21 = plotExpSummary(exp2Grafica,'DivKL');
plot22 = plotExpSummary(exp2Grafica,'DivKL',columns=5);


# titulo del report
report2 = exreport('Report for interpolation vs tMoPs (Beta in [0,1] - [0,2])');

# añadir las diferentes partes del report
report2 = exreportAdd(report2, exp2);
report2 = exreportAdd(report2, test2);
report2 = exreportAdd(report2, list(plot21,plot22));

# añade un resumen al report
table2 = tabularExpSummary(exp2,'DivKL',digits=4,format="f",boldfaceColumns='min',
                           tableSplit=2);
report2 = exreportAdd(report2, table2);

# report en html
rr2 = exreportRender(report2,target='html',visualize=T);





##########################################################################

# CONTRASTES DISTRIBUCIÓN EXPONENCIAL

##########################################################################


# TEST 3: TODOS LOS MODELOS

# crear el objeto exreport con % de aciertos
exp3 = expCreate(rExp,methods='Modelo',problems='Distribucion',
                 name='Interpolation vs tMoPs (Exp)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m3 = aggregate(rExp$DivKL,list(rExp$Distribucion),FUN=min);
m31 = rep(m3$x,each=2);
r3Grafica = rExp;
r3Grafica$DivKL = rExp$DivKL / m31;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp3Grafica = expCreate(r3Grafica,methods='Modelo',problems='Distribucion',
                        name='Interpolation vs tMoPs (Exp)');

# contraste sobre el error
test3 = testPaired(exp3,output='DivKL',rankOrder='min');

# graficas
plot31 = plotExpSummary(exp3Grafica,'DivKL');
plot32 = plotExpSummary(exp3Grafica,'DivKL',columns=5);


# titulo del report
report3 = exreport('Report for interpolation vs tMoPs (Exp)');

# añadir las diferentes partes del report
report3 = exreportAdd(report3, exp3);
report3 = exreportAdd(report3, test3);
report3 = exreportAdd(report3, list(plot31,plot32));

# añade un resumen al report
table3 = tabularExpSummary(exp3,'DivKL',digits=4,format="f",boldfaceColumns='min',
                           tableSplit=2);
report3 = exreportAdd(report3, table3);

# report en html
rr3 = exreportRender(report3,target='html',visualize=T);






##########################################################################

# CONTRASTES DISTRIBUCIÓN T-STUDENT

##########################################################################


# TEST 4: TODOS LOS MODELOS

# crear el objeto exreport con % de aciertos
exp4 = expCreate(rtStu,methods='Modelo',problems='Distribucion',
                 name='Interpolation vs tMoPs (t-Student)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m4 = aggregate(rtStu$DivKL,list(rtStu$Distribucion),FUN=min);
m41 = rep(m4$x,each=2);
r4Grafica = rtStu;
r4Grafica$DivKL = rtStu$DivKL / m41;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp4Grafica = expCreate(r4Grafica,methods='Modelo',problems='Distribucion',
                        name='Interpolation vs tMoPs (t-Student)');

# contraste sobre el error
test4 = testPaired(exp4,output='DivKL',rankOrder='min');

# graficas
plot41 = plotExpSummary(exp4Grafica,'DivKL');
plot42 = plotExpSummary(exp4Grafica,'DivKL',columns=5);


# titulo del report
report4 = exreport('Report for interpolation vs tMoPs (t-Student)');

# añadir las diferentes partes del report
report4 = exreportAdd(report4, exp4);
report4 = exreportAdd(report4, test4);
report4 = exreportAdd(report4, list(plot41,plot42));

# añade un resumen al report
table4 = tabularExpSummary(exp4,'DivKL',digits=4,format="f",boldfaceColumns='min',
                           tableSplit=2);
report4 = exreportAdd(report4, table4);

# report en html
rr4 = exreportRender(report4,target='html',visualize=T);







##########################################################################

# CONTRASTES DISTRIBUCIÓN LOGNORMAL

##########################################################################


# TEST 5: TODOS LOS MODELOS

# crear el objeto exreport con % de aciertos
exp5 = expCreate(rLogNorm,methods='Modelo',problems='Distribucion',
                 name='Interpolation vs tMoPs (LogNormal)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m5 = aggregate(rBeta$DivKL,list(rBeta$Distribucion),FUN=min);
m51 = rep(m5$x,each=2);
r5Grafica = rBeta;
r5Grafica$DivKL = rBeta$DivKL / m51;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp5Grafica = expCreate(r5Grafica,methods='Modelo',problems='Distribucion',
                        name='Interpolation vs tMoPs (LogNormal)');

# contraste sobre el error
test5 = testPaired(exp5,output='DivKL',rankOrder='min');

# graficas
plot51 = plotExpSummary(exp5Grafica,'DivKL');
plot52 = plotExpSummary(exp5Grafica,'DivKL',columns=5);


# titulo del report
report5 = exreport('Report for interpolation vs tMoPs (LogNormal)');

# añadir las diferentes partes del report
report5 = exreportAdd(report5, exp5);
report5 = exreportAdd(report5, test5);
report5 = exreportAdd(report5, list(plot51,plot52));

# añade un resumen al report
table5 = tabularExpSummary(exp5,'DivKL',digits=4,format="f",boldfaceColumns='min',
                           tableSplit=2);
report5 = exreportAdd(report5, table5);

# report en html
rr5 = exreportRender(report5,target='html',visualize=T);





