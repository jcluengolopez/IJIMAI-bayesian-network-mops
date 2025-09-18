
# analisis de resultados usando el paquete exreport
library(exreport)
library(ggplot2)


# poner el directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("Resultados");

# cargar los dataframses con los resultados
load('ResultadosExp.RData');               load('ResultadosNorm.RData');
load('ResultadosGam.RData');               load('ResultadosBeta.RData');
load('ResultadosLogNorm.RData');           load('ResultadosBiMod.RData');
load('ResultadosChi2.RData');              load('ResultadostStudent.RData');
load('ResultadosWeibull.RData');           load('ResultadosRealData.RData');

r = rbind(resultadosExp,resultadosNorm,resultadosGam,resultadosBeta,resultadosLogNorm,
          resultadosBiMod,resultadosChi2,resultadostStudent,resultadosWeibull);
rr = resultadosRealData;

# volver al directorio de este script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# acortar los nombres de las distribuciones
r$f = gsub('gamma','Ga',r$f);             r$f = gsub('beta','Be',r$f)
r$f = gsub('lognorm','LN',r$f);           r$f = gsub('bimod','biM',r$f);
r$f = gsub('t-student','tStu',r$f);       r$f = gsub('weibull','Wei',r$f);
r$f = gsub('chi-2','Chi2',r$f);



# ordenar la matriz de combinaciones por distribucion
r = r[order(r$f),];

# quitar el guion al nombre de la divergencia
colnames(r)[7] = 'DivKL'

# crear columna con número de intervalos y grado polinomio conjuntamente
r$nInter[r$Método=='Shenoy' & r$nInter==1 & r$Grado==7] = '1-7';
r$nInter[r$Método=='Shenoy' & r$nInter==1 & r$Grado==4] = '1-4';
r$nInter[r$Método=='Shenoy' & r$nInter==1 & r$Grado==3] = '1-3';
r$nInter[r$Método=='Shenoy' & r$nInter==2 & r$Grado==7] = '2-7';
r$nInter[r$Método=='Shenoy' & r$nInter==2 & r$Grado==4] = '2-4';
r$nInter[r$Método=='Shenoy' & r$nInter==2 & r$Grado==3] = '2-3';
r$nInter[r$Método=='Shenoy' & r$nInter==3 & r$Grado==7] = '3-7';
r$nInter[r$Método=='Shenoy' & r$nInter==3 & r$Grado==4] = '3-4';
r$nInter[r$Método=='Shenoy' & r$nInter==3 & r$Grado==3] = '3-3';
r$nInter[r$Método=='Shenoy' & r$nInter==4 & r$Grado==7] = '4-7';
r$nInter[r$Método=='Shenoy' & r$nInter==4 & r$Grado==4] = '4-4';
r$nInter[r$Método=='Shenoy' & r$nInter==4 & r$Grado==3] = '4-3';
r$nInter[r$Método=='Shenoy' & r$nInter==4 & r$Grado==2] = '4-2';






##########################################################################

# CONTRASTES MÚLTIPLES PARA TAMAÑO MUESTRA EN PÉREZ-BERNABÉ Y tMoPs

##########################################################################


# TEST 10: PÉREZ-BERNABÉ (GRADO=7)

r10 = r[r$Método=='Perez-Bernabe' & r$Grado==7,]

# crear el objeto exreport con % de aciertos
exp10 = expCreate(r10,methods='Muestra',problems='f',
                  name='Sample size in Perez-Bernabe (deg=7)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m10 = aggregate(r10$DivKL,list(r10$f),FUN=min);
m101 = rep(m10$x,each=3);
r10Grafica = r10;
r10Grafica$DivKL = r10$DivKL / m101;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp10Grafica = expCreate(r10Grafica,methods='Muestra',problems='f',
                         name='Sample size in Perez-Bernabe (deg=7)');


# contraste sobre el error
test10 = testMultipleControl(exp10,output='DivKL',rankOrder='min');
test10b = testMultiplePairwise(exp10,output='DivKL',rankOrder='min');
test10Grafica = testMultipleControl(exp10Grafica,output='DivKL',rankOrder='min');

# graficas
plot101 = plotExpSummary(exp10Grafica,'DivKL');
plot102 = plotExpSummary(exp10Grafica, 'DivKL', columns=5);
plot103 = plotCumulativeRank(test10Grafica);
plot104 = plotRankDistribution(test10Grafica);

# primera table
table101 = tabularTestSummary(test10,columns=c('pvalue','rank','wtl'));
table101b = tabularTestPairwise(test10b);


# titulo del report
report10 = exreport('Report for sample size in Perez-Bernabe deg=7 (50, 100 or 1000)');

# añadir las diferentes partes del report
report10 = exreportAdd(report10, exp10);
report10 = exreportAdd(report10, test10);
report10 = exreportAdd(report10, list(plot101,plot102,plot103,plot104,table101,table101b));

# añade un resumen al report
table102 = tabularExpSummary(exp10,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report10 = exreportAdd(report10, table102);

# report en html
r10 = exreportRender(report10,target='html',visualize=T);




##########################################################################


# TEST 11: tMoPs (GRADO=7)

r11 = r[r$Método=='tMoPs' & r$Grado==7,]

# crear el objeto exreport con % de aciertos
exp11 = expCreate(r11,methods='Muestra',problems='f',name='Sample size in tMoPs (deg=7)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m11 = aggregate(r11$DivKL,list(r11$f),FUN=min);
m111 = rep(m11$x,each=3);
r11Grafica = r11;
r11Grafica$DivKL = r11$DivKL / m111;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp11Grafica = expCreate(r11Grafica,methods='Muestra',problems='f',
                         name='Sample size in tMoPs (deg=7)');


# contraste sobre el error
test11 = testMultipleControl(exp11,output='DivKL',rankOrder='min');
test11b = testMultiplePairwise(exp11,output='DivKL',rankOrder='min');
test11Grafica = testMultipleControl(exp11Grafica,output='DivKL',rankOrder='min');

# graficas
plot111 = plotExpSummary(exp11Grafica,'DivKL');
plot112 = plotExpSummary(exp11Grafica, 'DivKL', columns=5);
plot113 = plotCumulativeRank(test11Grafica);
plot114 = plotRankDistribution(test11Grafica);

# primera table
table111 = tabularTestSummary(test11,columns=c('pvalue','rank','wtl'));
table111b = tabularTestPairwise(test11b);


# titulo del report
report11 = exreport('Report for sample size in tMoPs deg=7 (50, 100 or 1000)');

# añadir las diferentes partes del report
report11 = exreportAdd(report11, exp11);
report11 = exreportAdd(report11, test11);
report11 = exreportAdd(report11, list(plot111,plot112,plot113,plot114,table111,table111b));

# añade un resumen al report
table112 = tabularExpSummary(exp11,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report11 = exreportAdd(report11, table112);

# report en html
r11 = exreportRender(report11,target='html',visualize=T);




##########################################################################


# TEST 20: PÉREZ-BERNABÉ (GRADO=3)

r20 = r[r$Método=='Perez-Bernabe' & r$Grado==3,]

# crear el objeto exreport con % de aciertos
exp20 = expCreate(r20,methods='Muestra',problems='f',
                  name='Sample size in Perez-Bernabe (deg=3)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m20 = aggregate(r20$DivKL,list(r20$f),FUN=min);
m201 = rep(m20$x,each=3);
r20Grafica = r20;
r20Grafica$DivKL = r20$DivKL / m201;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp20Grafica = expCreate(r20Grafica,methods='Muestra',problems='f',
                         name='Sample size in Perez-Bernabe (deg=3)');


# contraste sobre el error
test20 = testMultipleControl(exp20,output='DivKL',rankOrder='min');
test20b = testMultiplePairwise(exp20,output='DivKL',rankOrder='min');
test20Grafica = testMultipleControl(exp20Grafica,output='DivKL',rankOrder='min');

# graficas
plot201 = plotExpSummary(exp20Grafica,'DivKL');
plot202 = plotExpSummary(exp20Grafica, 'DivKL', columns=5);
plot203 = plotCumulativeRank(test20Grafica);
plot204 = plotRankDistribution(test20Grafica);

# primera table
table201 = tabularTestSummary(test20,columns=c('pvalue','rank','wtl'));
table201b = tabularTestPairwise(test20b);


# titulo del report
report20 = exreport('Report for sample size in Perez-Bernabe deg=3 (50, 100 or 1000)');

# añadir las diferentes partes del report
report20 = exreportAdd(report20, exp20);
report20 = exreportAdd(report20, test20);
report20 = exreportAdd(report20, list(plot201,plot202,plot203,plot204,table201,table201b));

# añade un resumen al report
table202 = tabularExpSummary(exp20,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report20 = exreportAdd(report20, table202);

# report en html
r20 = exreportRender(report20,target='html',visualize=T);



##########################################################################


# TEST 21: tMoPs (GRADO=3)

r21 = r[r$Método=='tMoPs' & r$Grado==3,]

# crear el objeto exreport con % de aciertos
exp21 = expCreate(r21,methods='Muestra',problems='f',name='Sample size in tMoPs (deg=3)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m21 = aggregate(r21$DivKL,list(r21$f),FUN=min);
m211 = rep(m21$x,each=3);
r21Grafica = r21;
r21Grafica$DivKL = r21$DivKL / m211;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp21Grafica = expCreate(r21Grafica,methods='Muestra',problems='f',
                         name='Sample size in tMoPs (deg=3)');


# contraste sobre el error
test21 = testMultipleControl(exp21,output='DivKL',rankOrder='min');
test21b = testMultiplePairwise(exp21,output='DivKL',rankOrder='min');
test21Grafica = testMultipleControl(exp21Grafica,output='DivKL',rankOrder='min');

# graficas
plot211 = plotExpSummary(exp21Grafica,'DivKL');
plot212 = plotExpSummary(exp21Grafica, 'DivKL', columns=5);
plot213 = plotCumulativeRank(test21Grafica);
plot214 = plotRankDistribution(test21Grafica);

# primera table
table211 = tabularTestSummary(test21,columns=c('pvalue','rank','wtl'));
table211b = tabularTestPairwise(test21b);


# titulo del report
report21 = exreport('Report for sample size in tMoPs deg=3 (50, 100 or 1000)');

# añadir las diferentes partes del report
report21 = exreportAdd(report21, exp21);
report21 = exreportAdd(report21, test21);
report21 = exreportAdd(report21, list(plot211,plot212,plot213,plot214,table211,table211b));

# añade un resumen al report
table212 = tabularExpSummary(exp21,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report21 = exreportAdd(report21, table212);

# report en html
r21 = exreportRender(report21,target='html',visualize=T);





##########################################################################


# TEST 22: PÉREZ-BERNABÉ (GRADO=4)

r22 = r[r$Método=='Perez-Bernabe' & r$Grado==4,]

# crear el objeto exreport con % de aciertos
exp22 = expCreate(r22,methods='Muestra',problems='f',
                  name='Sample size in Perez-Bernabe (deg=4)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m22 = aggregate(r22$DivKL,list(r22$f),FUN=min);
m221 = rep(m22$x,each=3);
r22Grafica = r22;
r22Grafica$DivKL = r22$DivKL / m221;

# crear el objeto exreport con % de aciertos relativos para grafica
exp22Grafica = expCreate(r22Grafica,methods='Muestra',problems='f',
                         name='Sample size in Perez-Bernabe (deg=4)');


# contraste sobre el error
test22 = testMultipleControl(exp22,output='DivKL',rankOrder='min');
test22b = testMultiplePairwise(exp22,output='DivKL',rankOrder='min');
test22Grafica = testMultipleControl(exp22Grafica,output='DivKL',rankOrder='min');

# graficas
plot221 = plotExpSummary(exp22Grafica,'DivKL');
plot222 = plotExpSummary(exp22Grafica, 'DivKL', columns=5);
plot223 = plotCumulativeRank(test22Grafica);
plot224 = plotRankDistribution(test22Grafica);

# primera table
table221 = tabularTestSummary(test22,columns=c('pvalue','rank','wtl'));
table221b = tabularTestPairwise(test22b);


# titulo del report
report22 = exreport('Report for sample size in Perez-Bernabe deg=4 (50, 100 or 1000)');

# añadir las diferentes partes del report
report22 = exreportAdd(report22, exp22);
report22 = exreportAdd(report22, test22);
report22 = exreportAdd(report22, list(plot221,plot222,plot223,plot224,table221,table221b));

# añade un resumen al report
table222 = tabularExpSummary(exp22,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report22 = exreportAdd(report22, table222);

# report en html
r22 = exreportRender(report22,target='html',visualize=T);




##########################################################################


# TEST 23: tMoPs (GRADO=4)

r23 = r[r$Método=='tMoPs' & r$Grado==4,]

# crear el objeto exreport con % de aciertos
exp23 = expCreate(r23,methods='Muestra',problems='f',name='Sample size in tMoPs (deg=4)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m23 = aggregate(r23$DivKL,list(r23$f),FUN=min);
m231 = rep(m23$x,each=3);
r23Grafica = r23;
r23Grafica$DivKL = r23$DivKL / m231;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp23Grafica = expCreate(r23Grafica,methods='Muestra',problems='f',
                         name='Sample size in tMoPs (deg=4)');


# contraste sobre el error
test23 = testMultipleControl(exp23,output='DivKL',rankOrder='min');
test23b = testMultiplePairwise(exp23,output='DivKL',rankOrder='min');
test23Grafica = testMultipleControl(exp23Grafica,output='DivKL',rankOrder='min');

# graficas
plot231 = plotExpSummary(exp23Grafica,'DivKL');
plot232 = plotExpSummary(exp23Grafica, 'DivKL', columns=5);
plot233 = plotCumulativeRank(test23Grafica);
plot234 = plotRankDistribution(test23Grafica);

# primera table
table231 = tabularTestSummary(test23,columns=c('pvalue','rank','wtl'));
table231b = tabularTestPairwise(test23b);


# titulo del report
report23 = exreport('Report for sample size in tMoPs deg=4 (50, 100 or 1000)');

# añadir las diferentes partes del report
report23 = exreportAdd(report23, exp23);
report23 = exreportAdd(report23, test23);
report23 = exreportAdd(report23, list(plot231,plot232,plot233,plot234,table231,table231b));

# añade un resumen al report
table232 = tabularExpSummary(exp23,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report23 = exreportAdd(report23, table232);

# report en html
r23 = exreportRender(report23,target='html',visualize=T);





##########################################################################

# CONTRASTES MÚLTIPLES PARA GRADO EN tMoPs

##########################################################################


# TEST 24: tMoP (TAMAÑO DE MUESTRA=1000)

r24 = r[r$Método=='tMoPs' & r$Muestra==1000,]

# crear el objeto exreport con % de aciertos
exp24 = expCreate(r24,methods='Grado',problems='f',name='Degree tMoPs (sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m24 = aggregate(r24$DivKL,list(r24$f),FUN=min);
m241 = rep(m24$x,each=3);
r24Grafica = r24;
r24Grafica$DivKL = r24$DivKL / m241;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp24Grafica = expCreate(r24Grafica,methods='Grado',problems='f',
                         name='Degree tMoPs (sample=1000)');


# contraste sobre el error
test24 = testMultipleControl(exp24,output='DivKL',rankOrder='min');
test24b = testMultiplePairwise(exp24,output='DivKL',rankOrder='min');
test24Grafica = testMultipleControl(exp24Grafica,output='DivKL',rankOrder='min');

# graficas
plot241 = plotExpSummary(exp24Grafica,'DivKL');
plot242 = plotExpSummary(exp24Grafica, 'DivKL', columns=5);
plot243 = plotCumulativeRank(test24Grafica);
plot244 = plotRankDistribution(test24Grafica);

# primera table
table241 = tabularTestSummary(test24,columns=c('pvalue','rank','wtl'));
table241b = tabularTestPairwise(test24b);


# titulo del report
report24 = exreport('Report degree tMoPs sample=1000 (3, 4 or 7)');

# añadir las diferentes partes del report
report24 = exreportAdd(report24, exp24);
report24 = exreportAdd(report24, test24);
report24 = exreportAdd(report24, list(plot241,plot242,plot243,plot244,table241,table241b));

# añade un resumen al report
table242 = tabularExpSummary(exp24,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report24 = exreportAdd(report24, table242);

# report en html
r24 = exreportRender(report24,target='html',visualize=T);



##########################################################################


# TEST 25: tMoP (TAMAÑO DE MUESTRA=100)

r25 = r[r$Método=='tMoPs' & r$Muestra==100,]

# crear el objeto exreport con % de aciertos
exp25 = expCreate(r25,methods='Grado',problems='f',name='Degree tMoPs (sample=100)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m25 = aggregate(r25$DivKL,list(r25$f),FUN=min);
m251 = rep(m25$x,each=3);
r25Grafica = r25;
r25Grafica$DivKL = r25$DivKL / m251;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp25Grafica = expCreate(r25Grafica,methods='Grado',problems='f',
                         name='Degree tMoPs (sample=100)');


# contraste sobre el error
test25 = testMultipleControl(exp25,output='DivKL',rankOrder='min');
test25b = testMultiplePairwise(exp25,output='DivKL',rankOrder='min');
test25Grafica = testMultipleControl(exp25Grafica,output='DivKL',rankOrder='min');

# graficas
plot251 = plotExpSummary(exp25Grafica,'DivKL');
plot252 = plotExpSummary(exp25Grafica, 'DivKL', columns=5);
plot253 = plotCumulativeRank(test25Grafica);
plot254 = plotRankDistribution(test25Grafica);

# primera table
table251 = tabularTestSummary(test25,columns=c('pvalue','rank','wtl'));
table251b = tabularTestPairwise(test25b);


# titulo del report
report25 = exreport('Report degree tMoPs sample=100 (3, 4 or 7)');

# añadir las diferentes partes del report
report25 = exreportAdd(report25, exp25);
report25 = exreportAdd(report25, test25);
report25 = exreportAdd(report25, list(plot251,plot252,plot253,plot254,table251,table251b));

# añade un resumen al report
table252 = tabularExpSummary(exp25,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report25 = exreportAdd(report25, table252);

# report en html
r25 = exreportRender(report25,target='html',visualize=T);



##########################################################################


# TEST 26: tMoP (TAMAÑO DE MUESTRA=50)

r26 = r[r$Método=='tMoPs' & r$Muestra==50,]

# crear el objeto exreport con % de aciertos
exp26 = expCreate(r26,methods='Grado',problems='f',name='Degree tMoPs (sample=50)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m26 = aggregate(r26$DivKL,list(r26$f),FUN=min);
m261 = rep(m26$x,each=3);
r26Grafica = r26;
r26Grafica$DivKL = r26$DivKL / m261;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp26Grafica = expCreate(r26Grafica,methods='Grado',problems='f',
                         name='Degree tMoPs (sample=50)');


# contraste sobre el error
test26 = testMultipleControl(exp26,output='DivKL',rankOrder='min');
test26b = testMultiplePairwise(exp26,output='DivKL',rankOrder='min');
test26Grafica = testMultipleControl(exp26Grafica,output='DivKL',rankOrder='min');

# graficas
plot261 = plotExpSummary(exp26Grafica,'DivKL');
plot262 = plotExpSummary(exp26Grafica, 'DivKL', columns=5);
plot263 = plotCumulativeRank(test26Grafica);
plot264 = plotRankDistribution(test26Grafica);

# primera table
table261 = tabularTestSummary(test26,columns=c('pvalue','rank','wtl'));
table261b = tabularTestPairwise(test26b);


# titulo del report
report26 = exreport('Report degree tMoPs sample=50 (3, 4 or 7)');

# añadir las diferentes partes del report
report26 = exreportAdd(report26, exp26);
report26 = exreportAdd(report26, test26);
report26 = exreportAdd(report26, list(plot261,plot262,plot263,plot264,table261,table261b));

# añade un resumen al report
table262 = tabularExpSummary(exp26,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report26 = exreportAdd(report26, table262);

# report en html
r26 = exreportRender(report26,target='html',visualize=T);






##########################################################################

# CONTRASTES MÚLTIPLES PARA NÚMERO DE INTERVALOS Y GRADO POLINOMIO EN MÉTODO SHENOY

##########################################################################


# TEST 12: NÚMERO DE INTERVALOS Y GRADO POLINOMIO PARA MÉTODO SHENOY

r12 = r[r$Método=='Shenoy',]

# crear el objeto exreport con % de aciertos
exp12 = expCreate(r12,methods='nInter',problems='f',name='nIntervals and degree in Shenoy');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m12 = aggregate(r12$DivKL,list(r12$f),FUN=min);
m121 = rep(m12$x,each=7);
r12Grafica = r12;
r12Grafica$DivKL = r12$DivKL / m121;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp12Grafica = expCreate(r12Grafica,methods='nInter',problems='f',
                         name='nIntervals and degree in Shenoy');


# contraste sobre el error
test12 = testMultipleControl(exp12,output='DivKL',rankOrder='min');
test12b = testMultiplePairwise(exp12,output='DivKL',rankOrder='min');
test12Grafica = testMultipleControl(exp12Grafica,output='DivKL',rankOrder='min');

# graficas
plot121 = plotExpSummary(exp12Grafica,'DivKL');
plot122 = plotExpSummary(exp12Grafica, 'DivKL', columns=5);
plot123 = plotCumulativeRank(test12Grafica);
plot124 = plotRankDistribution(test12Grafica);

# primera table
table121 = tabularTestSummary(test12,columns=c('pvalue','rank','wtl'));
table121b = tabularTestPairwise(test12b);


# titulo del report
report12 = exreport('Report for nIntervals and degree in Shenoy');

# añadir las diferentes partes del report
report12 = exreportAdd(report12, exp12);
report12 = exreportAdd(report12, test12);
report12 = exreportAdd(report12, list(plot121,plot122,plot123,plot124,table121,table121b));

# añade un resumen al report
table122 = tabularExpSummary(exp12,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report12 = exreportAdd(report12, table122);

# report en html
r12 = exreportRender(report12,target='html',visualize=T);







##########################################################################

# CONTRASTES MÚLTPLES ENTRE EL MEJOR DE CADA MODELO

##########################################################################


# TEST 13: PEREZ-BERNABE Y TMOP CON MUESTRA 1000 CON SHENOY 4 TROZOS GRADO 7

# obtener el mejor de cada modelo
m = aggregate(r$DivKL,list(r$f,r$Método),FUN=min);
a = c();
for (i in 1:nrow(m))   a[i] = which(r$DivKL==m$x[i] & r$f==m$Group.1[i] & 
                                      r$Método==m$Group.2[i])[1];

r13 = r[a,];

# ordenar la matriz de combinaciones por distribucion
r13 = r13[order(r13$f),];

# crear el objeto exreport con % de aciertos
exp13 = expCreate(r13,methods='Método',problems='f',name='Method (Shenoy, P-B, tMoPs)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m13 = aggregate(r13$DivKL,list(r13$f),FUN=min);
m131 = rep(m13$x,each=3);
r13Grafica = r13;
r13Grafica$DivKL = r13$DivKL / m131;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp13Grafica = expCreate(r13Grafica,methods='Método',problems='f',
                         name='Method (Shenoy, Pérez-B, tMoPs)');


# contraste sobre el error
test13 = testMultipleControl(exp13,output='DivKL',rankOrder='min');
test13b = testMultiplePairwise(exp13,output='DivKL',rankOrder='min');
test13Grafica = testMultipleControl(exp13Grafica,output='DivKL',rankOrder='min');

# graficas
plot131 = plotExpSummary(exp13Grafica,'DivKL');
plot132 = plotExpSummary(exp13Grafica, 'DivKL', columns=5);
plot133 = plotCumulativeRank(test13Grafica);
plot134 = plotRankDistribution(test13Grafica);

# primera table
table131 = tabularTestSummary(test13,columns=c('pvalue','rank','wtl'));
table131b = tabularTestPairwise(test13b);


# titulo del report
report13 = exreport('Report for method (Shenoy, Pérez-B, tMoPs)');

# añadir las diferentes partes del report
report13 = exreportAdd(report13, exp13);
report13 = exreportAdd(report13, test13);
report13 = exreportAdd(report13, list(plot131,plot132,plot133,plot134,table131,table131b));

# añade un resumen al report
table132 = tabularExpSummary(exp13,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report13 = exreportAdd(report13, table132);

# report en html
r13 = exreportRender(report13,target='html',visualize=T);







##########################################################################

# CONTRASTES MÚLTPLES ENTRE TODOS LOS MODELOS CON 1 TROZO

##########################################################################


# TEST 14: PEREZ-BERNABE, TMOP Y SHENOY (GRADO=3)

sh = (r$Método=='Shenoy' & r$nInter=='1-3')
pb = (r$Método=='Perez-Bernabe' & r$Grado==3 & r$Muestra==1000)
tm = (r$Método=='tMoPs' & r$Grado==3 & r$Muestra==1000)
r14 = r[sh | pb | tm,]


# crear el objeto exreport con % de aciertos
exp14 = expCreate(r14,methods='Método',problems='f',name='Method (Shenoy, P-B, tMoPs) deg=3')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m14 = aggregate(r14$DivKL,list(r14$f),FUN=min);
m141 = rep(m14$x,each=3);
r14Grafica = r14;
r14Grafica$DivKL = r14$DivKL / m141;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp14Grafica = expCreate(r14Grafica,methods='Método',problems='f',
                         name='Method (Shenoy, P-B, tMoPs) deg=3');


# contraste sobre el error
test14 = testMultipleControl(exp14,output='DivKL',rankOrder='min');
test14b = testMultiplePairwise(exp14,output='DivKL',rankOrder='min');
test14Grafica = testMultipleControl(exp14Grafica,output='DivKL',rankOrder='min');

# graficas
plot141 = plotExpSummary(exp14Grafica,'DivKL');
plot142 = plotExpSummary(exp14Grafica, 'DivKL', columns=5);
plot143 = plotCumulativeRank(test14Grafica);
plot144 = plotRankDistribution(test14Grafica);

# primera table
table141 = tabularTestSummary(test14,columns=c('pvalue','rank','wtl'));
table141b = tabularTestPairwise(test14b);


# titulo del report
report14 = exreport('Report for Method (Shenoy, P-B, tMoPs) degree=3');

# añadir las diferentes partes del report
report14 = exreportAdd(report14, exp14);
report14 = exreportAdd(report14, test14);
report14 = exreportAdd(report14, list(plot141,plot142,plot143,plot144,table141,table141b));

# añade un resumen al report
table142 = tabularExpSummary(exp14,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report14 = exreportAdd(report14, table142);

# report en html
r14 = exreportRender(report14,target='html',visualize=T);


##########################################################################




# TEST 15: PEREZ-BERNABE, TMOP Y SHENOY (GRADO=4)

sh = (r$Método=='Shenoy' & r$nInter=='1-4')
pb = (r$Método=='Perez-Bernabe' & r$Grado==4 & r$Muestra==1000)
tm = (r$Método=='tMoPs' & r$Grado==4 & r$Muestra==1000)
r15 = r[sh | pb | tm,]


# crear el objeto exreport con % de aciertos
exp15 = expCreate(r15,methods='Método',problems='f',name='Method (Shenoy, P-B, tMoPs) deg=4')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m15 = aggregate(r15$DivKL,list(r15$f),FUN=min);
m151 = rep(m15$x,each=3);
r15Grafica = r15;
r15Grafica$DivKL = r15$DivKL / m151;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp15Grafica = expCreate(r15Grafica,methods='Método',problems='f',
                         name='Method (Shenoy, P-B, tMoPs) deg=4');


# contraste sobre el error
test15 = testMultipleControl(exp15,output='DivKL',rankOrder='min');
test15b = testMultiplePairwise(exp15,output='DivKL',rankOrder='min');
test15Grafica = testMultipleControl(exp15Grafica,output='DivKL',rankOrder='min');

# graficas
plot151 = plotExpSummary(exp15Grafica,'DivKL');
plot152 = plotExpSummary(exp15Grafica, 'DivKL', columns=5);
plot153 = plotCumulativeRank(test15Grafica);
plot154 = plotRankDistribution(test15Grafica);

# primera table
table151 = tabularTestSummary(test15,columns=c('pvalue','rank','wtl'));
table151b = tabularTestPairwise(test15b);


# titulo del report
report15 = exreport('Report for Method (Shenoy, P-B, tMoPs) degree=4');

# añadir las diferentes partes del report
report15 = exreportAdd(report15, exp15);
report15 = exreportAdd(report15, test15);
report15 = exreportAdd(report15, list(plot151,plot152,plot153,plot154,table151,table151b));

# añade un resumen al report
table152 = tabularExpSummary(exp15,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report15 = exreportAdd(report15, table152);

# report en html
r15 = exreportRender(report15,target='html',visualize=T);




##########################################################################


# TEST 16: PEREZ-BERNABE, TMOP Y SHENOY (GRADO=7)

sh = (r$Método=='Shenoy' & r$nInter=='1-7')
pb = (r$Método=='Perez-Bernabe' & r$Grado==7 & r$Muestra==1000)
tm = (r$Método=='tMoPs' & r$Grado==7 & r$Muestra==1000)
r16 = r[sh | pb | tm,]


# crear el objeto exreport con % de aciertos
exp16 = expCreate(r16,methods='Método',problems='f',name='Method (Shenoy, P-B, tMoPs) deg=7')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m16 = aggregate(r16$DivKL,list(r16$f),FUN=min);
m161 = rep(m16$x,each=3);
r16Grafica = r16;
r16Grafica$DivKL = r16$DivKL / m161;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp16Grafica = expCreate(r16Grafica,methods='Método',problems='f',
                         name='Method (Shenoy, P-B, tMoPs) deg=7');


# contraste sobre el error
test16 = testMultipleControl(exp16,output='DivKL',rankOrder='min');
test16b = testMultiplePairwise(exp16,output='DivKL',rankOrder='min');
test16Grafica = testMultipleControl(exp16Grafica,output='DivKL',rankOrder='min');

# graficas
plot161 = plotExpSummary(exp16Grafica,'DivKL');
plot162 = plotExpSummary(exp16Grafica, 'DivKL', columns=5);
plot163 = plotCumulativeRank(test16Grafica);
plot164 = plotRankDistribution(test16Grafica);

# primera table
table161 = tabularTestSummary(test16,columns=c('pvalue','rank','wtl'));
table161b = tabularTestPairwise(test16b);


# titulo del report
report16 = exreport('Report for Method (Shenoy, P-B, tMoPs) degree=7');

# añadir las diferentes partes del report
report16 = exreportAdd(report16, exp16);
report16 = exreportAdd(report16, test16);
report16 = exreportAdd(report16, list(plot161,plot162,plot163,plot164,table161,table161b));

# añade un resumen al report
table162 = tabularExpSummary(exp16,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report16 = exreportAdd(report16, table162);

# report en html
r16 = exreportRender(report16,target='html',visualize=T);





##########################################################################

# CONTRASTES MÚLTPLES ENTRE TODOS LOS MODELOS CON 3 TROZOS

##########################################################################



# TEST 17: PEREZ-BERNABE, TMOP Y SHENOY (GRADO=3)

sh = (r$Método=='Shenoy' & r$nInter=='3-3')
pb = (r$Método=='Perez-Bernabe' & r$Grado==3 & r$Muestra==1000)
tm = (r$Método=='tMoPs' & r$Grado==3 & r$Muestra==1000)
r17 = r[sh | pb | tm,]


# crear el objeto exreport con % de aciertos
exp17 = expCreate(r17,methods='Método',problems='f',name='Method (Shenoy, P-B, tMoPs) deg=3')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m17 = aggregate(r17$DivKL,list(r17$f),FUN=min);
m171 = rep(m17$x,each=3);
r17Grafica = r17;
r17Grafica$DivKL = r17$DivKL / m171;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp17Grafica = expCreate(r17Grafica,methods='Método',problems='f',
                         name='Method (Shenoy, P-B, tMoPs) deg=3');


# contraste sobre el error
test17 = testMultipleControl(exp17,output='DivKL',rankOrder='min');
test17b = testMultiplePairwise(exp17,output='DivKL',rankOrder='min');
test17Grafica = testMultipleControl(exp17Grafica,output='DivKL',rankOrder='min');

# graficas
plot171 = plotExpSummary(exp17Grafica,'DivKL');
plot172 = plotExpSummary(exp17Grafica, 'DivKL', columns=5);
plot173 = plotCumulativeRank(test17Grafica);
plot174 = plotRankDistribution(test17Grafica);

# primera table
table171 = tabularTestSummary(test17,columns=c('pvalue','rank','wtl'));
table171b = tabularTestPairwise(test17b);


# titulo del report
report17 = exreport('Report for Method (Shenoy, P-B, tMoPs) degree=3');

# añadir las diferentes partes del report
report17 = exreportAdd(report17, exp17);
report17 = exreportAdd(report17, test17);
report17 = exreportAdd(report17, list(plot171,plot172,plot173,plot174,table171,table171b));

# añade un resumen al report
table172 = tabularExpSummary(exp17,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report17 = exreportAdd(report17, table172);

# report en html
r17 = exreportRender(report17,target='html',visualize=T);


##########################################################################




# TEST 18: PEREZ-BERNABE, TMOP Y SHENOY (GRADO=4)

sh = (r$Método=='Shenoy' & r$nInter=='3-4')
pb = (r$Método=='Perez-Bernabe' & r$Grado==4 & r$Muestra==1000)
tm = (r$Método=='tMoPs' & r$Grado==4 & r$Muestra==1000)
r18 = r[sh | pb | tm,]


# crear el objeto exreport con % de aciertos
exp18 = expCreate(r18,methods='Método',problems='f',name='Method (Shenoy, P-B, tMoPs) deg=4')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m18 = aggregate(r18$DivKL,list(r18$f),FUN=min);
m181 = rep(m18$x,each=3);
r18Grafica = r18;
r18Grafica$DivKL = r18$DivKL / m181;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp18Grafica = expCreate(r18Grafica,methods='Método',problems='f',
                         name='Method (Shenoy, P-B, tMoPs) deg=4');


# contraste sobre el error
test18 = testMultipleControl(exp18,output='DivKL',rankOrder='min');
test18b = testMultiplePairwise(exp18,output='DivKL',rankOrder='min');
test18Grafica = testMultipleControl(exp18Grafica,output='DivKL',rankOrder='min');

# graficas
plot181 = plotExpSummary(exp18Grafica,'DivKL');
plot182 = plotExpSummary(exp18Grafica, 'DivKL', columns=5);
plot183 = plotCumulativeRank(test18Grafica);
plot184 = plotRankDistribution(test18Grafica);

# primera table
table181 = tabularTestSummary(test18,columns=c('pvalue','rank','wtl'));
table181b = tabularTestPairwise(test18b);


# titulo del report
report18 = exreport('Report for Method (Shenoy, P-B, tMoPs) degree=4');

# añadir las diferentes partes del report
report18 = exreportAdd(report18, exp18);
report18 = exreportAdd(report18, test18);
report18 = exreportAdd(report18, list(plot181,plot182,plot183,plot184,table181,table181b));

# añade un resumen al report
table182 = tabularExpSummary(exp18,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report18 = exreportAdd(report18, table182);

# report en html
r18 = exreportRender(report18,target='html',visualize=T);


##########################################################################



# TEST 19: PEREZ-BERNABE, TMOP Y SHENOY (GRADO=7)

sh = (r$Método=='Shenoy' & r$nInter=='3-7')
pb = (r$Método=='Perez-Bernabe' & r$Grado==7 & r$Muestra==1000)
tm = (r$Método=='tMoPs' & r$Grado==7 & r$Muestra==1000)
r19 = r[sh | pb | tm,]


# crear el objeto exreport con % de aciertos
exp19 = expCreate(r19,methods='Método',problems='f',name='Method (Shenoy, P-B, tMoPs) deg=7')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m19 = aggregate(r19$DivKL,list(r19$f),FUN=min);
m191 = rep(m19$x,each=3);
r19Grafica = r19;
r19Grafica$DivKL = r19$DivKL / m191;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp19Grafica = expCreate(r19Grafica,methods='Método',problems='f',
                         name='Method (Shenoy, P-B, tMoPs) deg=7');


# contraste sobre el error
test19 = testMultipleControl(exp19,output='DivKL',rankOrder='min');
test19b = testMultiplePairwise(exp19,output='DivKL',rankOrder='min');
test19Grafica = testMultipleControl(exp19Grafica,output='DivKL',rankOrder='min');

# graficas
plot191 = plotExpSummary(exp19Grafica,'DivKL');
plot192 = plotExpSummary(exp19Grafica, 'DivKL', columns=5);
plot193 = plotCumulativeRank(test19Grafica);
plot194 = plotRankDistribution(test19Grafica);

# primera table
table191 = tabularTestSummary(test19,columns=c('pvalue','rank','wtl'));
table191b = tabularTestPairwise(test19b);


# titulo del report
report19 = exreport('Report for Method (Shenoy, P-B, tMoPs) degree=7');

# añadir las diferentes partes del report
report19 = exreportAdd(report19, exp19);
report19 = exreportAdd(report19, test19);
report19 = exreportAdd(report19, list(plot191,plot192,plot193,plot194,table191,table191b));

# añade un resumen al report
table192 = tabularExpSummary(exp19,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report19 = exreportAdd(report19, table192);

# report en html
r19 = exreportRender(report19,target='html',visualize=T);










##########################################################################

# CONTRASTES 2 A 2 EN TAMAÑO DE MUESTRA EN TMOP

##########################################################################


# TEST 30: TMOP 50 vs 100  (GRADO=3)

r30 = r[r$Método=='tMoPs' & r$Grado==3 & r$Muestra!=1000,];

# crear el objeto exreport con % de aciertos
exp30 = expCreate(r30,methods='Muestra',problems='f',
                  name='Sample size in tMoPs 50 vs 100 (deg=3)')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m30 = aggregate(r30$DivKL,list(r30$f),FUN=min);
m301 = rep(m30$x,each=2);
r30Grafica = r30;
r30Grafica$DivKL = r30$DivKL / m301;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp30Grafica = expCreate(r30Grafica,methods='Muestra',problems='f',
                         name='Sample size in tMoPs 50 vs 100 (deg=3)');

# contraste sobre el error
test30 = testPaired(exp30,output='DivKL',rankOrder='min');

# graficas
plot301 = plotExpSummary(exp30Grafica,'DivKL');
plot302 = plotExpSummary(exp30Grafica,'DivKL',columns=5);


# titulo del report
report30 = exreport('Report for sample size in tMoPs 50 vs 100 (deg=3)');

# añadir las diferentes partes del report
report30 = exreportAdd(report30, exp30);
report30 = exreportAdd(report30, test30);
report30 = exreportAdd(report30, list(plot301,plot302));

# añade un resumen al report
table30 = tabularExpSummary(exp30,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report30 = exreportAdd(report30, table30);

# report en html
r30 = exreportRender(report30,target='html',visualize=T);


##########################################################################



# TEST 31: TMOP 50 vs 100  (GRADO=4)

r31 = r[r$Método=='tMoPs' & r$Grado==4 & r$Muestra!=1000,];

# crear el objeto exreport con % de aciertos
exp31 = expCreate(r31,methods='Muestra',problems='f',
                  name='Sample size in tMoPs 50 vs 100 (deg=4)')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m31 = aggregate(r31$DivKL,list(r31$f),FUN=min);
m311 = rep(m31$x,each=2);
r31Grafica = r31;
r31Grafica$DivKL = r31$DivKL / m311;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp31Grafica = expCreate(r31Grafica,methods='Muestra',problems='f',
                         name='Sample size in tMoPs 50 vs 100 (deg=4)');

# contraste sobre el error
test31 = testPaired(exp31,output='DivKL',rankOrder='min');

# graficas
plot311 = plotExpSummary(exp31Grafica,'DivKL');
plot312 = plotExpSummary(exp31Grafica,'DivKL',columns=5);


# titulo del report
report31 = exreport('Report for sample size in tMoPs 50 vs 100 (deg=4)');

# añadir las diferentes partes del report
report31 = exreportAdd(report31, exp31);
report31 = exreportAdd(report31, test31);
report31 = exreportAdd(report31, list(plot311,plot312));

# añade un resumen al report
table31 = tabularExpSummary(exp31,'DivKL',digits=4,format="f",boldfaceColumns='min',
                             tableSplit=2);
report31 = exreportAdd(report31, table31);

# report en html
r31 = exreportRender(report31,target='html',visualize=T);


##########################################################################



# TEST 32: TMOP 50 vs 100  (GRADO=7)

r32 = r[r$Método=='tMoPs' & r$Grado==7 & r$Muestra!=1000,];

# crear el objeto exreport con % de aciertos
exp32 = expCreate(r32,methods='Muestra',problems='f',
                  name='Sample size in tMoPs 50 vs 100 (deg=7)')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m32 = aggregate(r32$DivKL,list(r32$f),FUN=min);
m321 = rep(m32$x,each=2);
r32Grafica = r32;
r32Grafica$DivKL = r32$DivKL / m321;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp32Grafica = expCreate(r32Grafica,methods='Muestra',problems='f',
                         name='Sample size in tMoPs 50 vs 100 (deg=7)');

# contraste sobre el error
test32 = testPaired(exp32,output='DivKL',rankOrder='min');

# graficas
plot321 = plotExpSummary(exp32Grafica,'DivKL');
plot322 = plotExpSummary(exp32Grafica,'DivKL',columns=5);


# titulo del report
report32 = exreport('Report for sample size in tMoPs 50 vs 100 (deg=7)');

# añadir las diferentes partes del report
report32 = exreportAdd(report32, exp32);
report32 = exreportAdd(report32, test32);
report32 = exreportAdd(report32, list(plot321,plot322));

# añade un resumen al report
table32 = tabularExpSummary(exp32,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report32 = exreportAdd(report32, table32);

# report en html
r32 = exreportRender(report32,target='html',visualize=T);



##########################################################################


# TEST 33: TMOP 100 vs 1000  (GRADO=3)

r33 = r[r$Método=='tMoPs' & r$Grado==3 & r$Muestra!=50,];

# crear el objeto exreport con % de aciertos
exp33 = expCreate(r33,methods='Muestra',problems='f',
                  name='Sample size in tMoPs 100 vs 1000 (deg=3)')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m33 = aggregate(r33$DivKL,list(r33$f),FUN=min);
m331 = rep(m33$x,each=2);
r33Grafica = r33;
r33Grafica$DivKL = r33$DivKL / m331;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp33Grafica = expCreate(r33Grafica,methods='Muestra',problems='f',
                         name='Sample size in tMoPs 100 vs 1000 (deg=3)');

# contraste sobre el error
test33 = testPaired(exp33,output='DivKL',rankOrder='min');

# graficas
plot331 = plotExpSummary(exp33Grafica,'DivKL');
plot332 = plotExpSummary(exp33Grafica,'DivKL',columns=5);


# titulo del report
report33 = exreport('Report for sample size in tMoPs 100 vs 1000 (deg=3)');

# añadir las diferentes partes del report
report33 = exreportAdd(report33, exp33);
report33 = exreportAdd(report33, test33);
report33 = exreportAdd(report33, list(plot331,plot332));

# añade un resumen al report
table33 = tabularExpSummary(exp33,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report33 = exreportAdd(report33, table33);

# report en html
r33 = exreportRender(report33,target='html',visualize=T);



##########################################################################


# TEST 34: TMOP 100 vs 1000  (GRADO=4)

r34 = r[r$Método=='tMoPs' & r$Grado==4 & r$Muestra!=50,];

# crear el objeto exreport con % de aciertos
exp34 = expCreate(r34,methods='Muestra',problems='f',
                  name='Sample size in tMoPs 100 vs 1000 (deg=4)')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m34 = aggregate(r34$DivKL,list(r34$f),FUN=min);
m341 = rep(m34$x,each=2);
r34Grafica = r34;
r34Grafica$DivKL = r34$DivKL / m341;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp34Grafica = expCreate(r34Grafica,methods='Muestra',problems='f',
                         name='Sample size in tMoPs 100 vs 1000 (deg=4)');

# contraste sobre el error
test34 = testPaired(exp34,output='DivKL',rankOrder='min');

# graficas
plot341 = plotExpSummary(exp34Grafica,'DivKL');
plot342 = plotExpSummary(exp34Grafica,'DivKL',columns=5);


# titulo del report
report34 = exreport('Report for sample size in tMoPs 100 vs 1000 (deg=4)');

# añadir las diferentes partes del report
report34 = exreportAdd(report34, exp34);
report34 = exreportAdd(report34, test34);
report34 = exreportAdd(report34, list(plot341,plot342));

# añade un resumen al report
table34 = tabularExpSummary(exp34,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report34 = exreportAdd(report34, table34);

# report en html
r34 = exreportRender(report34,target='html',visualize=T);



##########################################################################


# TEST 35: TMOP 100 vs 1000  (GRADO=7)

r35 = r[r$Método=='tMoPs' & r$Grado==7 & r$Muestra!=50,];

# crear el objeto exreport con % de aciertos
exp35 = expCreate(r35,methods='Muestra',problems='f',
                  name='Sample size in tMoPs 100 vs 1000 (deg=7)')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m35 = aggregate(r35$DivKL,list(r35$f),FUN=min);
m351 = rep(m35$x,each=2);
r35Grafica = r35;
r35Grafica$DivKL = r35$DivKL / m351;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp35Grafica = expCreate(r35Grafica,methods='Muestra',problems='f',
                         name='Sample size in tMoPs 100 vs 1000 (deg=7)');

# contraste sobre el error
test35 = testPaired(exp35,output='DivKL',rankOrder='min');

# graficas
plot351 = plotExpSummary(exp35Grafica,'DivKL');
plot352 = plotExpSummary(exp35Grafica,'DivKL',columns=5);


# titulo del report
report35 = exreport('Report for sample size in tMoPs 100 vs 1000 (deg=7)');

# añadir las diferentes partes del report
report35 = exreportAdd(report35, exp35);
report35 = exreportAdd(report35, test35);
report35 = exreportAdd(report35, list(plot351,plot352));

# añade un resumen al report
table35 = tabularExpSummary(exp35,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report35 = exreportAdd(report35, table35);

# report en html
r35 = exreportRender(report35,target='html',visualize=T);




##########################################################################


# TEST 36: TMOP 50 vs 1000  (GRADO=3)

r36 = r[r$Método=='tMoPs' & r$Grado==3 & r$Muestra!=100,];

# crear el objeto exreport con % de aciertos
exp36 = expCreate(r36,methods='Muestra',problems='f',
                  name='Sample size in tMoPs 50 vs 1000 (deg=3)')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m36 = aggregate(r36$DivKL,list(r36$f),FUN=min);
m361 = rep(m36$x,each=2);
r36Grafica = r36;
r36Grafica$DivKL = r36$DivKL / m361;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp36Grafica = expCreate(r36Grafica,methods='Muestra',problems='f',
                         name='Sample size in tMoPs 50 vs 1000 (deg=3)');

# contraste sobre el error
test36 = testPaired(exp36,output='DivKL',rankOrder='min');

# graficas
plot361 = plotExpSummary(exp36Grafica,'DivKL');
plot362 = plotExpSummary(exp36Grafica,'DivKL',columns=5);


# titulo del report
report36 = exreport('Report for sample size in tMoPs 50 vs 1000 (deg=3)');

# añadir las diferentes partes del report
report36 = exreportAdd(report36, exp36);
report36 = exreportAdd(report36, test36);
report36 = exreportAdd(report36, list(plot361,plot362));

# añade un resumen al report
table36 = tabularExpSummary(exp36,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report36 = exreportAdd(report36, table36);

# report en html
r36 = exreportRender(report36,target='html',visualize=T);



##########################################################################


# TEST 37: TMOP 50 vs 1000  (GRADO=4)

r37 = r[r$Método=='tMoPs' & r$Grado==4 & r$Muestra!=100,];

# crear el objeto exreport con % de aciertos
exp37 = expCreate(r37,methods='Muestra',problems='f',
                  name='Sample size in tMoPs 50 vs 1000 (deg=4)')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m37 = aggregate(r37$DivKL,list(r37$f),FUN=min);
m371 = rep(m37$x,each=2);
r37Grafica = r37;
r37Grafica$DivKL = r37$DivKL / m371;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp37Grafica = expCreate(r37Grafica,methods='Muestra',problems='f',
                         name='Sample size in tMoPs 50 vs 1000 (deg=4)');

# contraste sobre el error
test37 = testPaired(exp37,output='DivKL',rankOrder='min');

# graficas
plot371 = plotExpSummary(exp37Grafica,'DivKL');
plot372 = plotExpSummary(exp37Grafica,'DivKL',columns=5);


# titulo del report
report37 = exreport('Report for sample size in tMoPs 50 vs 1000 (deg=4)');

# añadir las diferentes partes del report
report37 = exreportAdd(report37, exp37);
report37 = exreportAdd(report37, test37);
report37 = exreportAdd(report37, list(plot371,plot372));

# añade un resumen al report
table37 = tabularExpSummary(exp37,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report37 = exreportAdd(report37, table37);

# report en html
r37 = exreportRender(report37,target='html',visualize=T);



##########################################################################


# TEST 38: TMOP 50 vs 1000  (GRADO=7)

r38 = r[r$Método=='tMoPs' & r$Grado==7 & r$Muestra!=100,];

# crear el objeto exreport con % de aciertos
exp38 = expCreate(r38,methods='Muestra',problems='f',
                  name='Sample size in tMoPs 50 vs 1000 (deg=7)')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m38 = aggregate(r38$DivKL,list(r38$f),FUN=min);
m381 = rep(m38$x,each=2);
r38Grafica = r38;
r38Grafica$DivKL = r38$DivKL / m381;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp38Grafica = expCreate(r38Grafica,methods='Muestra',problems='f',
                         name='Sample size in tMoPs 50 vs 1000 (deg=7)');

# contraste sobre el error
test38 = testPaired(exp38,output='DivKL',rankOrder='min');

# graficas
plot381 = plotExpSummary(exp38Grafica,'DivKL');
plot382 = plotExpSummary(exp38Grafica,'DivKL',columns=5);


# titulo del report
report38 = exreport('Report for sample size in tMoPs 50 vs 1000 (deg=7)');

# añadir las diferentes partes del report
report38 = exreportAdd(report38, exp38);
report38 = exreportAdd(report38, test38);
report38 = exreportAdd(report38, list(plot381,plot382));

# añade un resumen al report
table38 = tabularExpSummary(exp38,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report38 = exreportAdd(report38, table38);

# report en html
r38 = exreportRender(report38,target='html',visualize=T);








##########################################################################

# CONTRASTES 2 A 2 EN GRADO EN TMOP CON TAMAÑO DE MUESTRA=1000 (IGUALDAD EN LOS OTROS)

##########################################################################


# TEST 40: TMOP 3 vs 4

r40 = r[r$Método=='tMoPs' & r$Grado!=7 & r$Muestra==1000,];

# crear el objeto exreport con % de aciertos
exp40 = expCreate(r40,methods='Grado',problems='f',
                  name='Degree in tMoPs 3 vs 4 (sample=1000)')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m40 = aggregate(r40$DivKL,list(r40$f),FUN=min);
m401 = rep(m40$x,each=2);
r40Grafica = r40;
r40Grafica$DivKL = r40$DivKL / m401;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp40Grafica = expCreate(r40Grafica,methods='Grado',problems='f',
                         name='Degree in tMoPs 3 vs 4 (sample=1000)');

# contraste sobre el error
test40 = testPaired(exp40,output='DivKL',rankOrder='min');

# graficas
plot401 = plotExpSummary(exp40Grafica,'DivKL');
plot402 = plotExpSummary(exp40Grafica,'DivKL',columns=5);


# titulo del report
report40 = exreport('Report for degree in tMoPs 3 vs 4 (sample=1000)');

# añadir las diferentes partes del report
report40 = exreportAdd(report40, exp40);
report40 = exreportAdd(report40, test40);
report40 = exreportAdd(report40, list(plot401,plot402));

# añade un resumen al report
table40 = tabularExpSummary(exp40,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report40 = exreportAdd(report40, table40);

# report en html
r40 = exreportRender(report40,target='html',visualize=T);


##########################################################################



# TEST 41: TMOP 3 vs 7

r41 = r[r$Método=='tMoPs' & r$Grado!=4 & r$Muestra==1000,];

# crear el objeto exreport con % de aciertos
exp41 = expCreate(r41,methods='Grado',problems='f',
                  name='Degree in tMoPs 3 vs 7 (sample=1000)')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m41 = aggregate(r41$DivKL,list(r41$f),FUN=min);
m411 = rep(m41$x,each=2);
r41Grafica = r41;
r41Grafica$DivKL = r41$DivKL / m411;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp41Grafica = expCreate(r41Grafica,methods='Grado',problems='f',
                         name='Degree in tMoPs 3 vs 7 (sample=1000)');

# contraste sobre el error
test41 = testPaired(exp41,output='DivKL',rankOrder='min');

# graficas
plot411 = plotExpSummary(exp41Grafica,'DivKL');
plot412 = plotExpSummary(exp41Grafica,'DivKL',columns=5);


# titulo del report
report41 = exreport('Report for degree in tMoPs 3 vs 7 (sample=1000)');

# añadir las diferentes partes del report
report41 = exreportAdd(report41, exp41);
report41 = exreportAdd(report41, test41);
report41 = exreportAdd(report41, list(plot411,plot412));

# añade un resumen al report
table41 = tabularExpSummary(exp41,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report41 = exreportAdd(report41, table41);

# report en html
r41 = exreportRender(report41,target='html',visualize=T);


##########################################################################



# TEST 42: TMOP 4 vs 7

r42 = r[r$Método=='tMoPs' & r$Grado!=3 & r$Muestra==1000,];

# crear el objeto exreport con % de aciertos
exp42 = expCreate(r42,methods='Grado',problems='f',
                  name='Degree in tMoPs 4 vs 7 (sample=1000)')

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m42 = aggregate(r42$DivKL,list(r42$f),FUN=min);
m421 = rep(m42$x,each=2);
r42Grafica = r42;
r42Grafica$DivKL = r42$DivKL / m421;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp42Grafica = expCreate(r42Grafica,methods='Grado',problems='f',
                         name='Degree in tMoPs 4 vs 7 (sample=1000)');

# contraste sobre el error
test42 = testPaired(exp42,output='DivKL',rankOrder='min');

# graficas
plot421 = plotExpSummary(exp42Grafica,'DivKL');
plot422 = plotExpSummary(exp42Grafica,'DivKL',columns=5);


# titulo del report
report42 = exreport('Report for degree in tMoPs 4 vs 7 (sample=1000)');

# añadir las diferentes partes del report
report42 = exreportAdd(report42, exp42);
report42 = exreportAdd(report42, test42);
report42 = exreportAdd(report42, list(plot421,plot422));

# añade un resumen al report
table42 = tabularExpSummary(exp42,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report42 = exreportAdd(report42, table42);

# report en html
r42 = exreportRender(report42,target='html',visualize=T);






##########################################################################

# CONTRASTES 2 A 2 EN MEJOR DE CADA MÉTODO (SHENOY, P-B, tMoPs)

##########################################################################



# TEST 50: PEREZ-BERNABE vs tMoPs

# obtener el mejor de cada modelo
m = aggregate(r$DivKL,list(r$f,r$Método),FUN=min);
a = c();
for (i in 1:nrow(m))   a[i] = which(r$DivKL==m$x[i] & r$f==m$Group.1[i] & 
                                    r$Método==m$Group.2[i])[1];

r50 = r[a,];
r50 = r50[r50$Método!='Shenoy',];

# ordenar la matriz de combinaciones por distribucion
r50 = r50[order(r50$f),];

# crear el objeto exreport con % de aciertos
exp50 = expCreate(r50,methods='Método',problems='f',name='Method P-B vs tMoPs (best)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m50 = aggregate(r50$DivKL,list(r50$f),FUN=min);
m501 = rep(m50$x,each=2);
r50Grafica = r50;
r50Grafica$DivKL = r50$DivKL / m501;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp50Grafica = expCreate(r50Grafica,methods='Método',problems='f',
                         name='Method P-B vs tMoPs (best)');

# contraste sobre el error
test50 = testPaired(exp50,output='DivKL',rankOrder='min');

# graficas
plot501 = plotExpSummary(exp50Grafica,'DivKL');
plot502 = plotExpSummary(exp50Grafica,'DivKL',columns=5);


# titulo del report
report50 = exreport('Report for method P-B vs tMoPs (best)');

# añadir las diferentes partes del report
report50 = exreportAdd(report50, exp50);
report50 = exreportAdd(report50, test50);
report50 = exreportAdd(report50, list(plot501,plot502));

# añade un resumen al report
table50 = tabularExpSummary(exp50,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report50 = exreportAdd(report50, table50);

# report en html
r50 = exreportRender(report50,target='html',visualize=T);


##########################################################################



# TEST 51: PEREZ-BERNABE vs SHENOY

# obtener el mejor de cada modelo
m = aggregate(r$DivKL,list(r$f,r$Método),FUN=min);
a = c();
for (i in 1:nrow(m))   a[i] = which(r$DivKL==m$x[i] & r$f==m$Group.1[i] & 
                                      r$Método==m$Group.2[i])[1];

r51 = r[a,];
r51 = r51[r51$Método!='tMoPs',];

# ordenar la matriz de combinaciones por distribucion
r51 = r51[order(r51$f),];

# crear el objeto exreport con % de aciertos
exp51 = expCreate(r51,methods='Método',problems='f',name='Method P-B vs Shenoy (best)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m51 = aggregate(r51$DivKL,list(r51$f),FUN=min);
m511 = rep(m51$x,each=2);
r51Grafica = r51;
r51Grafica$DivKL = r51$DivKL / m511;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp51Grafica = expCreate(r51Grafica,methods='Método',problems='f',
                         name='Method P-B vs Shenoy (best)');

# contraste sobre el error
test51 = testPaired(exp51,output='DivKL',rankOrder='min');

# graficas
plot511 = plotExpSummary(exp51Grafica,'DivKL');
plot512 = plotExpSummary(exp51Grafica,'DivKL',columns=5);


# titulo del report
report51 = exreport('Report for method P-B vs Shenoy (best)');

# añadir las diferentes partes del report
report51 = exreportAdd(report51, exp51);
report51 = exreportAdd(report51, test51);
report51 = exreportAdd(report51, list(plot511,plot512));

# añade un resumen al report
table51 = tabularExpSummary(exp51,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report51 = exreportAdd(report51, table51);

# report en html
r51 = exreportRender(report51,target='html',visualize=T);


##########################################################################



# TEST 52: tMoPs vs SHENOY

# obtener el mejor de cada modelo
m = aggregate(r$DivKL,list(r$f,r$Método),FUN=min);
a = c();
for (i in 1:nrow(m))   a[i] = which(r$DivKL==m$x[i] & r$f==m$Group.1[i] & 
                                      r$Método==m$Group.2[i])[1];

r52 = r[a,];
r52 = r52[r52$Método!='Perez-Bernabe',];

# ordenar la matriz de combinaciones por distribucion
r52 = r52[order(r52$f),];

# crear el objeto exreport con % de aciertos
exp52 = expCreate(r52,methods='Método',problems='f',name='Method tMoPs vs Shenoy (best)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m52 = aggregate(r52$DivKL,list(r52$f),FUN=min);
m521 = rep(m52$x,each=2);
r52Grafica = r52;
r52Grafica$DivKL = r52$DivKL / m521;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp52Grafica = expCreate(r52Grafica,methods='Método',problems='f',
                         name='Method tMoPs vs Shenoy (best)');

# contraste sobre el error
test52 = testPaired(exp52,output='DivKL',rankOrder='min');

# graficas
plot521 = plotExpSummary(exp52Grafica,'DivKL');
plot522 = plotExpSummary(exp52Grafica,'DivKL',columns=5);


# titulo del report
report52 = exreport('Report for method tMoPs vs Shenoy (best)');

# añadir las diferentes partes del report
report52 = exreportAdd(report52, exp52);
report52 = exreportAdd(report52, test52);
report52 = exreportAdd(report52, list(plot521,plot522));

# añade un resumen al report
table52 = tabularExpSummary(exp52,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report52 = exreportAdd(report52, table52);

# report en html
r52 = exreportRender(report52,target='html',visualize=T);





##########################################################################

# CONTRASTES 2 A 2 CON 1 TROZO, MUESTRA=1000 Y GRADO=3 (SHENOY, P-B, tMoPs)

##########################################################################



# TEST 53: PEREZ-BERNABE vs tMoPs

r53 = r[r$Método!='Shenoy' & r$Grado==3 & r$Muestra==1000,]

# crear el objeto exreport con % de aciertos
exp53 = expCreate(r53,methods='Método',problems='f',
                  name='Method P-B vs tMoPs (1 piece, deg=3, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m53 = aggregate(r53$DivKL,list(r53$f),FUN=min);
m531 = rep(m53$x,each=2);
r53Grafica = r53;
r53Grafica$DivKL = r53$DivKL / m531;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp53Grafica = expCreate(r53Grafica,methods='Método',problems='f',
                         name='Method P-B vs tMoPs (1 piece, deg=3, sample=1000)');

# contraste sobre el error
test53 = testPaired(exp53,output='DivKL',rankOrder='min');

# graficas
plot531 = plotExpSummary(exp53Grafica,'DivKL');
plot532 = plotExpSummary(exp53Grafica,'DivKL',columns=5);


# titulo del report
report53 = exreport('Report for method P-B vs tMoPs (1 piece, deg=3, sample=1000)');

# añadir las diferentes partes del report
report53 = exreportAdd(report53, exp53);
report53 = exreportAdd(report53, test53);
report53 = exreportAdd(report53, list(plot531,plot532));

# añade un resumen al report
table53 = tabularExpSummary(exp53,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report53 = exreportAdd(report53, table53);

# report en html
r53 = exreportRender(report53,target='html',visualize=T);


##########################################################################




# TEST 54: PEREZ-BERNABE vs SHENOY

sh = (r$Método=='Shenoy' & r$nInter=='1-3')
pb = (r$Método=='Perez-Bernabe' & r$Grado==3 & r$Muestra==1000)
r54 = r[sh | pb,]

# crear el objeto exreport con % de aciertos
exp54 = expCreate(r54,methods='Método',problems='f',
                  name='Method P-B vs Shenoy (1 piece, deg=3, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m54 = aggregate(r54$DivKL,list(r54$f),FUN=min);
m541 = rep(m54$x,each=2);
r54Grafica = r54;
r54Grafica$DivKL = r54$DivKL / m541;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp54Grafica = expCreate(r54Grafica,methods='Método',problems='f',
                         name='Method P-B vs Shenoy (1 piece, deg=3, sample=1000)');

# contraste sobre el error
test54 = testPaired(exp54,output='DivKL',rankOrder='min');

# graficas
plot541 = plotExpSummary(exp54Grafica,'DivKL');
plot542 = plotExpSummary(exp54Grafica,'DivKL',columns=5);


# titulo del report
report54 = exreport('Report for method P-B vs Shenoy (1 piece, deg=3, sample=1000)');

# añadir las diferentes partes del report
report54 = exreportAdd(report54, exp54);
report54 = exreportAdd(report54, test54);
report54 = exreportAdd(report54, list(plot541,plot542));

# añade un resumen al report
table54 = tabularExpSummary(exp54,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report54 = exreportAdd(report54, table54);

# report en html
r54 = exreportRender(report54,target='html',visualize=T);


##########################################################################




# TEST 55: tMoPs vs SHENOY

sh = (r$Método=='Shenoy' & r$nInter=='1-3')
tm = (r$Método=='tMoPs' & r$Grado==3 & r$Muestra==1000)
r55 = r[sh | tm,]

# crear el objeto exreport con % de aciertos
exp55 = expCreate(r55,methods='Método',problems='f',
                  name='Method tMoPs vs Shenoy (1 piece, deg=3, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m55 = aggregate(r55$DivKL,list(r55$f),FUN=min);
m551 = rep(m55$x,each=2);
r55Grafica = r55;
r55Grafica$DivKL = r55$DivKL / m551;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp55Grafica = expCreate(r55Grafica,methods='Método',problems='f',
                         name='Method tMoPs vs Shenoy (1 piece, deg=3, sample=1000)');

# contraste sobre el error
test55 = testPaired(exp55,output='DivKL',rankOrder='min');

# graficas
plot551 = plotExpSummary(exp55Grafica,'DivKL');
plot552 = plotExpSummary(exp55Grafica,'DivKL',columns=5);


# titulo del report
report55 = exreport('Report for method tMoPs vs Shenoy (1 piece, deg=3, sample=1000)');

# añadir las diferentes partes del report
report55 = exreportAdd(report55, exp55);
report55 = exreportAdd(report55, test55);
report55 = exreportAdd(report55, list(plot551,plot552));

# añade un resumen al report
table55 = tabularExpSummary(exp55,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report55 = exreportAdd(report55, table55);

# report en html
r55 = exreportRender(report55,target='html',visualize=T);






##########################################################################

# CONTRASTES 2 A 2 CON 1 TROZO, MUESTRA=1000 Y GRADO=4 (SHENOY, P-B, tMoPs)

##########################################################################



# TEST 56: PEREZ-BERNABE vs tMoPs

r56 = r[r$Método!='Shenoy' & r$Grado==4 & r$Muestra==1000,]

# crear el objeto exreport con % de aciertos
exp56 = expCreate(r56,methods='Método',problems='f',
                  name='Method P-B vs tMoPs (1 piece, deg=4, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m56 = aggregate(r56$DivKL,list(r56$f),FUN=min);
m561 = rep(m56$x,each=2);
r56Grafica = r56;
r56Grafica$DivKL = r56$DivKL / m561;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp56Grafica = expCreate(r56Grafica,methods='Método',problems='f',
                         name='Method P-B vs tMoPs (1 piece, deg=4, sample=1000)');

# contraste sobre el error
test56 = testPaired(exp56,output='DivKL',rankOrder='min');

# graficas
plot561 = plotExpSummary(exp56Grafica,'DivKL');
plot562 = plotExpSummary(exp56Grafica,'DivKL',columns=5);


# titulo del report
report56 = exreport('Report for method P-B vs tMoPs (1 piece, deg=4, sample=1000)');

# añadir las diferentes partes del report
report56 = exreportAdd(report56, exp56);
report56 = exreportAdd(report56, test56);
report56 = exreportAdd(report56, list(plot561,plot562));

# añade un resumen al report
table56 = tabularExpSummary(exp56,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report56 = exreportAdd(report56, table56);

# report en html
r56 = exreportRender(report56,target='html',visualize=T);


##########################################################################





# TEST 57: PEREZ-BERNABE vs SHENOY

sh = (r$Método=='Shenoy' & r$nInter=='1-4')
pb = (r$Método=='Perez-Bernabe' & r$Grado==4 & r$Muestra==1000)
r57 = r[sh | pb,]

# crear el objeto exreport con % de aciertos
exp57 = expCreate(r57,methods='Método',problems='f',
                  name='Method P-B vs Shenoy (1 piece, deg=4, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m57 = aggregate(r57$DivKL,list(r57$f),FUN=min);
m571 = rep(m57$x,each=2);
r57Grafica = r57;
r57Grafica$DivKL = r57$DivKL / m571;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp57Grafica = expCreate(r57Grafica,methods='Método',problems='f',
                         name='Method P-B vs Shenoy (1 piece, deg=4, sample=1000)');

# contraste sobre el error
test57 = testPaired(exp57,output='DivKL',rankOrder='min');

# graficas
plot571 = plotExpSummary(exp57Grafica,'DivKL');
plot572 = plotExpSummary(exp57Grafica,'DivKL',columns=5);


# titulo del report
report57 = exreport('Report for method P-B vs Shenoy (1 piece, deg=4, sample=1000)');

# añadir las diferentes partes del report
report57 = exreportAdd(report57, exp57);
report57 = exreportAdd(report57, test57);
report57 = exreportAdd(report57, list(plot571,plot572));

# añade un resumen al report
table57 = tabularExpSummary(exp57,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report57 = exreportAdd(report57, table57);

# report en html
r57 = exreportRender(report57,target='html',visualize=T);


##########################################################################




# TEST 58: tMoPs vs SHENOY

sh = (r$Método=='Shenoy' & r$nInter=='1-4')
tm = (r$Método=='tMoPs' & r$Grado==4 & r$Muestra==1000)
r58 = r[sh | tm,]

# crear el objeto exreport con % de aciertos
exp58 = expCreate(r58,methods='Método',problems='f',
                  name='Method tMoPs vs Shenoy (1 piece, deg=4, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m58 = aggregate(r58$DivKL,list(r58$f),FUN=min);
m581 = rep(m58$x,each=2);
r58Grafica = r58;
r58Grafica$DivKL = r58$DivKL / m581;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp58Grafica = expCreate(r58Grafica,methods='Método',problems='f',
                         name='Method tMoPs vs Shenoy (1 piece, deg=4, sample=1000)');

# contraste sobre el error
test58 = testPaired(exp58,output='DivKL',rankOrder='min');

# graficas
plot581 = plotExpSummary(exp58Grafica,'DivKL');
plot582 = plotExpSummary(exp58Grafica,'DivKL',columns=5);


# titulo del report
report58 = exreport('Report for method tMoPs vs Shenoy (1 piece, deg=4, sample=1000)');

# añadir las diferentes partes del report
report58 = exreportAdd(report58, exp58);
report58 = exreportAdd(report58, test58);
report58 = exreportAdd(report58, list(plot581,plot582));

# añade un resumen al report
table58 = tabularExpSummary(exp58,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report58 = exreportAdd(report58, table58);

# report en html
r58 = exreportRender(report58,target='html',visualize=T);








##########################################################################

# CONTRASTES 2 A 2 CON 1 TROZO, MUESTRA=1000 Y GRADO=7 (SHENOY, P-B, tMoPs)

##########################################################################



# TEST 59: PEREZ-BERNABE vs tMoPs

r59 = r[r$Método!='Shenoy' & r$Grado==7 & r$Muestra==1000,]

# crear el objeto exreport con % de aciertos
exp59 = expCreate(r59,methods='Método',problems='f',
                  name='Method P-B vs tMoPs (1 piece, deg=7, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m59 = aggregate(r59$DivKL,list(r59$f),FUN=min);
m591 = rep(m59$x,each=2);
r59Grafica = r59;
r59Grafica$DivKL = r59$DivKL / m591;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp59Grafica = expCreate(r59Grafica,methods='Método',problems='f',
                         name='Method P-B vs tMoPs (1 piece, deg=7, sample=1000)');

# contraste sobre el error
test59 = testPaired(exp59,output='DivKL',rankOrder='min');

# graficas
plot591 = plotExpSummary(exp59Grafica,'DivKL');
plot592 = plotExpSummary(exp59Grafica,'DivKL',columns=5);


# titulo del report
report59 = exreport('Report for method P-B vs tMoPs (1 piece, deg=7, sample=1000)');

# añadir las diferentes partes del report
report59 = exreportAdd(report59, exp59);
report59 = exreportAdd(report59, test59);
report59 = exreportAdd(report59, list(plot591,plot592));

# añade un resumen al report
table59 = tabularExpSummary(exp59,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report59 = exreportAdd(report59, table59);

# report en html
r59 = exreportRender(report59,target='html',visualize=T);


##########################################################################





# TEST 60: PEREZ-BERNABE vs SHENOY

sh = (r$Método=='Shenoy' & r$nInter=='1-7')
pb = (r$Método=='Perez-Bernabe' & r$Grado==7 & r$Muestra==1000)
r60 = r[sh | pb,]

# crear el objeto exreport con % de aciertos
exp60 = expCreate(r60,methods='Método',problems='f',
                  name='Method P-B vs Shenoy (1 piece, deg=7, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m60 = aggregate(r60$DivKL,list(r60$f),FUN=min);
m601 = rep(m60$x,each=2);
r60Grafica = r60;
r60Grafica$DivKL = r60$DivKL / m601;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp60Grafica = expCreate(r60Grafica,methods='Método',problems='f',
                         name='Method P-B vs Shenoy (1 piece, deg=7, sample=1000)');

# contraste sobre el error
test60 = testPaired(exp60,output='DivKL',rankOrder='min');

# graficas
plot601 = plotExpSummary(exp60Grafica,'DivKL');
plot602 = plotExpSummary(exp60Grafica,'DivKL',columns=5);


# titulo del report
report60 = exreport('Report for method P-B vs Shenoy (1 piece, deg=7, sample=1000)');

# añadir las diferentes partes del report
report60 = exreportAdd(report60, exp60);
report60 = exreportAdd(report60, test60);
report60 = exreportAdd(report60, list(plot601,plot602));

# añade un resumen al report
table60 = tabularExpSummary(exp60,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report60 = exreportAdd(report60, table60);

# report en html
r60 = exreportRender(report60,target='html',visualize=T);


##########################################################################




# TEST 61: tMoPs vs SHENOY

sh = (r$Método=='Shenoy' & r$nInter=='1-7')
tm = (r$Método=='tMoPs' & r$Grado==7 & r$Muestra==1000)
r61 = r[sh | tm,]

# crear el objeto exreport con % de aciertos
exp61 = expCreate(r61,methods='Método',problems='f',
                  name='Method tMoPs vs Shenoy (1 piece, deg=7, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m61 = aggregate(r61$DivKL,list(r61$f),FUN=min);
m611 = rep(m61$x,each=2);
r61Grafica = r61;
r61Grafica$DivKL = r61$DivKL / m611;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp61Grafica = expCreate(r61Grafica,methods='Método',problems='f',
                         name='Method tMoPs vs Shenoy (1 piece, deg=7, sample=1000)');

# contraste sobre el error
test61 = testPaired(exp61,output='DivKL',rankOrder='min');

# graficas
plot611 = plotExpSummary(exp61Grafica,'DivKL');
plot612 = plotExpSummary(exp61Grafica,'DivKL',columns=5);


# titulo del report
report61 = exreport('Report for method tMoPs vs Shenoy (1 piece, deg=7, sample=1000)');

# añadir las diferentes partes del report
report61 = exreportAdd(report61, exp61);
report61 = exreportAdd(report61, test61);
report61 = exreportAdd(report61, list(plot611,plot612));

# añade un resumen al report
table61 = tabularExpSummary(exp61,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report61 = exreportAdd(report61, table61);

# report en html
r61 = exreportRender(report61,target='html',visualize=T);








##########################################################################

# CONTRASTES 2 A 2 CON 3 TROZOS, MUESTRA=1000 Y GRADO=3 (SHENOY, P-B, tMoPs)

##########################################################################



# TEST 62: PEREZ-BERNABE vs SHENOY

sh = (r$Método=='Shenoy' & r$nInter=='3-3')
pb = (r$Método=='Perez-Bernabe' & r$Grado==3 & r$Muestra==1000)
r62 = r[sh | pb,]

# crear el objeto exreport con % de aciertos
exp62 = expCreate(r62,methods='Método',problems='f',
                  name='Method P-B vs Shenoy (3 pieces, deg=3, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m62 = aggregate(r62$DivKL,list(r62$f),FUN=min);
m621 = rep(m62$x,each=2);
r62Grafica = r62;
r62Grafica$DivKL = r62$DivKL / m621;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp62Grafica = expCreate(r62Grafica,methods='Método',problems='f',
                         name='Method P-B vs Shenoy (3 pieces, deg=3, sample=1000)');

# contraste sobre el error
test62 = testPaired(exp62,output='DivKL',rankOrder='min');

# graficas
plot621 = plotExpSummary(exp62Grafica,'DivKL');
plot622 = plotExpSummary(exp62Grafica,'DivKL',columns=5);


# titulo del report
report62 = exreport('Report for method P-B vs Shenoy (3 pieces, deg=3, sample=1000)');

# añadir las diferentes partes del report
report62 = exreportAdd(report62, exp62);
report62 = exreportAdd(report62, test62);
report62 = exreportAdd(report62, list(plot621,plot622));

# añade un resumen al report
table62 = tabularExpSummary(exp62,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report62 = exreportAdd(report62, table62);

# report en html
r62 = exreportRender(report62,target='html',visualize=T);


##########################################################################




# TEST 63: tMoPs vs SHENOY

sh = (r$Método=='Shenoy' & r$nInter=='3-3')
tm = (r$Método=='tMoPs' & r$Grado==3 & r$Muestra==1000)
r63 = r[sh | tm,]

# crear el objeto exreport con % de aciertos
exp63 = expCreate(r63,methods='Método',problems='f',
                  name='Method tMoPs vs Shenoy (3 pieces, deg=3, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m63 = aggregate(r63$DivKL,list(r63$f),FUN=min);
m631 = rep(m63$x,each=2);
r63Grafica = r63;
r63Grafica$DivKL = r63$DivKL / m631;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp63Grafica = expCreate(r63Grafica,methods='Método',problems='f',
                         name='Method tMoPs vs Shenoy (3 pieces, deg=3, sample=1000)');

# contraste sobre el error
test63 = testPaired(exp63,output='DivKL',rankOrder='min');

# graficas
plot631 = plotExpSummary(exp63Grafica,'DivKL');
plot632 = plotExpSummary(exp63Grafica,'DivKL',columns=5);


# titulo del report
report63 = exreport('Report for method tMoPs vs Shenoy (3 pieces, deg=3, sample=1000)');

# añadir las diferentes partes del report
report63 = exreportAdd(report63, exp63);
report63 = exreportAdd(report63, test63);
report63 = exreportAdd(report63, list(plot631,plot632));

# añade un resumen al report
table63 = tabularExpSummary(exp63,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report63 = exreportAdd(report63, table63);

# report en html
r63 = exreportRender(report63,target='html',visualize=T);





##########################################################################

# CONTRASTES 2 A 2 CON 3 TROZOS, MUESTRA=1000 Y GRADO=4 (SHENOY, P-B, tMoPs)

##########################################################################



# TEST 64: PEREZ-BERNABE vs SHENOY

sh = (r$Método=='Shenoy' & r$nInter=='3-4')
pb = (r$Método=='Perez-Bernabe' & r$Grado==4 & r$Muestra==1000)
r64 = r[sh | pb,]

# crear el objeto exreport con % de aciertos
exp64 = expCreate(r64,methods='Método',problems='f',
                  name='Method P-B vs Shenoy (3 pieces, deg=4, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m64 = aggregate(r64$DivKL,list(r64$f),FUN=min);
m641 = rep(m64$x,each=2);
r64Grafica = r64;
r64Grafica$DivKL = r64$DivKL / m641;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp64Grafica = expCreate(r64Grafica,methods='Método',problems='f',
                         name='Method P-B vs Shenoy (3 pieces, deg=4, sample=1000)');

# contraste sobre el error
test64 = testPaired(exp64,output='DivKL',rankOrder='min');

# graficas
plot641 = plotExpSummary(exp64Grafica,'DivKL');
plot642 = plotExpSummary(exp64Grafica,'DivKL',columns=5);


# titulo del report
report64 = exreport('Report for method P-B vs Shenoy (3 pieces, deg=4, sample=1000)');

# añadir las diferentes partes del report
report64 = exreportAdd(report64, exp64);
report64 = exreportAdd(report64, test64);
report64 = exreportAdd(report64, list(plot641,plot642));

# añade un resumen al report
table64 = tabularExpSummary(exp64,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report64 = exreportAdd(report64, table64);

# report en html
r64 = exreportRender(report64,target='html',visualize=T);


##########################################################################




# TEST 65: tMoPs vs SHENOY

sh = (r$Método=='Shenoy' & r$nInter=='3-4')
tm = (r$Método=='tMoPs' & r$Grado==4 & r$Muestra==1000)
r65 = r[sh | tm,]

# crear el objeto exreport con % de aciertos
exp65 = expCreate(r65,methods='Método',problems='f',
                  name='Method tMoPs vs Shenoy (3 pieces, deg=4, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m65 = aggregate(r65$DivKL,list(r65$f),FUN=min);
m651 = rep(m65$x,each=2);
r65Grafica = r65;
r65Grafica$DivKL = r65$DivKL / m651;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp65Grafica = expCreate(r65Grafica,methods='Método',problems='f',
                         name='Method tMoPs vs Shenoy (3 pieces, deg=4, sample=1000)');

# contraste sobre el error
test65 = testPaired(exp65,output='DivKL',rankOrder='min');

# graficas
plot651 = plotExpSummary(exp65Grafica,'DivKL');
plot652 = plotExpSummary(exp65Grafica,'DivKL',columns=5);


# titulo del report
report65 = exreport('Report for method tMoPs vs Shenoy (3 pieces, deg=4, sample=1000)');

# añadir las diferentes partes del report
report65 = exreportAdd(report65, exp65);
report65 = exreportAdd(report65, test65);
report65 = exreportAdd(report65, list(plot651,plot652));

# añade un resumen al report
table65 = tabularExpSummary(exp65,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report65 = exreportAdd(report65, table65);

# report en html
r65 = exreportRender(report65,target='html',visualize=T);






##########################################################################

# CONTRASTES 2 A 2 CON 3 TROZOS, MUESTRA=1000 Y GRADO=7 (SHENOY, P-B, tMoPs)

##########################################################################



# TEST 66: PEREZ-BERNABE vs SHENOY

sh = (r$Método=='Shenoy' & r$nInter=='3-7')
pb = (r$Método=='Perez-Bernabe' & r$Grado==7 & r$Muestra==1000)
r66 = r[sh | pb,]

# crear el objeto exreport con % de aciertos
exp66 = expCreate(r66,methods='Método',problems='f',
                  name='Method P-B vs Shenoy (3 pieces, deg=7, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m66 = aggregate(r66$DivKL,list(r66$f),FUN=min);
m661 = rep(m66$x,each=2);
r66Grafica = r66;
r66Grafica$DivKL = r66$DivKL / m661;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp66Grafica = expCreate(r66Grafica,methods='Método',problems='f',
                         name='Method P-B vs Shenoy (3 pieces, deg=7, sample=1000)');

# contraste sobre el error
test66 = testPaired(exp66,output='DivKL',rankOrder='min');

# graficas
plot661 = plotExpSummary(exp66Grafica,'DivKL');
plot662 = plotExpSummary(exp66Grafica,'DivKL',columns=5);


# titulo del report
report66 = exreport('Report for method P-B vs Shenoy (3 pieces, deg=7, sample=1000)');

# añadir las diferentes partes del report
report66 = exreportAdd(report66, exp66);
report66 = exreportAdd(report66, test66);
report66 = exreportAdd(report66, list(plot661,plot662));

# añade un resumen al report
table66 = tabularExpSummary(exp66,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report66 = exreportAdd(report66, table66);

# report en html
r66 = exreportRender(report66,target='html',visualize=T);


##########################################################################




# TEST 67: tMoPs vs SHENOY

sh = (r$Método=='Shenoy' & r$nInter=='3-7')
tm = (r$Método=='tMoPs' & r$Grado==7 & r$Muestra==1000)
r67 = r[sh | tm,]

# crear el objeto exreport con % de aciertos
exp67 = expCreate(r67,methods='Método',problems='f',
                  name='Method tMoPs vs Shenoy (3 pieces, deg=7, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m67 = aggregate(r67$DivKL,list(r67$f),FUN=min);
m671 = rep(m67$x,each=2);
r67Grafica = r67;
r67Grafica$DivKL = r67$DivKL / m671;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp67Grafica = expCreate(r67Grafica,methods='Método',problems='f',
                         name='Method tMoPs vs Shenoy (3 pieces, deg=7, sample=1000)');

# contraste sobre el error
test67 = testPaired(exp67,output='DivKL',rankOrder='min');

# graficas
plot671 = plotExpSummary(exp67Grafica,'DivKL');
plot672 = plotExpSummary(exp67Grafica,'DivKL',columns=5);


# titulo del report
report67 = exreport('Report for method tMoPs vs Shenoy (3 pieces, deg=7, sample=1000)');

# añadir las diferentes partes del report
report67 = exreportAdd(report67, exp67);
report67 = exreportAdd(report67, test67);
report67 = exreportAdd(report67, list(plot671,plot672));

# añade un resumen al report
table67 = tabularExpSummary(exp67,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report67 = exreportAdd(report67, table67);

# report en html
r67 = exreportRender(report67,target='html',visualize=T);





##########################################################################

# CONTRASTES 2 A 2 CON 2 TROZOS, MUESTRA=1000 (SHENOY, tMoPs)

##########################################################################



# TEST 68: tMoPs vs SHENOY (GRADO=3)

sh = (r$Método=='Shenoy' & r$nInter=='2-3')
tm = (r$Método=='tMoPs' & r$Grado==3 & r$Muestra==1000)
r68 = r[sh | tm,]

# crear el objeto exreport con % de aciertos
exp68 = expCreate(r68,methods='Método',problems='f',
                  name='Method tMoPs vs Shenoy (2 pieces, deg=3, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m68 = aggregate(r68$DivKL,list(r68$f),FUN=min);
m681 = rep(m68$x,each=2);
r68Grafica = r68;
r68Grafica$DivKL = r68$DivKL / m681;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp68Grafica = expCreate(r68Grafica,methods='Método',problems='f',
                         name='Method tMoPs vs Shenoy (2 pieces, deg=3, sample=1000)');

# contraste sobre el error
test68 = testPaired(exp68,output='DivKL',rankOrder='min');

# graficas
plot681 = plotExpSummary(exp68Grafica,'DivKL');
plot682 = plotExpSummary(exp68Grafica,'DivKL',columns=5);


# titulo del report
report68 = exreport('Report for method tMoPs vs Shenoy (2 pieces, deg=3, sample=1000)');

# añadir las diferentes partes del report
report68 = exreportAdd(report68, exp68);
report68 = exreportAdd(report68, test68);
report68 = exreportAdd(report68, list(plot681,plot682));

# añade un resumen al report
table68 = tabularExpSummary(exp68,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report68 = exreportAdd(report68, table68);

# report en html
r68 = exreportRender(report68,target='html',visualize=T);


##########################################################################



# TEST 69: tMoPs vs SHENOY (GRADO=4)

sh = (r$Método=='Shenoy' & r$nInter=='2-4')
tm = (r$Método=='tMoPs' & r$Grado==4 & r$Muestra==1000)
r69 = r[sh | tm,]

# crear el objeto exreport con % de aciertos
exp69 = expCreate(r69,methods='Método',problems='f',
                  name='Method tMoPs vs Shenoy (2 pieces, deg=4, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m69 = aggregate(r69$DivKL,list(r69$f),FUN=min);
m691 = rep(m69$x,each=2);
r69Grafica = r69;
r69Grafica$DivKL = r69$DivKL / m691;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp69Grafica = expCreate(r69Grafica,methods='Método',problems='f',
                         name='Method tMoPs vs Shenoy (2 pieces, deg=4, sample=1000)');

# contraste sobre el error
test69 = testPaired(exp69,output='DivKL',rankOrder='min');

# graficas
plot691 = plotExpSummary(exp69Grafica,'DivKL');
plot692 = plotExpSummary(exp69Grafica,'DivKL',columns=5);


# titulo del report
report69 = exreport('Report for method tMoPs vs Shenoy (2 pieces, deg=4, sample=1000)');

# añadir las diferentes partes del report
report69 = exreportAdd(report69, exp69);
report69 = exreportAdd(report69, test69);
report69 = exreportAdd(report69, list(plot691,plot692));

# añade un resumen al report
table69 = tabularExpSummary(exp69,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report69 = exreportAdd(report69, table69);

# report en html
r69 = exreportRender(report69,target='html',visualize=T);


##########################################################################




# TEST 70: tMoPs vs SHENOY (GRADO=7)

sh = (r$Método=='Shenoy' & r$nInter=='2-7')
tm = (r$Método=='tMoPs' & r$Grado==7 & r$Muestra==1000)
r70 = r[sh | tm,]

# crear el objeto exreport con % de aciertos
exp70 = expCreate(r70,methods='Método',problems='f',
                  name='Method tMoPs vs Shenoy (2 pieces, deg=7, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m70 = aggregate(r70$DivKL,list(r70$f),FUN=min);
m701 = rep(m70$x,each=2);
r70Grafica = r70;
r70Grafica$DivKL = r70$DivKL / m701;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp70Grafica = expCreate(r70Grafica,methods='Método',problems='f',
                         name='Method tMoPs vs Shenoy (2 pieces, deg=7, sample=1000)');

# contraste sobre el error
test70 = testPaired(exp70,output='DivKL',rankOrder='min');

# graficas
plot701 = plotExpSummary(exp70Grafica,'DivKL');
plot702 = plotExpSummary(exp70Grafica,'DivKL',columns=5);


# titulo del report
report70 = exreport('Report for method tMoPs vs Shenoy (2 pieces, deg=7, sample=1000)');

# añadir las diferentes partes del report
report70 = exreportAdd(report70, exp70);
report70 = exreportAdd(report70, test70);
report70 = exreportAdd(report70, list(plot701,plot702));

# añade un resumen al report
table70 = tabularExpSummary(exp70,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report70 = exreportAdd(report70, table70);

# report en html
r70 = exreportRender(report70,target='html',visualize=T);







##########################################################################

# CONTRASTES 2 A 2 CON 4 TROZOS, MUESTRA=1000 (SHENOY, tMoPs)

##########################################################################



# TEST 71: tMoPs vs SHENOY (GRADO=3)

sh = (r$Método=='Shenoy' & r$nInter=='4-3')
tm = (r$Método=='tMoPs' & r$Grado==3 & r$Muestra==1000)
r71 = r[sh | tm,]

# crear el objeto exreport con % de aciertos
exp71 = expCreate(r71,methods='Método',problems='f',
                  name='Method tMoPs vs Shenoy (4 pieces, deg=3, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m71 = aggregate(r71$DivKL,list(r71$f),FUN=min);
m711 = rep(m71$x,each=2);
r71Grafica = r71;
r71Grafica$DivKL = r71$DivKL / m711;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp71Grafica = expCreate(r71Grafica,methods='Método',problems='f',
                         name='Method tMoPs vs Shenoy (4 pieces, deg=3, sample=1000)');

# contraste sobre el error
test71 = testPaired(exp71,output='DivKL',rankOrder='min');

# graficas
plot711 = plotExpSummary(exp71Grafica,'DivKL');
plot712 = plotExpSummary(exp71Grafica,'DivKL',columns=5);


# titulo del report
report71 = exreport('Report for method tMoPs vs Shenoy (4 pieces, deg=3, sample=1000)');

# añadir las diferentes partes del report
report71 = exreportAdd(report71, exp71);
report71 = exreportAdd(report71, test71);
report71 = exreportAdd(report71, list(plot711,plot712));

# añade un resumen al report
table71 = tabularExpSummary(exp71,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report71 = exreportAdd(report71, table71);

# report en html
r71 = exreportRender(report71,target='html',visualize=T);


##########################################################################




# TEST 72: tMoPs vs SHENOY (GRADO=4)

sh = (r$Método=='Shenoy' & r$nInter=='4-4')
tm = (r$Método=='tMoPs' & r$Grado==4 & r$Muestra==1000)
r72 = r[sh | tm,]

# crear el objeto exreport con % de aciertos
exp72 = expCreate(r72,methods='Método',problems='f',
                  name='Method tMoPs vs Shenoy (4 pieces, deg=4, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m72 = aggregate(r72$DivKL,list(r72$f),FUN=min);
m721 = rep(m72$x,each=2);
r72Grafica = r72;
r72Grafica$DivKL = r72$DivKL / m721;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp72Grafica = expCreate(r72Grafica,methods='Método',problems='f',
                         name='Method tMoPs vs Shenoy (4 pieces, deg=4, sample=1000)');

# contraste sobre el error
test72 = testPaired(exp72,output='DivKL',rankOrder='min');

# graficas
plot721 = plotExpSummary(exp72Grafica,'DivKL');
plot722 = plotExpSummary(exp72Grafica,'DivKL',columns=5);


# titulo del report
report72 = exreport('Report for method tMoPs vs Shenoy (4 pieces, deg=4, sample=1000)');

# añadir las diferentes partes del report
report72 = exreportAdd(report72, exp72);
report72 = exreportAdd(report72, test72);
report72 = exreportAdd(report72, list(plot721,plot722));

# añade un resumen al report
table72 = tabularExpSummary(exp72,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report72 = exreportAdd(report72, table72);

# report en html
r72 = exreportRender(report72,target='html',visualize=T);


##########################################################################




# TEST 73: tMoPs vs SHENOY (GRADO=7)

sh = (r$Método=='Shenoy' & r$nInter=='4-7')
tm = (r$Método=='tMoPs' & r$Grado==7 & r$Muestra==1000)
r73 = r[sh | tm,]

# crear el objeto exreport con % de aciertos
exp73 = expCreate(r73,methods='Método',problems='f',
                  name='Method tMoPs vs Shenoy (4 pieces, deg=7, sample=1000)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m73 = aggregate(r73$DivKL,list(r73$f),FUN=min);
m731 = rep(m73$x,each=2);
r73Grafica = r73;
r73Grafica$DivKL = r73$DivKL / m731;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp73Grafica = expCreate(r73Grafica,methods='Método',problems='f',
                         name='Method tMoPs vs Shenoy (4 pieces, deg=7, sample=1000)');

# contraste sobre el error
test73 = testPaired(exp73,output='DivKL',rankOrder='min');

# graficas
plot731 = plotExpSummary(exp73Grafica,'DivKL');
plot732 = plotExpSummary(exp73Grafica,'DivKL',columns=5);


# titulo del report
report73 = exreport('Report for method tMoPs vs Shenoy (4 pieces, deg=7, sample=1000)');

# añadir las diferentes partes del report
report73 = exreportAdd(report73, exp73);
report73 = exreportAdd(report73, test73);
report73 = exreportAdd(report73, list(plot731,plot732));

# añade un resumen al report
table73 = tabularExpSummary(exp73,'DivKL',digits=4,format="f",boldfaceColumns='min',
                            tableSplit=2);
report73 = exreportAdd(report73, table73);

# report en html
r73 = exreportRender(report73,target='html',visualize=T);






##########################################################################

# CONTRASTES 2 A 2 CON DE DATOS REALES SEGÚN EL GRADO (P-B, tMoPs)

##########################################################################


# TEST 80: PEREZ-BERNABE vs tMoPs (degree=3)

rr80 = rr[rr$Grado==3,]

# crear el objeto exreport con % de aciertos
exp80 = expCreate(rr80,methods='Método',problems='Variable',
                  name='Method P-B vs tMoPs (deg=3)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m80 = aggregate(rr80$BIC,list(rr80$Variable),FUN=max);
m801 = rep(m80$x,each=2);
r80Grafica = rr80;
r80Grafica$DivKL = rr80$BIC / m801;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp80Grafica = expCreate(r80Grafica,methods='Método',problems='Variable',
                         name='Method P-B vs tMoPs (deg=3)');

# contraste sobre el error
test80 = testPaired(exp80,output='BIC',rankOrder='max');

# graficas
plot801 = plotExpSummary(exp80Grafica,'BIC');
plot802 = plotExpSummary(exp80Grafica,'BIC',columns=5);


# titulo del report
report80 = exreport('Report for method P-B vs tMoPs (deg=3)');

# añadir las diferentes partes del report
report80 = exreportAdd(report80, exp80);
report80 = exreportAdd(report80, test80);
report80 = exreportAdd(report80, list(plot801,plot802));

# añade un resumen al report
table80 = tabularExpSummary(exp80,'BIC',digits=4,format="f",boldfaceColumns='max',
                            tableSplit=2);
report80 = exreportAdd(report80, table80);

# report en html
rr80 = exreportRender(report80,target='html',visualize=T);


##########################################################################


# TEST 81: PEREZ-BERNABE vs tMoPs (degree=4)

rr81 = rr[rr$Grado==4,]

# crear el objeto exreport con % de aciertos
exp81 = expCreate(rr81,methods='Método',problems='Variable',
                  name='Method P-B vs tMoPs (deg=4)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m81 = aggregate(rr81$BIC,list(rr81$Variable),FUN=max);
m811 = rep(m81$x,each=2);
r81Grafica = rr81;
r81Grafica$DivKL = rr81$BIC / m811;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp81Grafica = expCreate(r81Grafica,methods='Método',problems='Variable',
                         name='Method P-B vs tMoPs (deg=4)');

# contraste sobre el error
test81 = testPaired(exp81,output='BIC',rankOrder='max');

# graficas
plot811 = plotExpSummary(exp81Grafica,'BIC');
plot812 = plotExpSummary(exp81Grafica,'BIC',columns=5);


# titulo del report
report81 = exreport('Report for method P-B vs tMoPs (deg=4)');

# añadir las diferentes partes del report
report81 = exreportAdd(report81, exp81);
report81 = exreportAdd(report81, test81);
report81 = exreportAdd(report81, list(plot811,plot812));

# añade un resumen al report
table81 = tabularExpSummary(exp81,'BIC',digits=4,format="f",boldfaceColumns='max',
                            tableSplit=2);
report81 = exreportAdd(report81, table81);

# report en html
rr81 = exreportRender(report81,target='html',visualize=T);


##########################################################################



# TEST 82: PEREZ-BERNABE vs tMoPs (degree=7)

rr82 = rr[rr$Grado==7,]

# crear el objeto exreport con % de aciertos
exp82 = expCreate(rr82,methods='Método',problems='Variable',
                  name='Method P-B vs tMoPs (deg=7)');

# calcular el minimo % de cada base de datos y calcular el resultado porcentual
m82 = aggregate(rr82$BIC,list(rr82$Variable),FUN=max);
m821 = rep(m82$x,each=2);
r82Grafica = rr82;
r82Grafica$DivKL = rr82$BIC / m821;

# crear el objeto exreport con % de aciertos relativos para gr?fica
exp82Grafica = expCreate(r82Grafica,methods='Método',problems='Variable',
                         name='Method P-B vs tMoPs (deg=7)');

# contraste sobre el error
test82 = testPaired(exp82,output='BIC',rankOrder='max');

# graficas
plot821 = plotExpSummary(exp82Grafica,'BIC');
plot822 = plotExpSummary(exp82Grafica,'BIC',columns=5);


# titulo del report
report82 = exreport('Report for method P-B vs tMoPs (deg=7)');

# añadir las diferentes partes del report
report82 = exreportAdd(report82, exp82);
report82 = exreportAdd(report82, test82);
report82 = exreportAdd(report82, list(plot821,plot822));

# añade un resumen al report
table82 = tabularExpSummary(exp82,'BIC',digits=4,format="f",boldfaceColumns='max',
                            tableSplit=2);
report82 = exreportAdd(report82, table82);

# report en html
rr82 = exreportRender(report82,target='html',visualize=T);




