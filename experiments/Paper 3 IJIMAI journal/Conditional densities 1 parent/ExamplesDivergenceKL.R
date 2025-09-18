


f = function(x,y)   2*x + y + 1
g = function(x,y)   y + 3
p = c(3,1)

# function of the first integral
fToInt = function(x,y,f1=f,g1=g,coef=p)    f1(x,y) +g1(y,coef)


g = function(y,coef=p)
{
  s = 0;
  for (i in 1:length(coef))
  {
    s = s + coef[i]*y^(i-1);
  }
  return(s)
}

integral2(fToInt,1,3,0,2)





############################################################

fY = function(y)   2*y
fXY = function(x,y)   3*x^2 + 2*y + 4
p = c(1,2,1)

# function of the first integral
fToInt = function(x,y,f1Y=fY,f1XY=fXY,g1=g,coef=p)    
  f1Y(y) * f1XY(x,y) * log(f1XY(x,y)/g1(x,coef))


g = function(y,coef=p)
{
  s = 0;
  for (i in 1:length(coef))
  {
    s = s + coef[i]*y^(i-1);
  }
  return(s)
}


# function of each integral
fToInt1 = function(x,y)  (6*x^2*y + 4*y^2+8*y) * log(3*x^2+2*y+4)
fToInt2 = function(x,y)  (6*x^2*y + 4*y^2+8*y) * log(x^2+2*x+1)

integral2(fToInt1,0,2,0,2)
integral2(fToInt2,0,2,0,2)

DivKLpolyXY(fY,fXY,p,0,2,0,2)





############################################################

fY = function(y)   2*y
fXY = function(x,y)   3*x^2 + 2*y + 4

v1 = cbind(c(0,3),c(0,4),c(0,3));            p1 = c(1,2,1)
v2 = cbind(c(1,3),c(0,4),c(2,5));            p2 = c(3,1,NA)

v = cbind(v1,v2);            p = rbind(p1,p2)

tMoP = list(v,p,c(0,2,4))


# functions of each integral
fToInt1 = function(x,y)  (6*x^2*y + 4*y^2+8*y) * log(3*x^2+2*y+4)

fToInt21 = function(x,y)  (6*x^2*y + 4*y^2+8*y) * log(x^2+2*x+1)
fToInt22 = function(x,y)  (6*x^2*y + 4*y^2+8*y) * log(3)

fToInt23 = function(x,y)  (6*x^2*y + 4*y^2+8*y) * log(2)
fToInt24 = function(x,y)  (6*x^2*y + 4*y^2+8*y) * log(x+3)
fToInt25 = function(x,y)  (6*x^2*y + 4*y^2+8*y) * log(5)

# step by step
i1 = integral2(fToInt1,0,3,0,2)
i2 = integral2(fToInt21,0,3,0,2)
i3 = integral2(fToInt1,3,4,0,2)
i4 = integral2(fToInt22,3,4,0,2)

i5 = integral2(fToInt1,0,1,2,4)
i6 = integral2(fToInt23,0,1,2,4)
i7 = integral2(fToInt1,1,3,2,4)
i8 = integral2(fToInt24,1,3,2,4)
i9 = integral2(fToInt1,3,4,2,4)
i10 = integral2(fToInt25,3,4,2,4)

i1$Q - i2$Q + i3$Q - i4$Q + i5$Q - i6$Q + i7$Q - i8$Q + i9$Q - i10$Q



# using new function
DivKL1parent(fY,fXY,tMoP,0,4,0,4)





############################################################

fY = function(y)   2*y
fZ = function(z)   z^2+z+3
fXYZ = function(x,y,z)   3*x^2*z + z + 2*y + 4

v1 = cbind(c(0,3),c(0,4),c(0,3));            p1 = c(1,2,1)
v2 = cbind(c(0,4),c(0,4),c(0,0));            p2 = c(1,2,NA)
v3 = cbind(c(1,3),c(0,4),c(2,5));            p3 = c(3,1,NA)
v4 = cbind(c(0,4),c(0,4),c(0,0));            p4 = c(1,-1,1)
v5 = cbind(c(0,2),c(0,4),c(0,1));            p5 = c(0,1,NA)
v6 = cbind(c(0,4),c(0,4),c(0,0));            p6 = c(6,-1,NA)
v7 = cbind(c(1,4),c(0,4),c(2,0));            p7 = c(5,2,NA)

v = cbind(v1,v2,v3,v4,v5,v6,v7);            p = rbind(p1,p2,p3,p4,p5,p6,p7);

d = c(rep(c(1,1,1,1,2,2,3,3),4),rep(c(4,4,5,5,7,7,7,7),2),rep(c(6,6,6,6,7,7,7,7),2))
discretYZ = t(matrix(d,ncol=8))

tMoP = list(v,p,discretYZ)


# functions of each integral
fToInt1 = function(x,y,z)  2*y * (z^2+z+3) * (3*x^2*z+z+2*y+4) * log(3*x^2*z+z+2*y+4)

fToInt21 = function(x,y,z)  2*y * (z^2+z+3) * (3*x^2*z+z+2*y+4) * log(x^2+2*x+1)
fToInt22 = function(x,y,z)  2*y * (z^2+z+3) * (3*x^2*z+z+2*y+4) * log(3)

fToInt23 = function(x,y,z)  2*y * (z^2+z+3) * (3*x^2*z+z+2*y+4) * log(2*x+1)

fToInt24 = function(x,y,z)  2*y * (z^2+z+3) * (3*x^2*z+z+2*y+4) * log(2)
fToInt25 = function(x,y,z)  2*y * (z^2+z+3) * (3*x^2*z+z+2*y+4) * log(x+3)
fToInt26 = function(x,y,z)  2*y * (z^2+z+3) * (3*x^2*z+z+2*y+4) * log(5)

fToInt27 = function(x,y,z)  2*y * (z^2+z+3) * (3*x^2*z+z+2*y+4) * log(x^2-x+1)

fToInt28 = function(x,y,z)  2*y * (z^2+z+3) * (3*x^2*z+z+2*y+4) * log(x)
fToInt29 = function(x,y,z)  2*y * (z^2+z+3) * (3*x^2*z+z+2*y+4) * log(1)

fToInt210 = function(x,y,z)  2*y * (z^2+z+3) * (3*x^2*z+z+2*y+4) * log(-x+6)

fToInt211 = function(x,y,z)  2*y * (z^2+z+3) * (3*x^2*z+z+2*y+4) * log(2)
fToInt212 = function(x,y,z)  2*y * (z^2+z+3) * (3*x^2*z+z+2*y+4) * log(2*x+5)


# step by step
i1 = integral3(fToInt1,0,3,0,4,0,8);            i2 = integral3(fToInt21,0,3,0,4,0,8);
i3 = integral3(fToInt1,3,4,0,4,0,8);            i4 = integral3(fToInt22,3,4,0,4,0,8);

i5 = integral3(fToInt1,0,4,0,4,8,12);           i6 = integral3(fToInt23,0,4,0,4,8,12);

i7 = integral3(fToInt1,0,1,0,4,12,16);          i8 = integral3(fToInt24,0,1,0,4,12,16);
i9 = integral3(fToInt1,1,3,0,4,12,16);          i10 = integral3(fToInt25,1,3,0,4,12,16);
i11 = integral3(fToInt1,3,4,0,4,12,16);         i12 = integral3(fToInt26,3,4,0,4,12,16);

i13 = integral3(fToInt1,0,4,4,6,0,4);           i14 = integral3(fToInt27,0,4,4,6,0,4);

i15 = integral3(fToInt1,0,2,4,6,4,8);           i16 = integral3(fToInt28,0,2,4,6,4,8);
i17 = integral3(fToInt1,2,4,4,6,4,8);           i18 = integral3(fToInt29,2,4,4,6,4,8);

i19 = integral3(fToInt1,0,4,6,8,0,8);           i20 = integral3(fToInt210,0,4,6,8,0,8);

i21 = integral3(fToInt1,0,1,4,8,8,16);          i22 = integral3(fToInt211,0,1,4,8,8,16);
i23 = integral3(fToInt1,1,4,4,8,8,16);          i24 = integral3(fToInt212,1,4,4,8,8,16);


KL1 = i1-i2+i3-i4+i5-i6+i7-i8+i9-i10+i11-i12+i13-i14+i15-i16+i17-i18+i19-i20+i21-i22+i23-i24



# using new function
KL2 = DivKL2parent(fY,fZ,fXYZ,tMoP,0,4,0,8,0,16)
KL1==KL2


