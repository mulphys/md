Xa0=[] # initial coordinates of particle a
Xb0=[] # initial coordinates of particle b
ta # inittal time of a
tb # initial time of b
Xa[i] = Xa0[i] + Va[i]*(t-ta)
Xb[i] = Xb0[i] + Vb[i]*(t-tb)
dX[i]=Xb[i]-Xa[i]
d = dX[i]**2
= Xa0[i]*Xb0[i] + (t-ta)*Va[i]*Xb0[i] + (t-tb)*Vb[i]*Xa0[i] + (t-ta)*(t-tb)*Va[i]*Vb[i]
= 
d < Ra + Rb
dX0 = 0
for i in range(DIM)
	dX0+=dX[i]**2
dX0=sqrt(dX0)
