% http://www.dayah.com/periodic/
% http://en.wikipedia.org/wiki/Molecular_mass
% http://en.wikipedia.org/wiki/Hydrogen_sulfide
Mu=1.66053886e-27 % kg or 6.022e-23g (?)   molecular unit mass
MH=1.0079     % u: Mass of hydrogen atom
MH2=2*MH
mH2=Mu*MH2
MS=32.065     % u: Mass of sulfur
MH2S=2*MH+MS %=: mass of hydrogen sulfide
mH2S=Mu*MH2S %=: mass of hydrogen sulfide
MNi=58.693
mNi=MNi*Mu %=: mass of nickel
mH2=3.34731424e-27  % kg: mass of H2 molecule
mH2S=5.659428616e-26 % kg: mass of H2S molecule

k=1.38066e-23 %J/K: =8.617385e-5V/K  Boltzman constant
T=300 % K: temperature
% molecular velocity:

m=mH2
v=(8.0*k*T/(pi*m))^.5

m=mH2S
v=(8.0*k*T/(pi*m))^.5

%Number of molecules: PV = NkT: N=PV/kT
P=1e5 % Pa
Dnm=10 % nm: domain size in nano-meters
D=Dnm*1e-9 % m: domain size in meters
V=(D)^3 % m^3: domain volume 

N=P*V/(k*T) %=2.4e7: 20 million in 1000nm cube

% molecule sizes:
d=3.0e-10 %m: average molecule size 

% Avarage separation:
l=(V/N)^(1/3) %=3.5e-9=35.0e-10

l/d %=11.533

% COLLISION: 
% MOLECULS BOUNCE:
% In the center of mass frame:
% Momentum:
% mA2*vA2 + mB2*vB2 = 0  
% Energy:
% mA1*vA1^2 + mB1*vB1^2 = 2*K1 = mA2*vA2^2 + mB2*vB2^2
% from momentum: 
% vB2=-vA2*mA2/mB2
% from energy:
% mA2*vA2^2 + mB2*(-vA2*mA2/mB2)^2 = 2*K1
% vA2 = (2*K1/(mA2 + mA2^2/mB2))^0.5 = (2*K1/(mA2*(1+mA2/mB2)))^0.5
% STICK TOGETHER: 
% mA1+mA2=m2
% Momentum:
% m2*v2 = 0
% Energy:
% K1 = 0
% Internal energy of the molecule is increased by: dE=K1
%
% http://en.wikipedia.org/wiki/Kinetic_theory
% the number of atomic or molecular collisions with a wall of a container per unit area per unit time:
% n=0.25*N/V*Vav = 0.25*rho/m*Vav  (1)
% Vav=(8*k*T/(pi*m))^0.5 = average molecular velocity
% Also:
% http://www.seas.upenn.edu/~waddavi/old/education/courses/chem101/kineticmoleculartheory.pdf
% PV = 1/3*N*m*Vav^2
% P= 1/3*rho*Vav^2
% We can rewrite (1) in terms of pressure as:
% n=0.25*3*P*Vav^(-2)*Vav/m = 0.75*P/(Vav*m)

% Distribution on the energy levels:
X=-5.0:0.1:5.0
n=length(X)
Y=zeros(1,n);
for i=1:n
	Y(i)=1.0-exp(-abs(X(i)));
	if X(i)<0.0
		Y(i)=-Y(i);
	end
end
plot(X,Y)
Y=0.0:0.05:1.0
n=length(Y)
X=zeros(1,n);
for i=1:n
	if Y(i)<=0.5
		y=2*Y(i);
		X(i)=log(y);
	else
		y=1.0-2.0*(Y(i)-0.5);
		X(i)=-log(y);
	end
end
plot(Y,X)
