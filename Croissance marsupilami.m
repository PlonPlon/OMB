Lmax=150;% Linf=130 ou 3597
KTref=0.001029;% K celui de la station 101 car berceau de l'espèce
Tsim=365*1;% durée de la simulation
L=7;% longueur initiale du marsupilami
dt=1;
data=[];
TA=9048;% c'est un paramètre, il ne change pas
Tref=310;% température prise pour 37°C car température du berceau de l'espèce
ampT=10.366;
Xk=0.14; %micromol d'azote par litre
poisson=Z;
met=poisson(:,1);
Temp=load('Temp2011.txt');
T=[Temp];
Tmoy=mean(T(:,2))
x=1

                                                                                                   
for t=1:dt:Tsim
	Tt=273+Tmoy+0.5*ampT*cos((2*pi*t)/365);
	N=met(x,:);
	f=N/(Xk+N);
	cT=exp((TA/Tref)-(TA/Tt));
	dLdt=cT*KTref*(f*Lmax-L);
	L=L+dLdt*dt;
	data=[data;t L];
	x=x+1;
end
