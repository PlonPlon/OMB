clear all

Linf=120; # Longueur théorique max en cm
L=7; # Longueur initiale en cm
K=0.001; # Constante en jour^-1
T=[0:3650]; # Vecteur temps
Long=[]; # Tableau taille du marsu
Data=[];
dt=1;
t=0;

for i=0:3650
  dLdt=K*(Linf - L);
  L=L+dLdt*dt;
  t=t+dt;
  Data=[Data t];
  Long=[Long L];
end

Lex=Linf*(1-exp(-K*Data));

figure(1)
plot(T,Long,'-k',T,Lex,'r')
ylabel("Longueur (cm)")
xlabel("Temps (jour)")
legend('dLdt=K(Linf-L)','L(t)=Linf*(1-exp(-K*(t-t0)))')
title("Evolution du taux de croissance du Marsupilami à la Station 301 en fonction du temps")
grid

# Savoir comparer la solution exacte et la solution approximative