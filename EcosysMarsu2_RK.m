clear all ,close all

#Initialisation :
tic
#  Variables communes :

t0 = 1;
tf = 365*12;
dt = 0.1;
n = tf/dt;
t = t0;
T = t; #Matrice dans laquelle on stocke tous les temps


# Modèle NPZ


N0 = 4; #μ.mole.N.l.^-1 
P0 = 2.5; #μ.mole.N.l.^-1
Z0 = 0.5; #μ.mole.N.l.^-1

N = N0; #Abondance d'azote (en μ.mole.N.l.^-1) et à t = 0
P = P0; #Abondance d'algue (en μ.mole.N.l.^-1) et à t = 0
Z = Z0; #Abondance de Métynnis (en μ.mole.N.l.^-1) et à t = 0

Nsol = N; #Matrice dans laquelle on stocke toutes les solutions de N
Psol = P; #Matrice dans laquelle on stocke toutes les solutions de P
Zsol = Z; #Matrice dans laquelle on stocke toutes les solutions de Z

IntensiteLum = load('IntensiteLum.txt');

# Modèle croissance Marsu

L0 = 7; #cm 

L = L0; #Taille du Marsupilami (en cm) et à T0

Lsol = L; #Matrice dans laquelle on stocke toutes les solutions de L

Data = load('Temperatures.txt'); #On stocke la température dans une matrice


for i = 2:n #On commence la boucle à 2 car la valeur du jour 1 est celui à t = 0
  
  
  [NRK1,PRK1,ZRK1] = NPZ(N,P,Z,IntensiteLum,t);
  [NRK2,PRK2,ZRK2] = NPZ(N + 1/2*dt*NRK1,P + 1/2*dt*PRK1,Z + 1/2*dt*ZRK1,IntensiteLum,t);
  [NRK3,PRK3,ZRK3] = NPZ(N + 1/2*dt*NRK2,P + 1/2*dt*PRK2,Z + 1/2*dt*ZRK2,IntensiteLum,t);
  [NRK4,PRK4,ZRK4] = NPZ(N + dt*NRK3,P + dt*PRK3,Z + dt*ZRK3,IntensiteLum,t);
  
  
  dNdt = 1/6*NRK1 + 1/3*NRK2 + 1/3*NRK3 + 1/6*NRK4; 
  dPdt = 1/6*PRK1 + 1/3*PRK2 + 1/3*PRK3 + 1/6*PRK4;
  dZdt = 1/6*ZRK1 + 1/3*ZRK2 + 1/3*ZRK3 + 1/6*ZRK4;
  
  # Résolution des équations du modèle NPZ
  
  
  N = N + dt*dNdt;    #On calcule les solutions de N,P,Z pour chaque valeur de dt
  P = P + dt*dPdt;
  Z = Z + dt*dZdt;
  
  
  
  # Résolution de l'équation différentielle du modèle de croissance du Marsu
  
  [LRK1] = CroissMarsu(t,Z,L,Data);
  [LRK2] = CroissMarsu(t,Z,L + 1/2*dt*LRK1,Data);
  [LRK3] = CroissMarsu(t,Z,L + 1/2*dt*LRK2,Data);
  [LRK4] = CroissMarsu(t,Z,L + 1/2*dt*LRK3,Data);
  
  dLdt = 1/6*LRK1 + 1/3*LRK2 + 1/3*LRK3 + 1/6*LRK4;
  
  if dLdt < 0 #Si le flux de croissance du Marsu est inférieur à 0, il devient nul. En effet, notre Marsu ne peut pas raptisser !
     dLdt = 0;  
     
  end
  
  L = L + dLdt*dt;
  
  
  Nsol = [Nsol ; N]; #On stocke toutes les solutions de N,P,Z,L dans des matrices
  Psol = [Psol ; P];
  Zsol = [Zsol ; Z];
  Lsol = [Lsol ; L];
  
   t = t + dt;
  T = [T ; t];
end




# Graphiques



figure(1) 

#subplot(2,1,1)

title('Modèle de production lacustre & de croisssance du Marsu')
plot(T,Nsol,'b',T,Psol,'g',T,Zsol,'r') # Toutes les courbes sur un même graphique
legend('Azote','Algue','Metynnis')
xlabel('Temps (jour)')
ylabel('Concentration en Azote (\mumole N.l.^{-1})')
grid on


figure (2)

#subplot(2,1,2)

[hax,h1,h2] = plotyy(T,Lsol,T,Zsol);
set(hax(1),'ycolor','k');
set(hax(2),'ycolor','r');
set (h1, 'color','k','linestyle','-');
set (h2, 'color','r');
#legend('Marsu','Metynnis');
xlabel('Temps (jour)');
ylabel(hax(1),'Taillen (cm)');
ylabel(hax(2),'Concentration en Azote (\mumole N.l.^{-1})');
grid on

toc