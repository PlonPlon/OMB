clear all, close all

#Initialisation :

#  Variables communes :

t0 = 1;
tf = 365*4;
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


# Modèle croissance Marsu

L0 = 7; #cm 

L = L0; #Taille du Marsupilami (en cm) et à T0

Lsol = L; #Matrice dans laquelle on stocke toutes les solutions de L


for i = 2:n #On commence la boucle à 2 car la valeur du jour 1 est celui à t = 0
  
  t = t + dt;
  
  [dNdt,dPdt,dZdt,dLdt] = Equadif(N,P,Z,L)
  
  # Résolution des équations du modèle NPZ
  
  
  N = N + dNdt*dt; #On calcul les solutions de N,P,Z pour chaque valeur de dt
  P = P + dPdt*dt;
  Z = Z + dZdt*dt;
  
  
  
  
  # Résolution de l'équation différentielle du modèle de croissance du Marsu
  
  
  if dLdt < 0 #Si le flux de croissance du Marsu est inférieur à 0, il devient nul. En effet, notre Marsu ne peut pas raptisser !
     dLdt = 0;  
     
  end
  
  L = L + dLdt*dt;
  
  
  Nsol = [Nsol ; N]; #On stocke toutes les solutions de N,P,Z,L dans des matrices
  Psol = [Psol ; P];
  Zsol = [Zsol ; Z];
  Lsol = [Lsol ; L];
  
  T = [T ; t];
end




# Graphiques



figure(1) 

#subplot(2,1,1)

plot(T,Nsol,'b',T,Psol,'g',T,Zsol,'r') # Toutes les courbes sur un même graphique
legend('Azote','Algue','Metynnis')
xlabel('Temps (jour)')
ylabel('Concentration en Azote (\mumole N.l.^{-1})')
grid on


#subplot(2,1,2)
figure (2)
[hax,h1,h2] = plotyy(T,Lsol,T,Zsol);
set(hax(1),'ycolor','k');
set(hax(2),'ycolor','r');
set (h1, 'color','k','linestyle','-');
set (h2, 'color','r');
legend('Marsu','Metynnis');
xlabel('Temps (jour)');
ylabel(hax(1),'Taillen (cm)');
ylabel(hax(2),'Concentration en Azote (\mumole N.l.^{-1})');
grid on
