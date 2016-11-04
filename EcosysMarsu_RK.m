clear all, close all

#Initialisation :
#       Partie commune :

t0 = 1;
tf = 365*1;
dt = 0.1;
n = tf/dt;
t = t0;
T = t; #Matrice dans laquelle on stocke tous les temps

#                                  Modèle NPZ

N0 = 4; #μ.mole.N.l.^-1
P0 = 2.5; #μ.mole.N.l.^-1
Z0 = 0.5; #μ.mole.N.l.^-1

N = N0; #Abondance d'azote (en μ.mole.N.l.^-1) et à t = 0
P = P0; #Abondance d'algue (en μ.mole.N.l.^-1) et à t = 0
Z = Z0; #Abondance de Métynnis (en μ.mole.N.l.^-1) et à t = 0

Nsol = N; #Matrice dans laquelle on stocke toutes les solutions de N
Psol = P; #Matrice dans laquelle on stocke toutes les solutions de P
Zsol = Z; #Matrice dans laquelle on stocke toutes les solutions de Z




#                             Modèle croissance Marsu
L0 = 7; #cm 

L = L0; #Taille du Marsupilami (en cm) et à T0

Lsol = L; #Matrice dans laquelle on stocke toutes les solutions de L





#                     Résolution numérique des équations difféntiellelles
#                            par la méthode de Runge-Kutta



for i = 2:n #On commence la boucle à 2 car la valeur du jour 1 est celui à t = 0
  
  t = t + dt;
  
  [NRK1,PRK1,ZRK1,LRK1] = Equadif(N,P,Z,L);
  [NRK2,PRK2,ZRK2,LRK2] = Equadif(1/2*dt*NRK1,1/2*dt*PRK1,1/2*dt*ZRK1,1/2*dt*LRK1);
  [NRK3,PRK3,ZRK3,LRK3] = Equadif(1/2*dt*NRK2,1/2*dt*PRK2,1/2*dt*ZRK2,1/2*dt*LRK2);
  [NRK4,PRK4,ZRK4,LRK4] = Equadif(dt*NRK3,dt*PRK3,dt*ZRK3,dt*LRK3);
  
  dNdt = 1/6*NRK1 + 1/3*NRK2 + 1/3*NRK3 + 1/6*NRK4; 
  dPdt = 1/6*PRK1 + 1/3*PRK2 + 1/3*PRK3 + 1/6*PRK4;
  dZdt = 1/6*ZRK1 + 1/3*ZRK2 + 1/3*ZRK3 + 1/6*ZRK4;
  
  
  ##         Résolution des équations du modèle NPZ
  
  
  N = N + dt*dNdt;    #On calcule les solutions de N,P,Z pour chaque valeur de dt
  P = P + dt*dPdt;
  Z = Z + dt*dZdt;
  
  
  
  ##         Résolution de l'équation différentielle du modèle de croissance du Marsu
  
  
  if dLdt < 0 #Si le flux de croissance du Marsu est inférieur à 0, il devient nul.
    #          En effet, notre Marsu ne peut pas raptisser !
     dLdt = 0;
   
  end
  
  
  L = L + dLdt*dt;
  
  
  
  
  ##          Stockage de toutes les solutions
  
  Nsol = [Nsol ; N]; #On stocke toutes les solutions de N,P,Z,L dans des matrices
  Psol = [Psol ; P];
  Zsol = [Zsol ; Z];
  Lsol = [Lsol ; L];
  
  T = [T ; t];
end




#                             Graphiques
