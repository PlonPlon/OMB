clear all, close all


#Initialisation :

t0 = 0;
tf = 365;
dt = 0.1;
n = tf/dt;

N0 = 4;      #μ.mole.N.l.^-1
P0 = 2.5;    #μ.mole.N.l.^-1
Z0 = 0.5;    #μ.mole.N.l.^-1

t = t0;
N = N0;     #Abondance d'azote (en μ.mole.N.l.^-1) et à t = 0
P = P0;     #Abondance d'algue (en μ.mole.N.l.^-1) et à t = 0
Z = Z0;     #Abondance de Métynnis (en μ.mole.N.l.^-1) et à t = 0

T = t;      #Matrice dans laquelle on stocke tous les temps
Nsol = N;   #Matrice dans laquelle on stocke toutes les solutions de N
Psol = P;   #Matrice dans laquelle on stocke toutes les solutions de P
Zsol = Z;   #Matrice dans laquelle on stocke toutes les solutions de Z

for i = 1:n
  
  t = t + dt;
  
  [NRK1,PRK1,ZRK1] = NPZ(N,P,Z);
  [NRK2,PRK2,ZRK2] = NPZ(1/2*dt*NRK1,1/2*dt*PRK1,1/2*dt*ZRK1);
  [NRK3,PRK3,ZRK3] = NPZ(1/2*dt*NRK2,1/2*dt*PRK2,1/2*dt*ZRK2);
  [NRK4,PRK4,ZRK4] = NPZ(dt*NRK3,dt*PRK3,dt*ZRK3);
  
  N = N + dt*(1/6*NRK1 + 1/3*NRK2 + 1/3*NRK3 + 1/6*NRK4);    #On calcule les solutions de N,P,Z pour chaque valeur de dt
  P = P + dt*(1/6*PRK1 + 1/3*PRK2 + 1/3*PRK3 + 1/6*PRK4);
  Z = Z + dt*(1/6*ZRK1 + 1/3*ZRK2 + 1/3*ZRK3 + 1/6*ZRK4);
  
  Nsol = [Nsol ; N];  #On stocke toutes les solutions de N,P,Z dans des matrices
  Psol = [Psol ; P];
  Zsol = [Zsol ; Z];
  
  T = [T ; t];
end


figure(1) #Toutes les courbes sur un même graphique

plot(T,Nsol,'b',T,Psol,'g',T,Zsol,'r')
legend('Azote','Algue','Metynnis')
xlabel('Temps (jour)')
ylabel('Concentration en Azote (u.mole.N.l.^{-1})')
grid on