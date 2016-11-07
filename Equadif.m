## Fonction équation differentielle
## Calcul les flux de différentes équations différentielles utilisée dans le modèle de production lacustre (NPZ) et dans celui de croissance du Marsupilami.


function [dNdt,dPdt,dZdt,dLdt] = Equadif(N,P,Z,L)

# Variables du modèle NPZ :


Vm = 1.0; #jour−^1 – taux maximum d’assimilation de l’azote par le phyton
ks = 1.0; #mole.N.l^−1 – constante de demi-saturation de l’assimilation de l’azote par le phyton
Rm = 1.0; #jour.^−1 – taux maximum de broûtage du phyton par les Métynnis
g = 0.2; #jour.^−1 – taux de mortalité des Métynnis
lambda = 0.2;  #(mole.N.l.^−1)^−1 constante d’Ivlev pour le broûtage ; c’est une constante empirique couramment utilisée dans les modèles pour figurer le broûtage.
e = 0.1; #jour.^−1 – taux de mortalité du phyton
gamma = 0.7; #(sans dimension) proportion d’azote assimilée par les Métynnis
fI0 = 0.25; #fonction de l’intensité de la lumière (habituellement une décroissance expo-nentielle en fonction de la profondeur) ; dans le cas présent, la profondeur de l’étang est constante et on fera l’hypothèse, dans un premier temps au moins, que fI0 ) est constante et égale à 0.25Ò

# Variables du modèle de croissance du Marsu

Linf = 120; #Longueur théorique Max du Marsu en cm

K = 0.001; # Constante de croissance en jour^-1

# Calcul du modèle NPZ:

dNdt = -(Vm*N)/(ks+N)*fI0*P+(1-gamma)*Z*Rm*(1-exp(-lambda*P))+e*P+g*Z; #Taux d'accroissement en azote

dPdt = (Vm*N)/(ks+N)*fI0*P-Z*Rm*(1-exp(-lambda*P))-e*P; #Taux d'accroisssement d'algue

dZdt = gamma*Z*Rm*(1-exp(-lambda*P))-g*Z; #Taux d'accroissement en Métynnis

# Calcul du modèle de croissance marsupilami

dLdt = K*(Linf*Z - L);


endfunction
