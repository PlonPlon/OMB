## Fonction Croissance Marsu
## Calcule le flux de l'équations différentielles utilisée dans le modèle de croissance du Marsupilami.


function [dLdt] = CroissMarsu(t,Z,L,Data)


# Variables du modèle de croissance du Marsu

Linf = 120; #Longueur théorique Max du Marsu en cm

Tref = 37 + 273; #Température de référence en °K

TA = 8000; #Température d'Arrhénius

Kref = 0.03*30; #En jour^(-1)



#Extraction des températures :

AmpT=abs(min(Data(:,2))-max(Data(:,2))); #Amplitude thermique des températures en°C

MoyT=mean(Data(:,2)); #Moyenne des températures en °C


# Calcul du modèle de la correction de température :

T = 273 + MoyT + 0.5*AmpT*cos((2*pi*t)/365); #Temperature calculée en K

cT = exp(TA/Tref - TA/T); #Correction de température sans unité


# Calcul du modèle de croissance marsupilami

dLdt = cT*Kref*(Linf*Z-L);



endfunction
