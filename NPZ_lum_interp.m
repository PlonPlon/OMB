 #                                                 Fonction equations différentielles

#{
Calcul les flux de différentes équations différentielles utilisée dans le modèle de production lacustre (NPZ)
et dans celui de croissance du Marsupilami.
#}


function [dNdt,dPdt,dZdt] = NPZ(N,P,Z,IntensiteLum,t)

# Variables du modèle NPZ :


Vm = 1.0; #jour−^1 – taux maximum d’assimilation de l’azote par le phyton
ks = 1.0; #mole.N.l^−1 – constante de demi-saturation de l’assimilation de l’azote par le phyton
Rm = 1.0; #jour.^−1 – taux maximum de broûtage du phyton par les Métynnis
g = 0.2; #jour.^−1 – taux de mortalité des Métynnis
lambda = 0.2;  #(mole.N.l.^−1)^−1 constante d’Ivlev pour le broûtage ; c’est une constante empirique couramment utilisée dans les modèles pour figurer le broûtage.
e = 0.1; #jour.^−1 – taux de mortalité du phyton
gamma = 0.7; #(sans dimension) proportion d’azote assimilée par les Métynnis
#fI0 = 0.25; #fonction de l’intensité de la lumière (habituellement une décroissance expo-nentielle en fonction de la profondeur) ; dans le cas présent, la profondeur de l’étang est constante et on fera l’hypothèse, dans un premier temps au moins, que f (I 0 ) est constante et égale à 0.25Ò

# Variation de l'intensité lumineuse

#{

Pour éviter les "sauts" dans nos données, on interpole l'intensité lumineuse entre deux mesures.

#}

Inorm = IntensiteLum(1,2); #On choisi arbitrairement de dire que l'intenstité lumineuse normale est celle du jour 1.

I0 = -Inorm/log(0.75); #Calcul du paramètre I0 de la fonction f(I0), comme décrite dans l'article "NPZ Models of Plankton Dynamics: Their Construction, Coupling to Physics, and Application". On utilise la fonction f(I0) = 1 - exp(-I0/I) avec I l'intensité lumin

x = IntensiteLum(:,1)'; #On récupère tous les jours de mesures de l'intensité lumineuse dans un vecteur horizontal.

if t <= 366
  
  fI0 = [1 - exp(-IntensiteLum(1,2)/I0),1 - exp(-IntensiteLum(2,2)/I0)]; #On calcul fI0 selon la formule donnée plus haut pour le jour 1 et 366
  fI0interp = interp1(x(1:2),fI0,t); #On interpole l'intensité lumineuse entre le jour 1 et 366 en fonction du temps
  
  elseif t <= 731
    
    fI0 = [1 - exp(-IntensiteLum(2,2)/I0),1 - exp(-IntensiteLum(3,2)/I0)]; #On calcul fI0 selon la formule donnée plus haut pour le jour 1 et 366
    fI0interp = interp1(x(2:3),fI0,t); #On interpole l'intensité lumineuse entre le jour 1 et 366 en fonction du temps
    
      elseif t <= 1096
      
      fI0 = [1 - exp(-IntensiteLum(3,2)/I0),1 - exp(-IntensiteLum(4,2)/I0)]; #On calcul fI0 selon la formule donnée plus haut pour le jour 1 et 366
      fI0interp = interp1(x(3:4),fI0,t);
      
        elseif t <= 1461
        
        fI0 = [1 - exp(-IntensiteLum(4,2)/I0),1 - exp(-IntensiteLum(5,2)/I0)];
        fI0interp = interp1(x(4:5),fI0,t);
        
          elseif t <= 1826
          
          fI0 = [1 - exp(-IntensiteLum(5,2)/I0),1 - exp(-IntensiteLum(6,2)/I0)];
          fI0interp = interp1(x(5:6),fI0,t);
          
            elseif t <= 2191
            
            fI0 = [1 - exp(-IntensiteLum(6,2)/I0),1 - exp(-IntensiteLum(7,2)/I0)];
            fI0interp = interp1(x(6:7),fI0,t);
            
              elseif t <= 2556
              
              fI0 = [1 - exp(-IntensiteLum(7,2)/I0),1 - exp(-IntensiteLum(8,2)/I0)];
              fI0interp = interp1(x(7:8),fI0,t);
              
                elseif t <= 2921
                
                fI0 = [1 - exp(-IntensiteLum(8,2)/I0),1 - exp(-IntensiteLum(9,2)/I0)];
                fI0interp = interp1(x(8:9),fI0,t);
                
                  else
                  
                  fI0 = 1 - exp(-IntensiteLum(9,2)/I0);
                  fI0interp = fI0;
end

#Calcul du modèle NPZ:

dNdt = -(Vm*N)/(ks+N)*fI0*P+(1-gamma)*Z*Rm*(1-exp(-lambda*P))+e*P+g*Z; #Taux d'accroissement en azote

dPdt = (Vm*N)/(ks+N)*fI0interp*P-Z*Rm*(1-exp(-lambda*P))-e*P; #Taux d'accroisssement d'algue

dZdt = gamma*Z*Rm*(1-exp(-lambda*P))-g*Z; #Taux d'accroissement en Métynnis

endfunction
