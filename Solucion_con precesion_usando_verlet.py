#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Script que utiliza Verlet para calcular y graficar aproximadamente
30 orbitas considerando alpha distinto de 0. Además grafica energia vs tiempo.
'''

import numpy as np
import matplotlib.pyplot as plt
from planeta import Planeta

condicion_inicial  =  [10, 0, 0, 0.3]
p  =  Planeta(condicion_inicial, alpha = 10**(-2.257))

# periodo aproximado por tercera ley de kepler
T = np.sqrt(4*np.pi**2*8**3) #  =  142 aproximadamente
t_total = 30*T
# se crean arreglos vacios de x, y, Energia
x = np.array([])
y = np.array([])
Energia = np.array([])
#Se agregan condiciones iniciales
x = np.append(x,10)
y = np.append(y,0)
Energia = np.append(Energia,p.energia_total())

dt = 1
t = 1
while t<= t_total:
    p.avanza_verlet(dt)
    xf, yf, vxf, vyf  =  p.y_actual
    x = np.append(x,xf)
    y = np.append(y,yf)
    Energia = np.append(Energia,p.energia_total())
    t+= dt
    
arreglo_de_tiempos = np.linspace(0,int(t_total),int(t_total)/dt+1)

"""calculo de velocidad angular"""

# se detecta inicio de la ultima orbita
t_inicio_ulima_orbita = int(t_total-T)

#se obtienen posiciones (x,y) desde t_inicio_ulima_orbita hasta el final
x_ultima_orbita = x[t_inicio_ulima_orbita+1:]
y_ultima_orbita = y[t_inicio_ulima_orbita+1:]

# calculo de distancia de cada punto anterior al foco
distancia = np.sqrt(x_ultima_orbita**2+y_ultima_orbita**2)

#se obtiene el indice de la maxima distancia (corresponde a afelio)
I=np.argmax(distancia)

#angulo entre vector foco-afelio inicial y foco-afelio final
angulo  =  np.arctan(np.abs(y_ultima_orbita[I]/x_ultima_orbita[I]))

#tiempo en que el planeta esta en afelio final
Tiempo_afelio_final= t_inicio_ulima_orbita + arreglo_de_tiempos[I]

# velocidad angular angulo / ( tiempo cuando se esta en afelio final )
velocidad_angular  =  angulo/Tiempo_afelio_final

# mostrar velocidad angular
print 'velocidad angular  =  ' + str(velocidad_angular)

plt.figure(1)
plt.clf
plt.plot(x,y)
plt.title('Orbita metodo verlet ($\\alpha \\neq 0$)')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.draw()
plt.show()

plt.figure(2)
plt.clf
plt.plot(arreglo_de_tiempos,Energia)
plt.title('Grafico Energia total vs tiempo')
plt.xlabel('tiempo')
plt.ylabel('Energia')
plt.grid(True)
plt.draw()
plt.show()