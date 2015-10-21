# -*- coding: utf-8 -*-
'''
Script que utiliza Verlet para calcular y graficar aproximadamente
5 orbitas. Además grafica energia vs tiempo.
'''

import numpy as np
import matplotlib.pyplot as plt
from Planeta import Planeta


condicion_inicial = [10, 0, 0, 0.3] #x, y, vx, vy

p = Planeta(condicion_inicial) #se establece condicio inicial como actuales

dt=1
t_total=900 #aproximadamente 5 orbitas
# se crean listas vacias de x, y, Energia
x=[]
y=[]
Energia=[]
#Se agregan condiciones iniciales
x.append(condicion_inicial[0])
y.append(condicion_inicial[1])
Energia.append(p.energia_total())
t=1 #ya se realizo el primer calculo
while t<=t_total:
    p.avanza_verlet(dt)
    x_siguiente, y_siguiente, vx_siguiente, vy_siguiente = p.y_actual
    x.append(x_siguiente)
    y.append(y_siguiente)
    Energia.append(p.energia_total())
    t+=dt

arreglo_de_tiempos=np.linspace(0,t_total,t_total/dt+1)

#Plotea orbita (y vs x)
plt.figure(1)
plt.plot(x,y)
plt.title('Orbita metodo Verlet ($\\alpha=0$)')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.show()

#Plotea Energia vs tiempo
plt.figure(2)
plt.clf
plt.plot(arreglo_de_tiempos,Energia)
plt.title('Grafico Energia vs tiempo')
plt.xlabel('E')
plt.ylabel('tiempo')
plt.grid(True)
plt.show()