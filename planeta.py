#!/usr/bin/env python
# -*- coding: utf-8 -*-

class Planeta(object):
    '''
    Clase queimplementa ecuaciones de movimiento contiene los metodos de
    integracion Euler explícito, Runge-Kutta 4 y Verlet. 
    Además calcula energia total
    '''


    def __init__(self, condicion_inicial, alpha=0):
        '''
        inicia las instancias de la clase.
        Ej. de uso:
        >> mercurio = Planeta([x0, y0, vx0, vy0])
        >> print(mercurio.alpha)
        >> 0.
        '''
        self.y_actual = condicion_inicial
        self.t_actual = 0.
        self.alpha = alpha


    def ecuacion_de_movimiento(self):
        '''
        Implementa la ecuación de movimiento, como sistema de ecuaciónes de
        primer orden.
        '''
        x, y, vx, vy = self.y_actual
        GM=1
        fx =lambda x,y: -(GM)*x/((x**2 + y**2)**(3/2.)) + 2*self.alpha*GM*x/((x**2 + y**2)**2)
        fy =lambda x,y: -(GM)*y/((x**2 + y**2)**(3/2.)) + 2*self.alpha*GM*y/((x**2 + y**2)**2)
        return [vx, vy, fx, fy]


    def avanza_euler(self, dt):
        '''
        Toma la condición actual del planeta y avanza su posicion y velocidad
        en un intervalo de tiempo dt usando el método de Euler explícito. El
        método no retorna nada, pero re-setea los valores de self.y_actual.
        '''
        x, y, vx, vy = self.y_actual
        fx=self.ecuacion_de_movimiento()[2] #f(x)
        fy=self.ecuacion_de_movimiento()[3] #f(y)
        vx_siguiente= vx  +  dt*fx(x,y)
        x_siguiente= x  +  dt*vx
        vy_siguiente= vy  +  dt*fy(x,y)
        y_siguiente= y  +  dt*vy
        self.y_actual = x_siguiente, y_siguiente, vx_siguiente, vy_siguiente
        pass


    def avanza_rk4(self, dt):
        '''
        Similar a avanza_euler, pero usando Runge-Kutta 4.
        '''

        fx=self.ecuacion_de_movimiento()[2] #f(x)
        fy=self.ecuacion_de_movimiento()[3] #f(y)

        x0,y0,vx0,vy0=self.y_actual
            
        l1x=dt*vx0
        l1y=dt*vy0
        
        k1x=dt*fx(x0,y0)
        k1y=dt*fy(x0,y0)
        
        l2x=dt*(vx0 + 0.5*k1x)
        l2y=dt*(vy0 + 0.5*k1y)
        
        k2x=dt*fx(x0 + 0.5*l1x,y0 + 0.5*l1y)
        k2y=dt*fy(x0 + 0.5*l1x,y0 + 0.5*l1y)
        
        l3x=dt*(vx0 + 0.5*k2x)
        l3y=dt*(vy0 + 0.5*k2y)
        
        k3x=dt*fx(x0 + 0.5*l2x,y0 + 0.5*l2y)
        k3y=dt*fy(x0 + 0.5*l2x,y0 + 0.5*l2y)
        
        l4x=dt*(vx0 + k3x)
        l4y=dt*(vy0 + k3y)
        
        k4x=dt*fx(x0 + l3x,y0 + l3y)
        k4y=dt*fy(x0 + l3x,y0 + l3y)
        
        x_siguiente=x0+(l1x + 2*l2x + 2*l3x + l4x)/6.0
        y_siguiente=y0 + (l1y + 2*l2y + 2*l3y + l4y)/6.0
        
        vx_siguiente=vx0 + (k1x + 2*k2x + 2*k3x + k4x)/6.0
        vy_siguiente=vy0 + (k1y + 2*k2y + 2*k3y + k4y)/6.0
    
        
        self.y_actual=x_siguiente,y_siguiente,vx_siguiente,vy_siguiente
        pass


    def avanza_verlet(self, dt):
        '''
        Similar a avanza_euler, pero usando Verlet.
        '''
        fx=self.ecuacion_de_movimiento()[2] #f(x)
        fy=self.ecuacion_de_movimiento()[3] #f(y)

        x0,y0,vx0,vy0=self.y_actual

        x_siguiente = x0  +  vx0*dt  +  0.5*fx(x0,y0)*dt**2
        y_siguiente = y0  +  vy0*dt  +  0.5*fy(x0,y0)*dt**2
        vx_siguiente = vx0  +  0.5*(fx(x_siguiente,y_siguiente) + fx(x0,y0))*dt
        vy_siguiente = vy0  +  0.5*(fy(x_siguiente,y_siguiente) + fy(x0,y0))*dt
        self.y_actual = x_siguiente, y_siguiente, vx_siguiente, vy_siguiente
        pass


    def energia_total(self):
        '''
        Calcula la enérgía total del sistema en las condiciones actuales.
        '''
        x, y, vx, vy = self.y_actual
        m=1
        GM=1
        E=0.5*m*(vx**2 + vy**2)-GM*m/((x**2 + y**2)**(1/2)) + self.alpha*GM*m/(((x**2 + y**2)))
        return E
