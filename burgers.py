#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 13:34:11 2018

@author: bonar
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import animation

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

        
c = 1.1                                                                        #signal velocity 
dt = 0.001                                                                     #timestep, adjust framerate in animation if changing 
dx = 0.01                                                                      #x step 
x_min = 0
x_max = 2*math.pi

xs = np.arange(x_min,x_max,dx)
us = np.zeros(len(xs), float)
us2 = np.zeros(len(xs), float)
for i in range(len(xs)):                                  
    us[i] = math.sin(xs[i])+1.1                                                #initial conditions
for i in range(len(xs)):                                  
    us2[i] = math.sin(xs[i])+1.1                                               

   
fig = plt.figure()                                                             #setting up figure 
ax = plt.axes(xlim=(x_min, x_max), ylim=(-2, 4))
ax.set_title(r"1D Non-linear Burgers' Equation", fontsize=20)
ax.set_ylabel(r'$u(x,t)$', fontsize = 20)
ax.set_xlabel(r'$x$', fontsize=16)
plt.plot(xs, us, 'r--', alpha = 0.5)                                           #plot inital condition 
red_patch = mpatches.Patch(color='red', label='Conservative')
blue_patch = mpatches.Patch(color='blue', label='Advective')
plt.legend(handles=[red_patch, blue_patch])

plotlays, plotcols = [2], ["b-","r-"]
lines = []
for index in range(2):
    lobj = ax.plot([],[],plotcols[index],lw=2,alpha = 0.5)[0]
    lines.append(lobj)

def cent(u, i):
    return(-c*(dt/(2*dx))*(u[i+2]-u[i])+u[i+1])                                #1D linear advection dicretised equations
def forw(u, i):
    return(-c*(dt/dx)*(u[i+2]-u[i+1])+u[i+1])
def back(u, i):
    return(-c*(dt/dx)*(u[i+1]-u[i])+u[i+1])    
def cent_b(u, i):
    return(-u[i+1]*(dt/(2*dx))*(u[i+2]-u[i])+u[i+1])                           #1D burger equations with different approximations of spatial integral
def forw_b(u, i):
    return(-u[i+1]*(dt/dx)*(u[i+2]-u[i+1])+u[i+1])
def back_b(u, i):
    return(-u[i+1]*(dt/dx)*(u[i+1]-u[i])+u[i+1])
def back(u, i):
    return(-c*(dt/dx)*(u[i+1]-u[i])+u[i+1])                                    #backward difference with advective form 
def back2(u, i):
    return(-(1/2)*(dt/dx)*(u[i+1]**2-u[i]**2)+u[i+1])                          #backward difference with conservative form
    
def prop(j, dif):                                                              #function for propogating u array with given equation (dif) 
    for x in range(j):
        u = us                                                                 #set up array with current u values, and add points for boundary conditions
        u = np.insert(u, 0, u[-1])                                             #inflow as u = 1, outflow has u = 0 (maybe) or periodic for burgers
        u = np.append(u, u[1])
        for i in range(len(us)):  
            us[i]= dif(u, i)                                                   #current u value in u array is indexed by i+1 (u[0] is us[-1])
    return(us)
    
def prop2(j, dif):       
    for x in range(j):
        u2 = us2                                                               #set up array with current u values, and add points for boundary conditions
        u2 = np.insert(u2, 0, u2[-1])                                          #inflow as u = 1, outflow has u = 0 (maybe) or periodic for burgers
        u2 = np.append(u2, u2[1])
        for i in range(len(us2)):  
            us2[i]= dif(u2, i)                                                 #current u value in u array is indexed by i+1 (u[0] is us[-1])
    return(us2)
    
    
def init():
    for line in lines:
        line.set_data([],[])                                                   #initialisation of lines for animation 
    return lines

def animate(j):                                                                #animation function, calculates u arrays at j time step using prop 
    xlist = [xs, xs]                                                           #change dif in prop() and prop2() to see other approximations of spatial derivative 
    ylist = [prop(j, back_b), prop2(j, back2)]
    for lnum,line in enumerate(lines):
        line.set_data(xlist[lnum], ylist[lnum])
    return(line,)


anim= animation.FuncAnimation(fig, animate, init_func=init, frames = 1000, 
                               interval=100, blit = False, repeat=False)
#anim.save('line.gif', dpi=80, writer='imagemagick')
plt.show()


