"""
author: jtzadoyko
created on Tue Oct 17 11:30:39 2017
"""

import numpy as np
from scipy.integrate import ode
from math import sin, cos, sqrt

DIM = 8 #number of equations for system of ODEs

"""PARAMETERS OF DIFFERENTIAL EQUATIONS"""
gm1 = 4.0E-2 #damping coefficient mass 1
gm2 = 1.0E-7 #damping coefficient mass 2
m1 = 1 #mass 1
m2 = 1 #mass 2
k1 = 4 #spring constant mass 1
k2 = 4 #spring constant mass 2
w1 = sqrt(k1/m1) #natural frequency mass 1
w2 = sqrt(k2/m2) #natural frequency mass 2
K = 0.4 #spring constant coupling spring
OM1 = sqrt(K/m1) #coupling frequency mass 1
OM2 = sqrt(K/m2) #coupling frequency mass 2
F = 0.1 #forcing amplitude



def system(t, y, ws):
    
    fReturn = np.zeros((DIM), float)
    """This are the functions to be integrated: f(i) = dy(i)/dt"""
    fReturn[0] = y[1] #x1 dot real
    fReturn[1] = -2*gm1*y[1] - w1*w1*y[0] + OM1*OM1*y[4] + (F/m1) * cos(ws*t)
    fReturn[2] = y[3] #x1 dot imag
    fReturn[3] = -2*gm1*y[3] - w1*w1*y[2] + OM1*OM1*y[6] - (F/m1) * sin(ws*t)
    fReturn[4] = y[5] #x2 dot real
    fReturn[5] = -2*gm2*y[5] - w2*w2*y[4] + OM2*OM2*y[0]
    fReturn[6] = y[7] #x2 dot imag
    fReturn[7] = -2*gm2*y[7] - w2*w2*y[6] + OM2*OM2*y[2]
    
    return fReturn

"""runs the RK853 method on the system function
   gets passed the forcing function's frequency"""
def dop853(ws):
    # Create an `ode` instance to solve the system of differential
    # equations defined by `system`, and set the solver method to 'dop853'.
    solver = ode(system)
    solver.set_integrator('dop853')
    
    # Give the value of omega to the solver. This is passed to
    # `system` when the solver calls it.
    solver.set_f_params(ws)
    
    # Set the initial values, y is the inital conditions.
    t0 = 0.0
    y = [0]*DIM
    solver.set_initial_value(y, t0)
    
    # Create the array `t` of time values at which to compute
    # the solution, and create an array to hold the solution.
    # Put the initial value in the solution array.
    t1 = 500
    N = 3000
    t = np.linspace(t0, t1, N)
    sol = np.empty((N, 8))
    sol[0] = y
    
    # Repeatedly call the `integrate` method to advance the
# solution to time t[k], and save the solution in sol[k].
    k = 1
    while solver.successful() and solver.t < t1:
        solver.integrate(t[k])
        sol[k] = solver.y
        k += 1    
    length = len(sol)
    power = np.zeros(length,float)
    
    """calculates Re[power]"""
    for i in range(length):    
        power[i] = sol[i][1] * F * cos(ws*t[i])
        
    h = t1/N # step size    
    
    """numeric integration to find average power over region"""
    avgP = 0
    avgP = avgP + 3/8*h*(power[2799]+power[2999])+7/6*h*(power[2800]+power[2998])
           +23/24*h*(power[2801]+power[2997])
    for i in range(2802,2997):
        avgP = avgP + h*power[i]
        
    
       
    
    return  t, avgP

def main():

    wsCons = w1-0.3 # starting point for iterations of ws

    
    """Write the output into a files"""
    output = open('Numerical EIT.dat', 'w')
    """Changes the value of the forcing frequency"""
    for s in range(0,1800):
        ws = wsCons + (s/3000)
        
            
        avgP = dop853(ws)[1] #average power of steady state
        detuned = ws - w1 # gets the detuned frequency for the system
        
        print('{0: f}    {1: f}'.format(detuned, avgP), file=output) #prints tsv
    
    """Whatever is open should be closed"""
    output.close()

main()







