import numpy as np
from math import sin,cos,sqrt

# Error tolerance and endpoints
a = 0.0
b = 500.0
Tol = 1.0E-8
order = 8  # Order of the system of equations
n = 4000
DIM = 8

# Initialize some parameters
ydumb = np.zeros((order), float)
y = np.zeros((order), float)
fReturn = np.zeros((order), float)
err = np.zeros((order), float)
k1 = np.zeros((order), float)
k2 = np.zeros((order), float)
k3 = np.zeros((order), float)
k4 = np.zeros((order), float)
k5 = np.zeros((order), float)
k6 = np.zeros((order), float)

# Initialize
y[0] = 0.0
y[1] = 0.0

h = (b-a)/n
hmin = h/64
hmax = h*64  # Min and max step sizes
t = a
j = 0


def f(t, y): 
    ws = 2.0                # Return RHS's force function
    gm1 = 4.0E-2 #damping coefficient mass 1
    gm2 = 1.0E-7 #damping coefficient mass 2
    m1 = 1.0 #mass 1
    m2 = 1.0 #mass 2
    k1 = 4.0 #spring constant mass 1
    k2 = 4.0 #spring constant mass 2
    w1 = sqrt(k1/m1) #natural frequency mass 1
    w2 = sqrt(k2/m2) #natural frequency mass 2
    K = 0.25 #spring constant coupling spring
    OM1 = sqrt(K/m1) #coupling frequency mass 1
    OM2 = sqrt(K/m2) #coupling frequency mass 2
    F = 0.1 #forcing amplitude
    ReF = (F/m1) * cos(ws*t) #real part of forcing function
    ImF = (F/m1) * sin(ws*t) #imaginary part of forcing function
    fReturn = np.zeros((DIM), float)
    """This are the functions to be integrated: f(i) = dy(i)/dt"""
    fReturn[0] = y[1] #x1 dot real
    fReturn[1] = -gm1*y[1] - w1*w1*y[0] + OM1*OM1*y[4] + ReF
    fReturn[2] = y[3] #x1 dot imag
    fReturn[3] = -gm1*y[3] - w1*w1*y[2] + OM1*OM1*y[6] - ImF
    fReturn[4] = y[5] #x2 dot real
    fReturn[5] = -gm2*y[5] - w2*w2*y[4] + OM2*OM2*y[0]
    fReturn[6] = y[7] #x2 dot imag
    fReturn[7] = -gm2*y[7] - w2*w2*y[6] + OM2*OM2*y[2]
    
    return fReturn


"""Open a file to save data"""
output = open('plot2.dat', 'w')


while (t < b):  # Loop over time
    if ((t+h) > b):  # Last step
        h = b - t

    f(t, y)  # Evaluate RHS's, return in fReturn

    k1[0] = h*fReturn[0]
    k1[1] = h*fReturn[1]
    for i in range(0, order):
        ydumb[i] = y[i] + k1[i]/4.0

    f(t+h/4.0, ydumb)
    k2[0] = h*fReturn[0]
    k2[1] = h*fReturn[1]
    for i in range(0, order):
        ydumb[i] = y[i] + 3.0*k1[i]/32.0 + 9.0*k2[i]/32.0

    f(t+3.0*h/8.0, ydumb)
    k3[0] = h*fReturn[0]
    k3[1] = h*fReturn[1]
    for i in range(0, order):
        ydumb[i] = y[i] + 1932.0*k1[i]/2197.0 - 7200.0*k2[i]/2197.0 + \
          7296.0*k3[i]/2197.0

    f(t+12.0*h/13.0, ydumb)
    k4[0] = h*fReturn[0]
    k4[1] = h*fReturn[1]
    for i in range(0, order):
        ydumb[i] = y[i] + 439.0*k1[i]/216.0 - 8.0*k2[i] + \
          3680.0*k3[i]/513.0 - 845.0*k4[i]/4104.0

    f(t+h, ydumb)
    k5[0] = h*fReturn[0]
    k5[1] = h*fReturn[1]
    for i in range(0, order):
        ydumb[i] = y[i] - 8.0*k1[i]/27.0 + 2.0*k2[i] - 3544.0*k3[i]/2565.0 + \
          1859.0*k4[i]/4104.0 - 11.0*k5[i]/40.0

    f(t+h/2, ydumb)
    k6[0] = h*fReturn[0]
    k6[1] = h*fReturn[1]
    for i in range(0, order):
        err[i] = abs(k1[i]/360.0 - 128.0*k3[i]/4275.0 - 2197.0*k4[i]/75240.0 +
                     k5[i]/50.0 + 2.0*k6[i]/55.0)

    if (err[0] < Tol or err[1] < Tol or h <= 2.0*hmin):  # Accept step
        for i in range(0, order):
            y[i] = y[i] + 25.0*k1[i]/216.0 + 1408.0*k3[i]/2565.0 \
              + 2197.0*k4[i]/4104.0 - k5[i]/5.0
        t = t + h
        j = j + 1

    if (err[0] == 0 or err[1] == 0):
        s = 0                              # Trap division by 0
    else:
        s = 0.84*pow(Tol*h/err[0], 0.25)   # Reduce step

    if (s < 0.75 and h > 2.0*hmin):
        h /= 2.0                          # Increase step
    else:
        if (s > 1.5 and 2.0*h < hmax):
            h *= 2.0

    print('{0: f}    {1: f}    {2: f}'.format(t, y[0], y[1], y[2]), file=output)
