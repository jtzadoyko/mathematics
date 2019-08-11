# file:   expFit.py
# date:   Nov 15, 2016
#
# This program uses a SciPy routine to fit an exponential curve.

import numpy as np
from scipy.optimize import curve_fit

def func(x, a, b, c):
    """ Exponential function form to fit the data points """
    return a * np.exp(-b * x) + c

# read data file into numpy array
csv = np.genfromtxt('plot.dat', delimiter=',')

xdata = csv[:,0]  # first column of the data file
ydata = csv[:,1]  # second column of the data file

popt, pcov = curve_fit(func, xdata, ydata)

print("The exponential function parameters are:")
print("a =", popt[0])
print("b =", popt[1])
print("c =", popt[2])
