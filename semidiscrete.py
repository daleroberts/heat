# solve the one dimensional heat equation on a domain
# using a semidiscrete method.
#
# Dale Roberts <dale.o.roberts@gmail.com>

import numpy as np
from numpy import pi, cos, sin

# space and time mesh

x, h = np.linspace(0.0, 1.0, 256, retstep=True)
t, p = np.linspace(0.0, 2.0, 256, retstep=True)

# initial condition

u0 = x*sin(pi*x)

# Laplacian

n = x.shape[0]

A = np.zeros((n,n))
A.flat[::(n+1)] = -2.0/h**2
A.flat[1::(n+1)] = A.flat[n::(n+1)] = 1.0/h**2

# setup U'(t) + F(U(t)) = 0

def F(u,s):
    w = np.dot(A,u)
    # set boundary conditions
    w[0] = sin(2*pi*s)
    w[-1] = 0.0
    return w

# solve system of ODEs

import scipy.integrate as integrators

U = integrators.odeint(F,u0,t)

# 2d plot for matrix U

import pylab

pylab.plot(x,U[0],'r')
for i in range(1, U.shape[0]):
    pylab.plot(x,U[i],'g')

# 3d plot where x and t are vectors containing the time and space points

x, h = np.linspace(0.0, 1.0, 256, retstep=True)
t, p = np.linspace(0.0, 2.0, 256, retstep=True)

xm, tm = np.meshgrid(x, t)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

Axes3D(pylab.figure()).plot_surface(xm, tm, U, cmap=cm.jet, alpha=0.8)

plt.show()
