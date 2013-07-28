# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 20:18:08 2012

@author: sebastien

Create a Cp/Ct table interpolation with sensitivity generation

"""


import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import scipy.io
mat = scipy.io.loadmat('AeroData.mat')

beta    = mat['beta_Data']
lambda_ = mat['lambda_Data']

beta -= np.min(beta)
#beta /= np.max(beta)
lambda_ -= np.min(lambda_)
#lambda_ /= np.max(lambda_)

Cp      = np.maximum(-0.0,mat['Cp_Data'])
#Cp      = mat['Cp_Data']

beta    = beta.reshape(beta.shape[1],)
lambda_ = lambda_.reshape(lambda_.shape[1],)

def MakeGrid(argx,argy):
    gridx   = [[argx[i] for j in range(len(argy))] for i in range(len(argx))]
    gridy   = [[argy[i] for j in range(len(argx))] for i in range(len(argy))]
    
    return np.array(gridx).T, np.array(gridy)

betagrid, lambdagrid = MakeGrid(beta,lambda_)

tck = interpolate.bisplrep(betagrid,lambdagrid,Cp,s=1e-2)


#spline order
knots_x     = tck[0]
knots_y     = tck[1]
checkpoints = tck[2]

p           = tck[3]
q           = tck[4]

assert(p==3)
assert(q==3)

#P is size mxn, row major representation
m           = len(knots_x)-p-1
n           = len(knots_y)-q-1
P           = tck[2].reshape(m,n)

Pline       = tck[2]

#### Compute data for derivatives
def Dx(P,knots_x,m,n,p):
    Px = np.zeros([m-1,n])
    for k1 in range(m-1):
        for k2 in range(n):
            Px[k1,k2] = p*(P[k1+1,k2] - P[k1,k2])/(knots_x[k1+p+1]- knots_x[k1+1])
    Px = Px.reshape(n*(m-1),)
    
    Ux = [[0. for k in range(p)]]
    Ux.append([knots_x[k] for k in range(p+1,m+1)])
    Ux.append([1.0 for k in range(p)])
    Ux = np.concatenate(Ux)

    return Px, Ux

def Dy(P,knots_y,m,n,q):
    Py = np.zeros([m,n-1])
    for k1 in range(m):
        for k2 in range(n-1):
            Py[k1,k2] = q*(P[k1,k2+1] - P[k1,k2])/(knots_y[k2+q+1]- knots_y[k2+1]) 
    Py = Py.reshape(m*(n-1),)

    Uy = [[0. for k in range(q)]]
    Uy.append([knots_y[k] for k in range(q+1,n+1)])
    Uy.append([1.0 for k in range(q)])
    Uy = np.concatenate(Uy)

    return Py, Uy

Px, Ux = Dx(P,knots_x,m,n,p)
Py, Uy = Dy(P,knots_y,m,n,q)

Pxx, Uxx = Dx(Px.reshape(m-1,n),      Ux, m-1 ,n   , p-1)
Pyy, Uyy = Dy(Py.reshape(m,n-1),      Uy, m   ,n-1 , q-1)
Pxy, _   = Dy(Px.reshape(m-1,n), knots_y, m-1 ,n   ,   q)

#### Write SplineData.h ####

def writeData(varname, data, vartype, fileobj):
    # Function to write a given variable
    if isinstance(data,int):
        Line = 'const ' + vartype+' ' + varname + ' = ' + str(data) + ';'
    else:
        Line = 'const ' + vartype+' ' + varname + '[' + str(len(data)) + '] = {'
        lendata = len(data)
        for k in range(lendata):
           Line += str(data[k])
           if (k < lendata-1):
              Line += ','
   
        Line += '};'

    fileobj.write(Line+'\n')
    
    
fileobj = open('SplineData.h','w')

## write variables ##
varDictionary = {'n':              [n,              'int'],
                 'm':              [m,              'int'],
                 'p':              [p,              'int'],
                 'q':              [q,              'int'],
                 'knots_x':        [knots_x,      'float'],
                 'length_knots_x': [len(knots_x),   'int'],
                 'knots_y':        [knots_y,      'float'],
                 'length_knots_y': [len(knots_y),   'int'],
                 'P':              [Pline,        'float'],
                 'length_P':       [len(Pline),     'int'],
                 'Px':             [Px,           'float'],
                 'length_Px':      [len(Px),        'int'],
                 'Py':             [Py,           'float'],
                 'length_Py':      [len(Py),        'int'],
                 'Ux':             [Ux,           'float'],
                 'length_Ux':      [len(Ux),        'int'],
                 'Uy':             [Uy,           'float'],
                 'length_Uy':      [len(Uy),        'int'],
                 'Uxx':            [Uxx,          'float'],
                 'length_Uxx':     [len(Uxx),       'int'],
                 'Uyy':            [Uyy,          'float'],
                 'length_Uyy':     [len(Uyy),       'int'],
                 'Pxx':            [Pxx,          'float'],
                 'length_Pxx':     [len(Pxx),       'int'],
                 'Pyy':            [Pyy,          'float'],
                 'length_Pyy':     [len(Pyy),       'int'],
                 'Pxy':            [Pxy,          'float'],
                 'length_Pxy':     [len(Pxy),       'int'],}

for key in varDictionary.keys():
    writeData(key, varDictionary[key][0], varDictionary[key][1], fileobj)

fileobj.close()
############################### END OF PURE PYTHON CODE ###############################




x = 5.
y = 3.

z = interpolate.bisplev(x, y, tck)
dzdx  = interpolate.bisplev(x, y, tck, dx = 1, dy = 0)
dzdy  = interpolate.bisplev(x, y, tck, dx = 0, dy = 1)
dzdx2 = interpolate.bisplev(x, y, tck, dx = 2, dy = 0)
dzdy2 = interpolate.bisplev(x, y, tck, dx = 0, dy = 2)
dzdxy = interpolate.bisplev(x, y, tck, dx = 1, dy = 1)


print "Compare spline interp", z  
print "Compare x derivative", dzdx
print "Compare y derivative", dzdy

print "Compare x2 derivative", dzdx2
print "Compare y2 derivative", dzdy2
print "Compare xy derivative", dzdxy

x,y = np.mgrid[beta[0]:beta[-1]:0.25,lambda_[0]:lambda_[-1]:0.25]
z = interpolate.bisplev(x[:,0],y[0,:],tck)
fig = plt.figure(1)
ax = fig.gca(projection='3d')
ax.plot_wireframe(betagrid,lambdagrid,Cp)
plt.title("Cp table")
ax.plot_wireframe(x,y,z,color='r')
plt.hold
plt.show()