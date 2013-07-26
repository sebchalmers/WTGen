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

#P is size mxn, row major representation
m           = len(knots_x)-p-1
n           = len(knots_y)-q-1
P           = tck[2].reshape(m,n)

Pline       = tck[2]

#### Compute data for derivatives
Px = np.zeros([m-1,n])
for k1 in range(m-1):
    for k2 in range(n):
        Px[k1,k2] = p*(P[k1+1,k2] - P[k1,k2])/(knots_x[k1+p+1]- knots_x[k1+1])
Px = Px.reshape(n*(m-1),)

Py = np.zeros([m,n-1])
for k1 in range(m):
    for k2 in range(n-1):
        Py[k1,k2] = q*(P[k1,k2+1] - P[k1,k2])/(knots_y[k2+q+1]- knots_y[k2+1]) 
Py = Py.reshape(m*(n-1),)

Ux = [[0. for k in range(p)]]
Ux.append([knots_x[k] for k in range(p+1,m+1)])
Ux.append([1.0 for k in range(p)])
Ux = np.concatenate(Ux)

Uy = [[0. for k in range(q)]]
Uy.append([knots_y[k] for k in range(q+1,n+1)])
Uy.append([1.0 for k in range(q)])
Uy = np.concatenate(Uy)

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
                 'length_Uy':      [len(Uy),        'int']}

for key in varDictionary.keys():
    writeData(key, varDictionary[key][0], varDictionary[key][1], fileobj)

fileobj.close()
############################### END OF PURE PYTHON CODE ###############################




########## Spline eval code -> Code generated C ###########

def basisFuncs(xi, order, U, i):
    # Minimal implementation of the Cox-deBoor formula
                       
    # Preallocate the N array
    N    = np.zeros(order+1)  # N[n] = Ni-n,_ 
    
    #Compute the first step (special branch, p=1, Ni,0 = 1)
    Den  =  U[i+1] - U[i]
    N[1] = (U[i+1] -   xi) / Den
    N[0] = (xi     - U[i]) / Den
    
    #These loops are to be unrolled for the code-generated version
    for p in range(2,order+1):
        # Update of Ni-p,p is standalone (cross arrow):
        Den  = U[i+1] - U[i-p+1]
        N[p] = (U[i+1] - xi)*N[p-1] / Den
        
        for k in range(p-1,0,-1):
            #Flat arrow
            N[k]  = (xi - U[i-k])*N[k] / Den

            #Cross arrow            
            Den   =  U[i+p-k+1] - U[i-k+1]
            N[k] += (U[i+p-k+1] - xi)*N[k-1] / Den


        # Update of Ni,p is standalone (flat arrow):
        N[0] = (xi - U[i])*N[0] / Den
    return N

def findspan(xi,knots):
    # Binary search over the knots
    index_up   =  len(knots)
    index_low  =  0
    
    if (xi <  index_low) or (xi > index_up):
        print "Out of bounds, abort"
        return [],[]
    else:
        while (index_up - index_low > 1):
            index_middle = int((index_up+index_low)/2)
            if (xi < knots[index_middle]):
                index_up = index_middle
            else:
                index_low = index_middle
                                         
    return index_low

def Pijeval(P,i,j,n):
   # Line major "matrix" P -> take out i,j entry
   return P[i*n + j]


def Blend(x_basis, y_basis, i_x, i_y, P, n):
    S = 0
    for k1 in range(len(x_basis)):
        for k2 in range(len(y_basis)):
            Pij = Pijeval(P,i_x-k1,i_y-k2,n)
            S += Pij*x_basis[k1]*y_basis[k2]        
    return S
  
def EvalSpline(x,y, knots_x, knots_y, checkpoints, Ux, Uy, Px, Py, p, q, n):
    ix  = findspan(x,knots_x)
    iy  = findspan(y,knots_y)
    
    basis_x = basisFuncs( x, p, knots_x, ix)
    basis_y = basisFuncs( y, q, knots_y, iy)             
    
    s = Blend(basis_x, basis_y, ix, iy, checkpoints, n)
    
    ix_tilde = findspan(x, Ux)
    iy_tilde = findspan(y, Uy)

    basis_x_tilde = basisFuncs( x, p-1, Ux, ix_tilde)
    basis_y_tilde = basisFuncs( y, q-1, Uy, iy_tilde)

    dsdx = Blend(basis_x_tilde, basis_y,       ix_tilde, iy,       Px, n)
    dsdy = Blend(basis_x      , basis_y_tilde, ix,       iy_tilde ,Py, n-1)

    return s, dsdx, dsdy



x = 0.
y = 7.

ix = findspan(x,knots_x)
iy = findspan(y,knots_y)
print "ix = ",ix
print "iy = ",iy

Pij = Pijeval(checkpoints,2,3,n)
print Pij

basis_x = basisFuncs(x, p, knots_x, ix)
basis_y = basisFuncs(y, q, knots_y, iy)
print "basis_x = ", basis_x
print "basis_y = ", basis_y

s = Blend(basis_x, basis_y, ix, iy, checkpoints, n)

print "s = ",s

ix_tilde = findspan(x, Ux)
iy_tilde = findspan(y, Uy)
    
print "ix_tilde = ", ix_tilde
print "iy_tilde = ", iy_tilde

basis_x_tilde = basisFuncs( x, p-1, Ux, ix_tilde)
basis_y_tilde = basisFuncs( y, q-1, Uy, iy_tilde)
   
print "basis_x_tilde = ", basis_x_tilde
print "basis_y_tilde = ", basis_y_tilde

z = interpolate.bisplev(x, y, tck)
dzdx = interpolate.bisplev(x, y, tck, dx = 1, dy = 0)
dzdy = interpolate.bisplev(x, y, tck, dx = 0, dy = 1)

S, dSdx, dSdy = EvalSpline(x, y, knots_x, knots_y, checkpoints, Ux, Uy, Px, Py, p, q, n)

print "Compare spline interp", S, z  
print "Compare x derivative", dzdx, dSdx
print "Compare y derivative", dzdy, dSdy

x,y = np.mgrid[beta[0]:beta[-1]:0.25,lambda_[0]:lambda_[-1]:0.25]
z = interpolate.bisplev(x[:,0],y[0,:],tck)
fig = plt.figure(1)
ax = fig.gca(projection='3d')
ax.plot_wireframe(betagrid,lambdagrid,Cp)
plt.title("Cp table")
ax.plot_wireframe(x,y,z,color='r')
plt.hold
plt.show()