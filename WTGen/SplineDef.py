# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 20:18:08 2012

@author: sebastien

Create a Cp/Ct table interpolation with sensitivity generation

"""


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.io
from SplineGen import *

# Cp table
mat = scipy.io.loadmat('CpData.mat')

beta    = mat['beta_Data']
lambda_ = mat['lambda_Data']
beta    = beta.reshape(beta.shape[1],)
lambda_ = lambda_.reshape(lambda_.shape[1],)
    
Cp      = mat['Cp_Data']

tck, betagrid, lambdagrid = SplineGen(beta,lambda_,Cp,'Cp')

#Ct table
mat = scipy.io.loadmat('CtData.mat')

beta    = mat['beta_Data']
lambda_ = mat['lambda_Data']
beta    = beta.reshape(beta.shape[1],)
lambda_ = lambda_.reshape(lambda_.shape[1],)
    
Ct      = mat['Ct_Data']

tck, betagrid, lambdagrid = SplineGen(beta,lambda_,Ct,'Ct')

############################### END OF PURE PYTHON CODE ###############################

#x = 5.
#y = 3.
#
#z = interpolate.bisplev(x, y, tck)
#dzdx  = interpolate.bisplev(x, y, tck, dx = 1, dy = 0)
#dzdy  = interpolate.bisplev(x, y, tck, dx = 0, dy = 1)
#dzdx2 = interpolate.bisplev(x, y, tck, dx = 2, dy = 0)
#dzdy2 = interpolate.bisplev(x, y, tck, dx = 0, dy = 2)
#dzdxy = interpolate.bisplev(x, y, tck, dx = 1, dy = 1)
#
#
#print "Compare spline interp", z  
#print "Compare x derivative", dzdx
#print "Compare y derivative", dzdy
#
#print "Compare x2 derivative", dzdx2
#print "Compare y2 derivative", dzdy2
#print "Compare xy derivative", dzdxy
#
#x,y = np.mgrid[beta[0]:beta[-1]:0.25,lambda_[0]:lambda_[-1]:0.25]
#z = interpolate.bisplev(x[:,0],y[0,:],tck)
#fig = plt.figure(1)
#ax = fig.gca(projection='3d')
#ax.plot_wireframe(betagrid,lambdagrid,Cp)
#plt.title("Cp table")
#ax.plot_wireframe(x,y,z,color='r')
#plt.hold
#plt.show()