"""
Created on Fri Nov 16 20:18:08 2012

@author: sebastien

Implements the closed-lopp robust NMPC scheme
Covariance is computed by Lyapunov integration,
the feedback is the NE controller computed from the shifted trajectories

"""

from casadi import*
import rawe 
from GNData import*

def makeDae():
    dae = rawe.dae.Dae()

    #Define some constants
    Mgfactor = 1e4
   
    rho = 1.23
    R = 63.
    eta = 0.95
    hH = 90.
   
    mT = 347460.
    mN = 240000.
    mB = 17740.
    mH = 56780.
    xi = 0.7
    iratio = 0.010309278350515
    JH = 115926.
    JB = 11776047.
    JG = 534.116
    f0 = 0.32
    ds = 0.01
   
    mTe = mT/4. + mN + mH + 3*mB; 
    cT = 4*pi*mTe*ds*f0;
    J = JH + 3*JB + JG/iratio/iratio;
    kT = mTe*(2*pi*f0)*(2*pi*f0);

    [omega,x,dx,Tg,beta,dbeta]    = dae.addX( ["omega", "x", "dx", "Tg", "beta", "dbeta"] )
    [dTg,ddbeta, v0, slack_omega] = dae.addU( ["dTg", "ddbeta" , "v0", "slack_omega"] )
    
    # some extra outputs, why not
    dae['EPower'] = omega*Tg

    
    vrel    = v0 - dx
    lambda_ = omega*R/vrel
    
    lambda_scaled = (lambda_ - lambda_star)/lambda_scale;
    beta_scaled   = (   beta -   beta_star)/beta_scale;
                
    zeta = []
    for i in range(int(Nbasis/2.)):
        zeta.append( lambda_scaled**(i+1) )
        zeta.append(  beta_scaled**(i+1) )
    
    zeta = veccat(zeta)

    Cp = Cpmax - 0.25*mul(zeta.T,mul(M.T,mul(M,zeta)))

    p00Ct =      0.1318    
    p10Ct =     -0.1833      
    p01Ct =     -0.0552    
    p20Ct =      0.1019  
    p11Ct =     0.04685  
    p02Ct =    0.005713   
    p30Ct =    -0.01414  
    p21Ct =    -0.01131 
    p12Ct =   -0.002921 
    p03Ct =   -0.000385   
    p40Ct =    0.000816 
    p31Ct =   0.0008633   
    p22Ct =   0.0003525   
    p13Ct =   0.0001456   
    p04Ct =   9.492e-06   
    p50Ct =   -1.71e-05    
    p41Ct =  -2.306e-05    
    p32Ct =    -1.3e-05    
    p23Ct =  -1.105e-05   
    p14Ct =  -2.064e-06     
    p05Ct =  -6.932e-08  

    Ct = p00Ct + p10Ct*lambda_ + p01Ct*beta + p20Ct*lambda_*lambda_ + p11Ct*lambda_*beta + p02Ct*beta*beta + p30Ct*lambda_*lambda_*lambda_ + p21Ct*lambda_*lambda_*beta + p12Ct*lambda_*beta*beta + p03Ct*beta*beta*beta + p40Ct*lambda_*lambda_*lambda_*lambda_ + p31Ct*lambda_*lambda_*lambda_*beta + p22Ct*lambda_*lambda_*beta*beta + p13Ct*lambda_*beta*beta*beta + p04Ct*beta*beta*beta*beta + p50Ct*lambda_*lambda_*lambda_*lambda_*lambda_ + p41Ct*lambda_*lambda_*lambda_*lambda_*beta + p32Ct*lambda_*lambda_*lambda_*beta*beta + p23Ct*lambda_*lambda_*beta*beta*beta + p14Ct*lambda_*beta*beta*beta*beta + p05Ct*beta*beta*beta*beta*beta
 
    Taero = 0.5*rho*pi*R*R*R*Cp*vrel*vrel/lambda_
    Faero = 0.5*rho*pi*R*R*Ct*vrel*vrel;

    domega = (Taero - Mgfactor*Tg/iratio)/J
    ddx    = (Faero - cT*dx - kT*x)/mTe
    
    # residual
    dae.setResidual([dae.ddt('omega') -  domega,
                     dae.ddt(    'x') -  dx,
                     dae.ddt(   'dx') -  ddx,
                     dae.ddt(   'Tg') -  dTg,
                     dae.ddt( 'beta') -  dbeta,
                     dae.ddt('dbeta') -  ddbeta   ])

    return dae

print "Make DAE"
fuckit = makeDae()