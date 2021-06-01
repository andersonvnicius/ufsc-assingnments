# -*- coding: utf-8 -*-
"""
--
Calculating chemical equilibrium of the combustion of methanol(CH3OH) with stardard athmospheric air
--
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

r = 14
n_of_steps = 100

#Plot vectors
T = np.ones(n_of_steps)
P = np.ones(n_of_steps)
U = np.ones(n_of_steps)
S = np.ones(n_of_steps)
phi = np.ones(n_of_steps)
gibbs = np.ones(n_of_steps)
heimh = np.ones(n_of_steps)
XH2O = np.ones(n_of_steps)
XCO2 = np.ones(n_of_steps)
XCO = np.ones(n_of_steps)
XN2 = np.ones(n_of_steps)

#output srings
#leg = 'r_' +str(r) + '_FP1atm'; print(r)
leg = 'r_' +str(r) + '_FT25C'; print(r)

#input data
intake_temperature = 273. + 25.
intake_pressure = 101325.

#loading element
meth = ct.Solution('gri30.xml')
meth.TP = 298,101325
meth.X = {'CH3OH':1,'O2':1.5 * 1,'N2':1.5 *3.76}

#k constant value
k = meth.cp/meth.cv

#calculating temperature and pressure @spark
spark_pressure = intake_pressure*r**(k)
spark_temperature = intake_temperature*r**(k-1)

for j in range (0,4):
    
    #new pressure or temp value
    intake_temperature = intake_temperature #+ 10
    intake_pressure = intake_pressure + 15198.75
    
    #new pressure and temperature @spark
    spark_temperature = intake_temperature*r**(k-1)
    spark_pressure = intake_pressure*r**(k)
    
    for i in range (0,n_of_steps):
        
        #new value of phi and reseting the element
        phi[i] = 0.6 + i * 1. /n_of_steps
        meth.X = np.zeros(53)
        		
        #new element parameters
        meth.X = {'CH3OH':1,'O2':1.5/phi[i] * 1,'N2':1.5/phi[i] *3.76}
        meth.TP = spark_temperature , spark_pressure
        		
        #calculating equilibrium
        meth.equilibrate('UV')
        
        #ouput data
        T[i] = meth.T
        P[i] = meth.P /1E6
        U[i] = meth.u
        S[i] = meth.s
        gibbs[i] = meth.g /1E6
        heimh[i] = (U[i] - T[i]*S[i]) /1E6
        XH2O[i] = meth['H2O'].X
        XCO2[i] = meth['CO2'].X
        XCO[i] = meth['CO'].X
        XN2[i] = meth['N2'].X

    #max values
    print("FT=", intake_temperature - 273)
    print("FP=", intake_pressure)
    print("phi=", phi[np.argmax(T)], "T=", np.amax(T))
    print("phi=", phi[np.argmax(P)], "P=", np.amax(P))
    print("phi=", phi[np.argmin(gibbs)], "gibbs=", np.amin(gibbs))
    print("phi=", phi[np.argmin(gibbs)], "hholtz=", np.amin(gibbs))
    #print()
    
    #plot
    plt.figure(1, figsize = (6,4))
    plt.plot(phi,P)
    plt.xlabel('phi'); plt.ylabel('P[MPa]');
    plt.figure(2, figsize = (6,4)) 
    plt.plot(phi,T)
    plt.xlabel('phi'); plt.ylabel('T[K]');
    plt.figure(3, figsize = (6,4))
    plt.plot(phi,gibbs)
    plt.xlabel('phi'); plt.ylabel('Gibbs[MJ]');
    plt.figure(4, figsize = (6,4))
    plt.plot(phi,heimh)
    plt.xlabel('phi'); plt.ylabel('Heimholtz[MJ]');
    plt.figure(5, figsize = (6,4))
    plt.plot(phi,XCO2)
    plt.xlabel('phi'); plt.ylabel('XCO2');
    plt.figure(6, figsize = (6,4))
    plt.plot(phi,XCO)
    plt.xlabel('phi'); plt.ylabel('XCO');
    plt.figure(7, figsize = (6,4))
    plt.plot(phi,XN2)
    plt.xlabel('phi'); plt.ylabel('XN2');
    plt.figure(8, figsize = (6,4))
    plt.plot(phi,XH2O)
    plt.xlabel('phi'); plt.ylabel('XH2O');

"""    
#save plot FT
dpi = 500
plt.figure(1); plt.title('r = %s' %str(r)); plt.legend(['Ti=25C','Ti=35C','Ti=45C','Ti=55C'], loc=0); plt.savefig("P_%s.png" % leg , dpi=dpi)
plt.figure(2); plt.title('r = %s' %str(r)); plt.legend(['Ti=25C','Ti=35C','Ti=45C','Ti=55C'], loc=0); plt.savefig("T_%s.png" % leg , dpi=dpi)
plt.figure(3); plt.title('r = %s' %str(r)); plt.legend(['Ti=25C','Ti=35C','Ti=45C','Ti=55C'], loc=0); plt.savefig("gibbs_%s.png" % leg , dpi=dpi)
plt.figure(4); plt.title('r = %s' %str(r)); plt.legend(['Ti=25C','Ti=35C','Ti=45C','Ti=55C'], loc=0); plt.savefig("heimholtz_%s.png" % leg , dpi=dpi)
plt.figure(5); plt.title('r = %s' %str(r)); plt.legend(['Ti=25C','Ti=35C','Ti=45C','Ti=55C'], loc=0); plt.savefig("XCO2_%s.png" % leg , dpi=dpi)
plt.figure(6); plt.title('r = %s' %str(r)); plt.legend(['Ti=25C','Ti=35C','Ti=45C','Ti=55C'], loc=0); plt.savefig("XCO_%s.png" % leg , dpi=dpi)
plt.figure(7); plt.title('r = %s' %str(r)); plt.legend(['Ti=25C','Ti=35C','Ti=45C','Ti=55C'], loc=0); plt.savefig("XN2_%s.png" % leg , dpi=dpi)
plt.figure(8); plt.title('r = %s' %str(r)); plt.legend(['Ti=25C','Ti=35C','Ti=45C','Ti=55C'], loc=0); plt.savefig("XH2O_%s.png" % leg , dpi=dpi)

"""
#save plot FT
dpi = 500
plt.figure(1); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("P_%s.png" % leg , dpi=dpi)
plt.figure(2); plt.title('r = %s' %str(r)); plt.legend(['Pi=1bar','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("T_%s.png" % leg , dpi=dpi)
plt.figure(3); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("gibbs_%s.png" % leg , dpi=dpi)
plt.figure(4); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("Heimholtz_%s.png" % leg , dpi=dpi)
plt.figure(5); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("XCO2_%s.png" % leg , dpi=dpi)
plt.figure(6); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("XCO_%s.png" % leg , dpi=dpi)
plt.figure(7); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("XH2O_%s.png" % leg , dpi=dpi)
plt.figure(8); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("XN2_%s.png" % leg , dpi=dpi)