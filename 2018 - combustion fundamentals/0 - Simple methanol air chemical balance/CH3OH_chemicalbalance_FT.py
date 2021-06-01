# -*- coding: utf-8 -*-
"""
--
Calculating chemical equilibrium of the combustion of methanol(CH3OH) with stardard athmospheric air
--
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

n_of_steps = 100

#input data
r = 16 #engine compression ratio
track_temp_c = 25.
intake_pressure = 101325.

#calc spark pressure
meth = ct.Solution('gri30.xml')
meth.X = {'CH3OH':1,'O2':1.5 * 1,'N2':1.5 *3.76}
k = meth.cp/meth.cv;
spark_pressure = intake_pressure*r**(k-1)
intake_temperature = 273. + track_temp_c

#creating plot vectors
T = np.ones(n_of_steps)
P = np.ones(n_of_steps)
dP = np.ones(n_of_steps)
rho = np.ones(n_of_steps)
phi = np.ones(n_of_steps)
afr = np.ones(n_of_steps)
XCH3OH = np.ones(n_of_steps)
XH2O = np.ones(n_of_steps)
XCO2 = np.ones(n_of_steps)
XCO = np.ones(n_of_steps)
XN2 = np.ones(n_of_steps)
#leg = 'r_' +str(r) + '_FP1atm'; print(r)
leg = 'r_' +str(r) + '_FT25C'; print(r)

for j in range (0,4):
    fuel_temperature = intake_temperature #+ 10*j 
    fuel_pressure = intake_pressure + 15198.75*j 
    spark_pressure = fuel_pressure*r**(k-1)
    for i in range (0,n_of_steps):
        phi[i] = 0.6 + i * 1/n_of_steps
        meth.X = np.zeros(53)
        meth.X = {'CH3OH':1,'O2':1.5/phi[i] * 1,'N2':1.5/phi[i] *3.76}
        meth.TP = fuel_temperature , spark_pressure
        meth.equilibrate('UV')
        #output variables
        T[i] = meth.T
        P[i] = meth.P /1E6
        XCH3OH[i] = meth['CH3OH'].X
        XH2O[i] = meth['H2O'].X
        XCO2[i] = meth['CO2'].X
        XCO[i] = meth['CO'].X
        XN2[i] = meth['N2'].X
    
    #max values
    #print("FT=", fuel_temperature - 273)
    print("FP=", fuel_pressure)
    print("phi=", phi[np.argmax(T)], "T=", np.amax(T))
    print("phi=", phi[np.argmax(P)], "P=", np.amax(P))
    print()

"""    
    #plot
    plt.figure(1, figsize = (8,4))
    plt.plot(phi,P)
    plt.xlabel('phi'); plt.ylabel('P[MPa]');
    plt.figure(2, figsize = (8,4)) 
    plt.plot(phi,T)
    plt.xlabel('phi'); plt.ylabel('T[K]');
    plt.figure(3, figsize = (8,4))
    plt.plot(phi,XCH3OH)
    plt.xlabel('phi'); plt.ylabel('XCH3OH');
    plt.figure(4, figsize = (8,4))
    plt.plot(phi,XCO2)
    plt.xlabel('phi'); plt.ylabel('XCO2');
    plt.figure(5, figsize = (8,4))
    plt.plot(phi,XCO)
    plt.xlabel('phi'); plt.ylabel('XCO');
    plt.figure(6, figsize = (8,4))
    plt.plot(phi,XN2)
    plt.xlabel('phi'); plt.ylabel('XN2');
    plt.figure(7, figsize = (8,4))
    plt.plot(phi,XH2O)
    plt.xlabel('phi'); plt.ylabel('XH2O');


#save plot FT
dpi = 1000
plt.figure(1); plt.title('r = %s' %str(r)); plt.legend(['Ti=25ºC','Ti=35ºC','Ti=45ºC','Ti=55ºC'], loc=0); plt.savefig("P_%s.png" % leg , dpi=dpi)
plt.figure(2); plt.title('r = %s' %str(r)); plt.legend(['Ti=25ºC','Ti=35ºC','Ti=45ºC','Ti=55ºC'], loc=0); plt.savefig("T_%s.png" % leg , dpi=dpi)
plt.figure(3); plt.title('r = %s' %str(r)); plt.legend(['Ti=25ºC','Ti=35ºC','Ti=45ºC','Ti=55ºC'], loc=0); plt.savefig("XCH3OH_%s.png" % leg , dpi=dpi)
plt.figure(4); plt.title('r = %s' %str(r)); plt.legend(['Ti=25ºC','Ti=35ºC','Ti=45ºC','Ti=55ºC'], loc=0); plt.savefig("XCO2_%s.png" % leg , dpi=dpi)
plt.figure(5); plt.title('r = %s' %str(r)); plt.legend(['Ti=25ºC','Ti=35ºC','Ti=45ºC','Ti=55ºC'], loc=0); plt.savefig("XCO_%s.png" % leg , dpi=dpi)
plt.figure(6); plt.title('r = %s' %str(r)); plt.legend(['Ti=25ºC','Ti=35ºC','Ti=45ºC','Ti=55ºC'], loc=0); plt.savefig("XH2O_%s.png" % leg , dpi=dpi)
plt.figure(7); plt.title('r = %s' %str(r)); plt.legend(['Ti=25ºC','Ti=35ºC','Ti=45ºC','Ti=55ºC'], loc=0); plt.savefig("XN2_%s.png" % leg , dpi=dpi)
"""
"""
#save plot FP
dpi = 1000
plt.figure(1); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("P_%s.png" % leg , dpi=dpi)
plt.figure(2); plt.title('r = %s' %str(r)); plt.legend(['Pi=1bar','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("T_%s.png" % leg , dpi=dpi)
plt.figure(3); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("XCH3OH_%s.png" % leg , dpi=dpi)
plt.figure(4); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("XCO2_%s.png" % leg , dpi=dpi)
plt.figure(5); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("XCO_%s.png" % leg , dpi=dpi)
plt.figure(6); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("XH2O_%s.png" % leg , dpi=dpi)
plt.figure(7); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("XN2_%s.png" % leg , dpi=dpi)
"""
