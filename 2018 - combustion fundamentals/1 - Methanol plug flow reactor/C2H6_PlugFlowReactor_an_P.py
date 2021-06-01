# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 16:04:13 2018

@author: andersonvnicius
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

#set simulation parameters
n_of_steps = 1000
n_of_samples = 4
analysis = "P"

#creating graphs
plt.figure(1)
plt.xlabel('$z$ [m]')
plt.ylabel('$T$ [K]')
plt.figure(2)
plt.xlabel('$z$ [m]')
plt.ylabel('$Y_{C2H6}$')

for i  in range (0,n_of_samples):
    #fuel input data
    init_T = 1000.
    init_p = (0.2 +i*0.1) * ct.one_atm
    phi = 0.8
    a = 16/phi
    
    #reactor input data
    lenght = 10E-2
    u0 = 4.25E-2
    r = 3E-2
    A = np.pi*r**2
    
    #creating fuel
    eth = ct.Solution('gri30.xml')
    eth.TPX = init_T, init_p, {'C2H6':1, 'O2':a * 1, 'N2':a * 3.76}
    
    #creating reactor
    reactor = ct.IdealGasConstPressureReactor(eth)
    sim = ct.ReactorNet([reactor])
    
    #calculating fuel mass flow rate
    mass_flow_rate = u0*A*eth.density
    
    #approximate a time step tp achieve a similar resolution
    t_total = lenght/u0
    dt = t_total/n_of_steps
    
    #creating data output vectors
    t = (np.arange(n_of_steps) + 1)*dt
    z = np.zeros_like(t)
    u = np.zeros_like(t)
    states = ct.SolutionArray(reactor.thermo)
    
    #solving problem
    for n1, t_i in enumerate(t):
        #time integration
        sim.advance(t_i)
        #compute velocity
        u[n1] = mass_flow_rate/A/reactor.thermo.density
        #transform velocity into space
        z[n1] = z[n1 - 1] + u[n1]*dt
        states.append(reactor.thermo.state)
    
    #graph plot
    init_p = round(init_p, 4)
    analysis = 'P = %s' %init_p
    plt.figure(1)
    plt.plot(z, states.T, label=analysis)
    plt.figure(2)
    plt.plot(z, states.Y[:, eth.species_index('C2H6')], label=analysis)
    
#plot graphs
plt.figure(1)
plt.legend(loc=0); #plt.savefig('PFR_T_z_%s.png' %analysis)

plt.figure(2)
plt.legend(loc=0); #plt.savefig('PFR_XC2H6_t_%s.png' %analysis)

"""
#plot graphs
plt.figure()
plt.plot(z, states.T, label='temperature')
plt.xlabel('$z$ [m]')
plt.ylabel('$T$ [K]')
plt.legend(loc=0)
plt.show()
plt.savefig('pfr_T_z.png')

plt.figure()
plt.plot(z, states.X[:, eth.species_index('C2H6')], label='remanescent C2H6')
plt.plot(z, states.X[:, eth.species_index('CO2')], label='remanescent CO2')
plt.plot(z, states.X[:, eth.species_index('CO')], label='remanescent CO')
plt.plot(z, states.X[:, eth.species_index('H2O')], label='remanescent H2O')
plt.xlabel('$z$ [m]')
plt.ylabel('$X_{element}$')
plt.legend(loc=0)
plt.show()
plt.savefig('pfr_Xremanescent_t.png')
"""
