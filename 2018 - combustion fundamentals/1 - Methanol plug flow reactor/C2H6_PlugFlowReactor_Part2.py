# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 16:04:13 2018

@author: andersonvnicius
"""

import cantera as ct
import numpy as np

#set simulation parameters
n_of_steps = 10000
i=0
R=1

while R>0.01:
    #start
    i=i+1
    print('Iteration=',i)
    
    #fuel input data
    init_T = 1000.
    init_p = 0.2 * ct.one_atm
    phi = 1
    a = 16/phi
    
    #reactor input data
    lenght = 10E-3
    u0 = (4.18 - 0.00025*i)*1E-3 
    r = 5E-2
    A = np.pi*r**2
    print('Initial fluid velocity=',u0)
    
    #creating fuel
    eth = ct.Solution('gri30.xml')
    eth.TPX = init_T, init_p, {'C2H6':1,'O2':a * 1,'N2':a * 3.76}
    Yi = np.min(eth['C2H6'].Y)
    
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
        #compute velocity and transform into space
        u[n1] = mass_flow_rate/A/reactor.thermo.density
        z[n1] = z[n1 - 1] + u[n1]*dt
        states.append(reactor.thermo.state)
    
    Yf = np.min(states.Y[:, eth.species_index('C2H6')])
    R = Yf/Yi
    print('Remanescent fuel Y=',R)