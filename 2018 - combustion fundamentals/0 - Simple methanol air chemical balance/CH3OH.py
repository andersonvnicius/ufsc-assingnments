# -*- coding: utf-8 -*-
"""
--
Calculating chemical equilibrium of the combustion of methanol(CH3OH) with stardard athmospheric air
--
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

#criando vetores para armazenar dados
T = np.ones(n_of_steps)
P = np.ones(n_of_steps)
EG = np.ones(n_of_steps)
rho = np.ones(n_of_steps)
phi = np.ones(n_of_steps)
afr = np.ones(n_of_steps)
XH2O = np.ones(n_of_steps)
XCO2 = np.ones(n_of_steps)
XCO = np.ones(n_of_steps)
XN2 = np.ones(n_of_steps)

#definindo numero de passos de phi
n_of_steps = 100

#dados de entrada
r = 14                          	#razao de compressao do motor
intake_temperature = 273. + 25.     #temperatura ambiente
intake_pressure = 101325.       	#pressao na admissao

#criando objeto com propriedades dos reagentes
meth = ct.Solution('gri30.xml')
meth.X = {'CH3OH':1,'O2':1.5 * 1,'N2':1.5 *3.76}

#calculando constante k
k = meth.cp/meth.cv

#calculando pressao e temperatura na combustao
spark_temperature = intake_temperature*r**(k - 1)
spark_pressure = intake_pressure*r**(k)

#criando legendas e file outputs
legPc = 'r_' +str(r) + '_FP1atm'; print(r)
legTc = 'r_' +str(r) + '_FT25C'; print(r)

for j in range (0,4):
    
	#calculando nova iteracao de temperatura ou pressao
	intake_temperature = intake_temperature + 10*j
	intake_pressure = intake_pressure #+ 15198.75E3*j
	
	#calculando nova temp e pressao na combustao
    spark_temperature = intake_temperature*r^(k - 1)
    spark_pressure = intake_pressure*r**(k-1)
	
    for i in range (0,n_of_steps):
	
        #calculando novo valor de phi e resetando composicao quimica
		phi[i] = 0.6 + i * 1/n_of_steps
        meth.X = np.zeros(53)
		
		#definindo nova condicao quimica
        meth.X = {'CH3OH':1,'O2':1.5/phi[i] * 1,'N2':1.5/phi[i] *3.76}
        meth.TP = spark_temperature , spark_pressure
        meth.equilibrate('UV')
        
		#output variables
        T[i] = meth.T
        P[i] = meth.P /1E6
        XH2O[i] = meth['H2O'].X
        XCO2[i] = meth['CO2'].X
        XCO[i] = meth['CO'].X
        XN2[i] = meth['N2'].X
    
    #max values
    print("FT=", fuel_temperature - 273)
    #print("FP=", fuel_pressure)
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


#save plot FP
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
#save plot FT
dpi = 1000
plt.figure(1); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("P_%s.png" % leg , dpi=dpi)
plt.figure(2); plt.title('r = %s' %str(r)); plt.legend(['Pi=1bar','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("T_%s.png" % leg , dpi=dpi)
plt.figure(3); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("XCH3OH_%s.png" % leg , dpi=dpi)
plt.figure(4); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("XCO2_%s.png" % leg , dpi=dpi)
plt.figure(5); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("XCO_%s.png" % leg , dpi=dpi)
plt.figure(6); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("XH2O_%s.png" % leg , dpi=dpi)
plt.figure(7); plt.title('r = %s' %str(r)); plt.legend(['Pi=1atm','Pi=1.15atm','Pi=1.3atm','Pi=1.45atm'], loc=0); plt.savefig("XN2_%s.png" % leg , dpi=dpi)
"""
