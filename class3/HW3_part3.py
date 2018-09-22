#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
HW #3 Part 3

Created on Mon Sep 17 09:37:34 2018

@author: jeffk
"""

# Import necessary libraries to perform the experiment and plot results
import moose
import pylab as plt
import numpy as np
import util as u

# Define the time that the simulation should run for, as well as the timestep
simTime = 300e-3
simDt = 5e-5

# Define the variables needed for the creation of the compartments and pulse
numDends = 9
soma_l = 30e-6
soma_rad = 10e-6
dend_l = 100e-6
dend_rad = 1e-6
RM = 2
CM = 0.01
RA = 10
Em = -65e-3
pulse_dur = 100e-3
pulse_amp = 0.1e-9
pulse_delay = 50e-3

# Create a neutral directory for the neuron
neuron = moose.Neutral('/neuron')

# Create the soma under the neuron directory
soma = u.createCompartment(neuron,'swagginSoma',soma_l,soma_rad,RM,CM,RA,Em)

# Create the dendritic compartments under the neuron directory
dend_branch = u.discretize(neuron,numDends,dend_l,dend_rad,RM,CM,RA,Em)

# Connect the soma to the first element of the dendritic branch
moose.connect(soma,'axialOut',dend_branch[0],'handleAxial')

# Create a pulse and connect it to the soma
soma_pulse = u.createPulse(soma, 'rollingWave', pulse_dur, pulse_amp, pulse_delay)

# Create data tables to store the voltage for the soma and each compartment
# making up the dendritic branch
data = moose.Neutral('/data')
soma_Vm = u.createDataTables(soma,data)
dend_tables = []
for dend in dend_branch:
    dend_tables.append(u.createDataTables(dend,data))

# Store the location of the middle dendrite so that it can be accessed
# for plotting purposes
mid_dend = len(dend_branch)//2
    
# Run the experiment and plot the results
moose.setClock(4,simDt)
moose.reinit()
moose.start(simTime)
t = plt.linspace(0,simTime,len(soma_Vm.vector))
plt.plot(t,soma_Vm.vector, 'r',label = 'Vm_soma')
plt.plot(t,dend_tables[mid_dend].vector,'b',label = 'Vm_mid_dend')
plt.legend()
plt.show()
