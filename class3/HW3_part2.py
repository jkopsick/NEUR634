#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
HW #3 Part 2

Created on Sun Sep 16 13:41:14 2018

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

# Define variables needed for the creation of the compartments and pulse
RM = 0.5
CM = 0.01
RA = 8
soma_l = 30e-6
soma_rad = 10e-6
dend_l = 100e-6
dend_rad = 1e-6
Em = -65e-3
pulse_dur = 100e-3
pulse_amp = 0.1e-9
pulse_delay = 50e-3

# Create a neutral directory for the neuron and create the compartments
neuron = moose.Neutral('/neuron')
soma = u.createCompartment(neuron,'swagginSoma',soma_l,soma_rad,RM,CM,RA,Em)
dend = u.createCompartment(neuron,'swagginDend',dend_l,dend_rad,RM,CM,RA,Em)

# Create the pulse
soma_pulse = u.createPulse(soma, 'rollingWave', pulse_dur, pulse_amp, pulse_delay)

# Connect the soma to the dendrite
moose.connect(soma,'axialOut',dend,'handleAxial')

# Create a neutral directory for the data and create data tables for the soma
# and the dendrite
data = moose.Neutral('/data')
soma_Vm = u.createDataTables(soma,data)
dend_Vm = u.createDataTables(dend,data)

# Set the clock so that all components of the simulation run           appropriately
moose.setClock(4, simDt)

# Run the experiment
moose.reinit()
moose.start(simTime)

# Plot the voltages for both the soma and dendrite
t = plt.linspace(0, simTime, len(soma_Vm.vector))
plt.plot(t, soma_Vm.vector,'r', label='Vm_soma')
plt.plot(t,dend_Vm.vector,'b',label='Vm_dend')
plt.legend()
plt.show()
