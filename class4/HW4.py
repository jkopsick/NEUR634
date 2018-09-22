#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 22 08:03:47 2018

@author: jeffk
"""

# Import necessary libraries to perform the experiment and plot results
import moose
import pylab as plt
import numpy as np
import util as u

# Define the time that the simulation should run for, as well as the timestep
simTime = 1000e-3
simDt = 5e-5

# Define the variables needed for the creation of the compartments and pulse
RM = 2
CM = 0.01
RA = 2
Em = -65e-3
pulse_dur = 100e-3
pulse_amp = 0.5e-9
pulse_delay1 = 50e-3
pulse_delay2 = 500e-3

# Load the model files into moose
pfile = 'layer2.p'
swcfile = '33-RVT162-1-4e.CNG.swc'
container1 = 'l2_cell'
container2 = 'gran_cell'
l2_cell = moose.loadModel(pfile,container1)
gran_cell = moose.loadModel(swcfile,container2)                         

# Compute specific resistivity, capacitance, and axial resistance for swc
# (p files already have this computed upon load if specified at top of pfile)
u.setCompParameters(gran_cell,'Compartment',RM,CM,RA,Em)

l2_cell = moose.wildcardFind('/l2_cell/#[TYPE=SymCompartment]')
gran_cell = moose.wildcardFind('/gran_cell/#[TYPE=Compartment]')

# Create data hierarchies for the two different cells
l2_data = moose.Neutral('/l2_data')
gran_data = moose.Neutral('/gran_data')

# Create pulses and connect them to each neuron's soma
l2_soma_pulse = u.createPulse(l2_cell[0], 'rollingWave', pulse_dur, pulse_amp, 
                           pulse_delay1, pulse_delay2)
g_soma_pulse = u.createPulse(gran_cell[0], 'rollingWave', pulse_dur, pulse_amp, 
                           pulse_delay1, pulse_delay2)

# Create data tables for the membrane potential for each compartment
l2_tables = []
gran_tables = []
for comp in l2_cell:
    l2_tables.append(u.createDataTables(comp,l2_data))
for comp in gran_cell:
    gran_tables.append(u.createDataTables(comp,gran_data))
    
# Choose the soma and a representative dendrite to view how voltage changes
# for the layer 2 cell
l2_soma_Vm = l2_tables[0]
l2_ap2_Vm = l2_tables[2]

# Choose the soma and a representative dendrite to view how voltage changes
# for the granule cell
g_soma0_Vm = gran_tables[0]
g_dend_5_00_Vm = gran_tables[-11]

# Run the experiment and plot the results for both neurons
moose.setClock(4,simDt)
moose.reinit()
moose.start(simTime)

f1 = plt.figure(1)
t = plt.linspace(0,simTime,len(l2_soma_Vm.vector))
plt.plot(t,l2_soma_Vm.vector, 'r',label = 'l2_soma_Vm')
plt.plot(t,l2_ap2_Vm.vector,'b',label = 'l2_dend_Vm')
plt.legend()
f1.show()

f2 = plt.figure(2)
t = plt.linspace(0,simTime,len(g_soma0_Vm.vector))
plt.plot(t,g_soma0_Vm.vector, 'r',label = 'g_soma_Vm')
plt.plot(t,g_dend_5_00_Vm.vector,'b',label = 'g_dend_Vm')
plt.legend()
f2.show()