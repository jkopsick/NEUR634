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
RM = 0.8
CM = 0.01
RA = 0.9
Em = -0.072
pulse_dur = 10e-3
pulse_amp = 0.1e-9
pulse_delay1 = 10e-3
pulse_delay2 = 2500e-3

# Load the model files into moose
swcfile = '72-1.CNG.swc'
container = 'cell1'
sd_cell = moose.loadModel(swcfile,container)                         

# Compute specific resistivity, capacitance, and axial resistance for swc
# (p files already have this computed upon load if specified at top of pfile)
u.setCompParameters(sd_cell,'Compartment',RM,CM,RA,Em)

# Re-define the cell variables so that they just include compartments
# (i.e. exclude spines)
sd_cell = moose.wildcardFind('/cell1/#[TYPE=Compartment]')

# Create data hierarchies for the two different cells
sd_data = moose.Neutral('/gran_data')


sd_soma_pulse = u.createPulse(sd_cell[0], 'rollingWave', pulse_dur, pulse_amp, 
                           pulse_delay1, pulse_delay2)

# Create data tables for the membrane potential for each compartment
sd_tables = []
for comp in sd_cell:
    sd_tables.append(u.createDataTables(comp,sd_data))

# Choose the soma and a representative dendrite to view how voltage changes
# for the granule cell
sd_soma0_Vm = sd_tables[0]

# Run the experiment and plot the results for both neurons
moose.setClock(4,simDt)
moose.reinit()
moose.start(simTime)

f2 = plt.figure(2)
t = plt.linspace(0,simTime,len(sd_soma0_Vm.vector))
plt.plot(t,sd_soma0_Vm.vector, 'r',label = 'sd_soma_Vm')
plt.legend()
f2.show()
