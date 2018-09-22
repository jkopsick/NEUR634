#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
HW #3 Part 1

Created on Sat Sep 15 08:58:52 2018

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
soma_l = 30e-6
soma_rad = 10e-6
RM = 1
CM = 0.01
RA = 1
Em = -65e-3
pulse_dur = 100e-3
pulse_amp = 0.1e-9
pulse_delay1 = 50e-3
pulse_delay2 = 1e9

# Create a neutral directory for the neuron
neuron = moose.Neutral('/neuron')

# Create the soma and connect a pulse to it
l = u.createCompartment(neuron,'swagginWagon',soma_l,soma_rad,RM,CM,RA,Em)
n = u.createPulse(l,'rollingWave',pulse_dur,pulse_amp,pulse_delay1,pulse_delay2)

# Create a neutral directory for the data and create data tables for the soma
data = moose.Neutral('/data')
m = u.createDataTables(l,data)

# Set the clock so that all components of the simulation run           appropriately
moose.setClock(4, simDt)

# Run the experiment
moose.reinit()
moose.start(simTime)

# Compute delta VM at the steady state
delta_Vm = l.Rm*n.level[0]*(1-np.exp(-0.3/(l.Rm*l.Cm)))
print delta_Vm

# Compute actual delta Vm at the steady state
actual_delta_Vm = (np.max(m.vector) -
                   np.min(m.vector))
print actual_delta_Vm

# Plot the voltage for the soma, along with a line for tau
t = plt.linspace(0, simTime, len(m.vector))
plt.plot(t, m.vector)
plt.axvline(x=50e-3 + l.Rm*l.Cm)
plt.show()
