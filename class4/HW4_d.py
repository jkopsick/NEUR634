#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 22 16:39:19 2018

@author: jeffk
"""

# Import NEURON into the python terminal
from neuron import h,gui
import util as u
import pylab as plt

# Set values for the pulse and for properties of the compartments
pulse_dur = 15
pulse_amp = 0.5
pulse_delay = 25
sim_dur = 100
soma_L = 12.6157
soma_D = 12.6157
dend_L = 200
dend_D = 1
dend_nseg = 101
soma_Ra = 2
soma_gpas = 1e-3
dend_Ra = 2
dend_gpas = 1e-4

# Create a soma and dendrite
soma = h.Section(name='soma')
dend = h.Section(name='dend')

# Divide the dendrite up into the specifed number of segments
dend.nseg = dend_nseg

# Connect the soma to the dendrite
dend.connect(soma(1))

# Create a passive cell membrane and set parameters for each compartment
soma.insert('pas')
dend.insert('pas')
soma.Ra = soma_Ra
soma(0.5).g_pas = soma_gpas
dend.Ra = dend_Ra
dend(0.5).g_pas = dend_gpas

# Define the length and diameter for the soma and dendrite
soma.L = soma.diam = soma_L
dend.L = dend_L
dend.diam = dend_D

# Store the middle segment of each compartment to use for the pulse and for
# measurement in the experiment

mid_soma = soma(0.5)
mid_dend = soma(0.5)

# Inject current into the soma
shock = u.createNeuronPulse(mid_soma,'shock',pulse_dur,pulse_amp,pulse_delay)

# Define where to record from for each compartment
soma_vec = u.recordCompNeuron(mid_soma)
dend_vec = u.recordCompNeuron(mid_dend)
#soma_v_vec = h.Vector()
#soma_v_vec.record(mid_soma._ref_v)
#dend_v_vec = h.Vector()
#dend_v_vec.record(mid_dend._ref_v)
t_vec = u.recordTimeNeuron() 

# Run the simulation
h.tstop = sim_dur
h.run()

# Plot the simulation
plt.plot(t_vec, soma_vec, 'r')
plt.plot(t_vec, dend_vec, 'b')
plt.xlabel('time (ms)')
plt.ylabel('V (mV)')
plt.show()
