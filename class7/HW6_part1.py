#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 13:27:01 2018

@author: jeffk

This python file adapts ionchannel.py and uses the functions that I have
created.
"""

# Import all libraries necessary to run the file
import numpy as np
import matplotlib.pyplot as plt
import moose
import util as u
from chan_proto_part1 import chan_set, rateParams, CaParams

# Set the resting membrane potential to be used in creating the channels
# and defining their biophysics
EREST_ACT = -70e-3 #: Resting membrane potential  

# Define the variables needed to create the compartment and the channels
squid_l = 50e-6
squid_rad = 15e-6
SA = np.pi*squid_l*2*squid_rad
Em = EREST_ACT + 10.613e-3
initVm = EREST_ACT
RM = 1/(0.3e-3*1e4)
CM = 1e-6*1e4
RA = 1
pulse_dur = 40e-3
pulse_delay1 = 20e-3
pulse_amp = 1e-9
pulse_delay2 = 1e9

# Create the compartment
neuron = moose.Neutral('/neuron')
squid = u.createCompartment(neuron,'HHsquidswag',squid_l,squid_rad,RM,CM,
                            RA,initVm,Em)

# Set the axial resistance to match the one in the ionchannel.py file
# and set the container to be used for channel creation
squid.Ra = 1
container = squid.parent.path

# Create the external current to be applied to the compartment
squid_pulse = u.createPulse(squid, 'rollingWave', pulse_dur, pulse_amp,
                           pulse_delay1, pulse_delay2)

# Create the data tables necessary to run the model
data = moose.Neutral('/data')
squid_Vm, squid_current = u.createDataTables(squid,data,squid_pulse)

# Create the channels necessary for the giant squid axon and place them in a library
u.createChanLib('/library', chan_set, rateParams, CaParams)

# Create copies of the prototype channels and connect them to the
# compartment
nachan = moose.copy('/library/Na', container, 'na', 1)
nachan.Gbar = 120e-3 * SA * 1e4 # Gbar_Na = 120 mS/cm^2
moose.connect(nachan, 'channel', squid, 'channel', 'OneToOne')

kchan = moose.copy('/library/K', container, 'k', 1)
kchan.Gbar = 36e-3 * SA * 1e4 # Gbar_K = 36 mS/cm^2
moose.connect(kchan, 'channel', squid, 'channel', 'OneToOne')

# Plot the simulation
simtime = 0.1
simdt = 0.25e-5
plotdt = 0.25e-3
for i in range(10):
    moose.setClock(i, simdt)
moose.setClock(8, plotdt)
moose.reinit()
moose.start(simtime)
t = np.linspace(0, simtime, len(squid_Vm.vector))
plt.plot(t, squid_Vm.vector * 1e3, label='Vm (mV)')
plt.plot(t, squid_current.vector * 1e9, label='current (nA)')
plt.legend()
plt.show()
