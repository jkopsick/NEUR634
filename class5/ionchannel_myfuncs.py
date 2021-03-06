#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 13:27:01 2018

@author: jeffk

This python file adapts ionchannel.py and uses the functions that I have
created.
"""

# ionchannel.py --- 
# 
# Filename: ionchannel.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Wed Sep 17 10:33:20 2014 (+0530)
# Version: 
# Last-Updated: 
#           By: 
#     Update #: 0
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# 
# 
# 

# Change log:
# 
# 
# 
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth
# Floor, Boston, MA 02110-1301, USA.
# 
# 

# Code:
"""This demo shows how to set the parameters for a Hodgkin-Huxley type ion channel.

Hodgkin-Huxley type ion channels are composed of one or more gates
that allow ions to cross the membrane. The gates transition between
open and closed states and this, taken over a large population of
ion channels over a patch of membrane has first order kinetics, where
the rate of change of fraction of open gates (n) is given by::

    dn/dt = alpha(Vm) * (1 - n) - beta(Vm) * n

where alpha and beta are rate parameters for gate opening and
closing respectively that depend on the membrane potential.
The final channel conductance is computed as::

    Gbar * m^x * h^y ...

where m, n are the fraction of open gates of different types and x,
y are the number of such gates in each channel. We can define the
channel by specifying the alpha and beta parameters as functions of
membrane potential and the exponents for each gate.
The number gates is commonly one or two.

Gate opening/closing rates have the form::

    y(x) = (A + B * x) / (C + exp((x + D) / F))

where x is membrane voltage and y is the rate parameter for gate
closing or opening.


"""

# Import all libraries necessary to run the file
import numpy as np
import matplotlib.pyplot as plt
import moose
import util as u

# Set the resting membrane potential to be used in creating the channels
# and defining their biophysics
EREST_ACT = -70e-3 #: Resting membrane potential


#: The parameters for defining m as a function of Vm
Na_m_params = [1e5 * (25e-3 + EREST_ACT),   # 'A_A':
                -1e5,                       # 'A_B':
                -1.0,                       # 'A_C':
                -25e-3 - EREST_ACT,         # 'A_D':
               -10e-3,                      # 'A_F':
               4e3,                     # 'B_A':
                0.0,                        # 'B_B':
                0.0,                        # 'B_C':
                0.0 - EREST_ACT,            # 'B_D':
                18e-3                       # 'B_F':    
               ]

#: Parameters for defining h gate of Na+ channel
Na_h_params = [ 70.0,                        # 'A_A':
                0.0,                       # 'A_B':
                0.0,                       # 'A_C':
                0.0 - EREST_ACT,           # 'A_D':
                0.02,                     # 'A_F':
                1000.0,                       # 'B_A':
                0.0,                       # 'B_B':
                1.0,                       # 'B_C':
                -30e-3 - EREST_ACT,        # 'B_D':
                -0.01                    # 'B_F':       
                ]

#: K+ channel in Hodgkin-Huxley model has only one gate, n and these
#are the parameters for the same
K_n_params = [ 1e4 * (10e-3 + EREST_ACT),   #  'A_A':
               -1e4,                      #  'A_B':
               -1.0,                       #  'A_C':
               -10e-3 - EREST_ACT,         #  'A_D':
               -10e-3,                     #  'A_F':
               0.125e3,                   #  'B_A':
               0.0,                        #  'B_B':
               0.0,                        #  'B_C':
               0.0 - EREST_ACT,            #  'B_D':
               80e-3                       #  'B_F':  
               ]

#: We define the rate parameters, which are functions of Vm as
#: interpolation tables looked up by membrane potential.
#: Minimum x-value for the interpolation table
VMIN = -30e-3 + EREST_ACT
#: Maximum x-value for the interpolation table
VMAX = 120e-3 + EREST_ACT
#: Number of divisions in the interpolation table
VDIVS = 3000


def create_na_proto():
    """Create and return a Na+ channel prototype '/library/na'

    The Na+ channel conductance has the equation::

        g = Gbar * m^3 * h 

    For each gate, we use the HHChannel.setupAlpha function to set up
    the interpolation table.

    """
    lib = moose.Neutral('/library')
    na = moose.HHChannel('/library/na')
    na.tick = -1
    #: The exponent for m gate is 3
    na.Xpower = 3
    #: channel/gateX is the m gate
    #: setting Xpower to a positive number automatically creates this gate.
    xGate = moose.element(na.path + '/gateX')    
    xGate.setupAlpha(Na_m_params +
                      [VDIVS, VMIN, VMAX])
    #: channel/gateY is the h gate
    #: Exponent for h gate is 1
    na.Ypower = 1
    yGate = moose.element(na.path + '/gateY')
    yGate.setupAlpha(Na_h_params + 
                      [VDIVS, VMIN, VMAX])
    return na
    


def create_k_proto():
    """Create and return a K+ channel prototype '/library/k'.

    The K+ channel conductance has the equation::

        g = Gbar * n^4

    """
    lib = moose.Neutral('/library')
    k = moose.HHChannel('/library/k')
    k.tick = -1
    k.Xpower = 4.0
    xGate = moose.HHGate(k.path + '/gateX')    
    xGate.setupAlpha(K_n_params +
                      [VDIVS, VMIN, VMAX])
    return k
    

if __name__ == '__main__':
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
    
    # Create copies of the prototype channels and connect them to the
    # compartment
    nachan = moose.copy(create_na_proto(), container, 'na', 1)
    nachan.Gbar = 120e-3 * SA * 1e4 # Gbar_Na = 120 mS/cm^2
    nachan.Ek = 115e-3 + EREST_ACT
    moose.connect(nachan, 'channel', squid, 'channel', 'OneToOne')
    kchan = moose.copy(create_k_proto(), container, 'k', 1)
    kchan.Gbar = 36e-3 * SA * 1e4 # Gbar_K = 36 mS/cm^2
    kchan.Ek = -12e-3 + EREST_ACT
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
