#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 15 07:51:45 2018

@author: jeffk
"""

# Import necessary libraries to perform the experiment and plot results
import moose
from pylab import *
import numpy as np
import os
owd = os.getcwd()
os.chdir("..")
from util import *
os.chdir(owd)

# Main code
def main(args=None):
	l = createCompartment('swagginWagon',30e-6,10e-6,2,0.01,-65e-3)
	n = createPulse(l, 'rollingWave', 100e-3, 0.1e-9, 50e-3)
	out = moose.Neutral('/out')
	vmtab = moose.Table('/out/somaVm')
	moose.connect(vmtab, 'requestOut', l, 'getVm')
	# Set the current injection to zero to see sole effect of pulse
	l.inject = 0
	# Set the clock so that all components of the simulation run 		  appropriately
	#l.tick = 2
	#moose.setClock(2, 0.25e-3)
	# Run the experiment
	moose.reinit()
	moose.start(0.3)
	# Compute delta VM at the steady state
	delta_Vm = l.Rm*n.level[0]*(1-np.exp(-0.3/(l.Rm*l.Cm)))
	print delta_Vm
	# Compute actual delta Vm at the steady state
	actual_delta_Vm = np.max(vmtab.vector) - np.min(vmtab.vector)
	print actual_delta_Vm
	# Plot the Voltage for the compartment for length of the 		  experiment
	t = linspace(0, 300e-3, len(vmtab.vector))
	plot(t, vmtab.vector)
	show()

if __name__ == '__main__':
	main()
