# Import necessary libraries to perform the experiment and plot results
import moose
import pylab as plt
import numpy as np
import util as u
from chan_proto_part1 import chan_set, rateParams

# Define the time that the simulation should run for, as well as the timestep
simTime = 1000e-3
simDt = 5e-5

# Define the variables needed for the creation of the compartments and pulse
RM = 2
CM = 0.01
RA = 2
Em = -65e-3
cond_set = {'Na': 120, 'K': 36}
pulse_dur = 100e-3
pulse_amp = 0.5e-9
pulse_delay1 = 50e-3
pulse_delay2 = 500e-3

# Load a multi-compartment model into MOOSE and give it channels
swcfile = '33-RVT162-1-4e.CNG.swc'
container = 'gran_cell'
libraryName = '/library'
compType = 'SymCompartment'
gran_cell = u.createMultiCompCell(swcfile, container, libraryName, compType, chan_set, cond_set,
                                  rateParams, RM, CM, RA, Em)
