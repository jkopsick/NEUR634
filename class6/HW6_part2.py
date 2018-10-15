# Import necessary libraries to perform the experiment and plot results
import moose
import pylab as plt
import numpy as np
import util as u
from chan_proto_part1 import chan_set, rateParams

# Define the variables needed for the creation of the compartments and pulse
EREST_ACT = -70e-3 #: Resting membrane potential
RM = 1/(0.3e-3*1e4)
CM = 1e-6*1e4
RA = 1
Em = EREST_ACT + 10.613e-3
initVm = EREST_ACT
cond_set = {'Na': 10000, 'K': 2500}
pulse_dur = 100e-3
pulse_amp = 0.75e-9
pulse_delay1 = 20e-3
pulse_delay2 = 1e9

# Load a multi-compartment model into MOOSE and give it channels
swcfile = '33-RVT162-1-4e.CNG.swc'
container = 'gran_cell'
libraryName = '/library'
compType = 'Compartment'
gran_cell = u.createMultiCompCell(swcfile, container, libraryName, compType, chan_set, cond_set,
                                  rateParams, RM, CM, RA, initVm, Em)
moose.showfield('/gran_cell/soma')

# Re-create the python variable pointing to the gran_cell to limit results just to 
# type compartment (excludes the spines and allows for createDataTables function to
# work properly)
gran_cell = moose.wildcardFind('/gran_cell/#[TYPE=Compartment]')

# Create the pulse and apply it to the granule cell's soma
g_soma_pulse = u.createPulse(gran_cell[0], 'rollingWave', pulse_dur, pulse_amp, 
                           pulse_delay1, pulse_delay2)

# Create a neutral object to store the data in
gran_data = moose.Neutral('/gran_data')
gran_tables = []
for comp in gran_cell:
    gran_tables.append(u.createDataTables(comp,gran_data,g_soma_pulse))

# Choose the soma and a representative dendrite to view how voltage changes
# for the granule cell
g_soma0_Vm = gran_tables[0][0]
g_soma0_Iex = gran_tables[0][1]
g_dend_5_00_Vm = gran_tables[-11][0]
g_dend_5_00_Iex = gran_tables[-11][1]


# Plot the simulation
simTime = 0.5
simdt = 0.25e-5
plotdt = 0.25e-3
for i in range(10):
    moose.setClock(i, simdt)

moose.setClock(8, plotdt)
moose.reinit()
moose.start(simTime)
t = np.linspace(0,simTime,len(g_soma0_Vm.vector))
plt.plot(t,g_soma0_Vm.vector * 1e3, 'r',label = 'g_soma_Vm (mV)')
plt.plot(t,g_dend_5_00_Vm.vector * 1e3,'b',label = 'g_dend_Vm (mV)')
plt.plot(t, g_soma0_Iex.vector * 1e3, label='g_soma_Iex (nA)')
plt.plot(t, g_dend_5_00_Iex.vector * 1e9, label='g_dend_Iex (nA)')
plt.xlabel('time (ms)')
plt.legend()
plt.show()
