# Import necessary libraries to perform the experiment and plot results
import moose
import pylab as plt
import numpy as np
import util as u
from chan_proto_part1 import chan_set, rateParams, Ca_pool_params, CaParams

# Define the variables needed for the creation of the compartments and pulse
cond_set = {'Na': 10000, 'K': 2500, 'KaF': 0, 'SKCa': 0, 'CaL': 0}
pulse_dur = 100e-3
pulse_amp = 0.75e-9
pulse_delay1 = 20e-3
pulse_delay2 = 1e9

# Load a multi-compartment model into MOOSE and give it channels
pfile = 'layer2.p'
container = 'l2_cell'
libraryName = '/library'
compType = 'SymCompartment'
l2_cell = u.createMultiCompCell(file_name = pfile, container_name = container, library_name = libraryName,
				comp_type = compType, channelSet = chan_set, condSet = cond_set, 
				rateParams = rateParams, CaParams = CaParams, 
				CaPoolParams = Ca_pool_params)

# Re-create the python variable pointing to the gran_cell to limit results just to 
# type compartment (excludes the spines and allows for createDataTables function to
# work properly)
l2_cell = moose.wildcardFind(l2_cell.path + '/' + '#[TYPE=' + compType + ']')

# Create the pulse and apply it to the granule cell's soma
l2_soma_pulse = u.createPulse(l2_cell[0], 'rollingWave', pulse_dur, pulse_amp, 
                           pulse_delay1, pulse_delay2)

# Create a neutral object to store the data in
l2_data = moose.Neutral('/l2_data')
l2_tables = []
for comp in l2_cell:
    l2_tables.append(u.createDataTables(comp,l2_data,l2_soma_pulse))

# Choose the soma and a representative dendrite to view how voltage changes
# for the granule cell
l2_soma_Vm = l2_tables[0][0]
l2_soma_Iex = l2_tables[0][1]
l2_ap2_Vm = l2_tables[2][0]
l2_ap2_Iex = l2_tables[2][1]

# Plot the simulation
simTime = 0.5
simdt = 0.25e-5
plotdt = 0.25e-3
for i in range(10):
    moose.setClock(i, simdt)

moose.setClock(8, plotdt)
moose.reinit()
moose.start(simTime)
t = np.linspace(0,simTime,len(l2_soma_Vm.vector))
plt.plot(t,l2_soma_Vm.vector * 1e3, 'r',label = 'l2_soma_Vm (mV)')
plt.plot(t,l2_ap2_Vm.vector * 1e3,'b',label = 'l2_ap2_Vm (mV)')
plt.plot(t, l2_soma_Iex.vector * 1e3, label='l2_soma_Iex (nA)')
plt.plot(t, l2_ap2_Iex.vector * 1e9, label='l2_ap2_Iex (nA)')
plt.xlabel('time (ms)')
plt.legend()
plt.show()
