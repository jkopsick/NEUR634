# Load in the libraries and functions needed to create the model
import moose
import numpy as np
import pylab as plt
import util as u

# Set plt to show the graphs that will be created for each channel type
plt.ion()

# Read the NML model into MOOSE
filename = 'L23_NoHotSpot.cell.nml'
reader = moose.mooseReadNML2(filename)

# Define the variables needed for the creation of the pulse
pulse_dur = 200e-3
pulse_amp = 0.75e-9
pulse_delay1 = 20e-3
pulse_delay2 = 1e9
compType = 'SymCompartment'

# Store the l23 cell in a python variable
l23_cell = moose.element('/library/L23_NoHotSpot')

# Re-create the python variable pointing to the gran_cell to limit results just to 
# type compartment (excludes the spines and allows for createDataTables function to
# work properly)
l23_cell = moose.wildcardFind(l23_cell.path + '/' + '#[TYPE=' + compType + ']')

# Create the pulse and apply it to the l23 cell's soma
l23_soma_pulse = u.createPulse(l23_cell[0], 'rollingWave', pulse_dur, pulse_amp, 
                           pulse_delay1, pulse_delay2)

# Create a neutral object to store the data in
l23_data = moose.Neutral('/l23_data')
l23_tables = []
for comp in l23_cell:
    l23_tables.append(u.createDataTables(comp,l23_data,l23_soma_pulse))


# Choose the soma and a representative dendrite to view how voltage changes
# for the granule cell
l23_soma_Vm = l23_tables[0][0]
l23_soma_Iex = l23_tables[0][1]
l23_seg2_dend4_121_Vm = l23_tables[100][0]
l23_seg2_dend4_121_Iex = l23_tables[100][1]


# Plot the simulation
simTime = 0.5
simdt = 2.5e-5
plotdt = 0.25e-3
for i in range(10):
    moose.setClock(i, simdt)

moose.setClock(8, plotdt)

# Use the hsolve method for the experiment
u.hsolve(l23_cell[0].parent.path, l23_cell[0].path, simdt)
moose.reinit()
moose.start(simTime)
t = np.linspace(0,simTime,len(l23_soma_Vm.vector))
plt.plot(t,l23_soma_Vm.vector * 1e3, 'r',label = 'l23_soma_Vm (mV)')
plt.plot(t,l23_seg2_dend4_121_Vm.vector * 1e3,'b',label = 'l23_seg2_dend4_121_Vm (mV)')
plt.plot(t, l23_soma_Iex.vector * 1e9, label='l23_soma_Iex (nA)')
plt.plot(t, l23_seg2_dend4_121_Iex.vector * 1e9, label='l23_seg2_dend4_121_Iex (nA)')
plt.xlabel('time (ms)')
plt.legend()
plt.show()
