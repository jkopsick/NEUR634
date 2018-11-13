# Import necessary libraries to perform the experiment and plot results
import moose
import pylab as plt
import numpy as np
import util as u
import plot_channel as pc
from chan_proto_part1 import chan_set, rateParams, HCNParams
plt.ion()

# Import hyperpolarization activated current from Ascoli, which takes into account the form for the 
# channel that is utilized in the Golding/Spruston/Kath paper (Magee 1998 channel)
nmlchannelfile = 'CA1_cell'
reader = moose.mooseReadNML2(nmlchannelfile)

# Define the variables needed for the creation of the compartments and pulse
EREST_ACT = -70e-3 #: Resting membrane potential
RM = 1/(0.3e-3*1e4)
CM = 1e-6*1e4
RA = 1
Em = EREST_ACT + 10.613e-3
initVm = EREST_ACT
cond_set = {'Na': 0*10000, 'K': 0*2500, 'HCN' : 20000}
pulse_dur = 400e-3
pulse_amp = -50e-12
pulse_delay1 = 20e-3
pulse_delay2 = 1e9

# Load a multi-compartment model into MOOSE and give it channels
swcfile = 'ri04.CNG.swc'
container = 'CA1_cell'
libraryName = '/library'
compType = 'Compartment'
CA1_cell = u.createMultiCompCell(swcfile, container, libraryName, compType, 					 chan_set, cond_set, rateParams, CaParams = None,
				 CaPoolParams = None, HCNParams = HCNParams, cell_RM = RM, 
				 cell_CM = CM, cell_RA = RA, cell_initVm = initVm, cell_Em = Em)

# Define the variables needed to view the undelying curves for the channel kinetics
plot_powers = True
VMIN = -120e-3
VMAX = 50e-3
CAMIN = 0.01e-3
CAMAX = 40e-3
channelList = ('Na','K','HCN')

# Graph the curves
for chan in channelList:
        libchan=moose.element('/library/' + chan)
        pc.plot_gate_params(libchan,plot_powers, VMIN, VMAX, CAMIN, CAMAX)

# Re-create the python variable pointing to the gran_cell to limit results just to 
# type compartment (excludes the spines and allows for createDataTables function to
# work properly)
CA1_cell = moose.wildcardFind('/CA1_cell/#[TYPE=Compartment]')

# Acquire the distance and the name of each compartment relative to the soma, and place them into
# a list to be used in the distance-dependent conductance
nameList = []
distList = []

for comp in CA1_cell:
     dist, name = u.get_dist_name(comp)
     nameList.append(name)
     distList.append(dist)

distList2 = []
distList2.append(nameList)
distList2.append(distList)

minDist = min(distList2[1])
maxDist = max(distList2[1])
indexMinDist = distList2[1].index(minDist)
indexMaxDist = distList2[1].index(maxDist)
nameMinDist = distList2[0][indexMinDist]
nameMaxDist = distList2[0][indexMaxDist]

# Get all the compartments that have a name of dend
[x for x in distList2[0] if "dend" in x]

# Create the pulse and apply it to the granule cell's soma
CA1_soma_pulse = u.createPulse(CA1_cell[0], 'rollingWave', pulse_dur, pulse_amp, 
                           pulse_delay1, pulse_delay2)

# Create a neutral object to store the data in
CA1_data = moose.Neutral('/CA1_data')
CA1_tables = []
for comp in CA1_cell:
    CA1_tables.append(u.createDataTables(comp,CA1_data,CA1_soma_pulse))

# Choose the soma and a representative dendrite to view how voltage changes
# for the granule cell
CA1_soma_Vm = CA1_tables[0][0]
CA1_soma_Iex = CA1_tables[0][1]
CA1_ap_61_20_Vm = CA1_tables[2500][0]
CA1_ap_61_20_Iex = CA1_tables[2500][1]
CA1_ap_121_50_Vm = CA1_tables[5000][0]
CA1_ap_121_50_Iex = CA1_tables[5000][1]
CA1_ap_167_86_Vm = CA1_tables[7500][0]
CA1_ap_167_86_Iex = CA1_tables[7500][1]

# Plot the simulation
simTime = 1000e-3
simdt = 2.5e-5
plotdt = 0.25e-3
for i in range(10):
    moose.setClock(i, simdt)

moose.setClock(8, plotdt)

# Use the hsolve method for the experiment
u.hsolve(CA1_cell[0].parent.path, CA1_cell[0].path, simdt)
moose.reinit()
moose.start(simTime)
t = np.linspace(0,simTime,len(CA1_soma_Vm.vector))
plt.plot(t,CA1_soma_Vm.vector * 1e3, 'r',label = 'CA1_soma_Vm (mV)')
plt.plot(t,CA1_ap_61_20_Vm.vector * 1e3, 'k',label = 'CA1_ap_61_20_Vm (mV)')
plt.plot(t,CA1_ap_121_50_Vm.vector * 1e3, 'b',label = 'CA1_ap_121_50_Vm (mV)')
plt.plot(t,CA1_ap_167_86_Vm.vector * 1e3, 'g',label = 'CA1_ap_167_86_Vm (mV)')
plt.xlabel('time (ms)')
plt.legend()
plt.show()
