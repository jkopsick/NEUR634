# Import all libraries necessary to run the file
import numpy as np
import matplotlib.pyplot as plt
import moose
import util as u
from chan_proto_part1 import chan_set, rateParams

# Define the parameters used to create the soma and dendrite compartments
EREST_ACT = -70e-3 #: Resting membrane potential  
Em = EREST_ACT + 10.613e-3
initVm = EREST_ACT

soma_l = 50e-6
soma_rad = 15e-6
dend_l = 100e-6
dend_rad = 1e-6

RM = 1/(0.3e-3*1e4)
CM = 1e-6*1e4
RA = 2

# Create the soma and dendrite under the neutral neuron directory, and connect them to each other
neuron = moose.Neutral('/neuron')
soma = u.createCompartment(neuron,'swagginSoma',soma_l,soma_rad,RM,CM,RA,Em,Em)
dend = u.createCompartment(neuron,'swagginDend',dend_l,dend_rad,RM,CM,RA,Em,Em)
moose.connect(soma,'axial',dend,'raxial')

# Define dictionaries for the excitatory and inhibitory channels, and the pre-synaptic neurons
glu = {'name': 'glu', 'Gbar' : 5e-9, 'tau1' : 1e-3, 'tau2' : 5e-3, 'erev' : 0}
GABA = {'name': 'GABA', 'Gbar' : 1e-6, 'tau1' : 3e-3, 'tau2' : 5e-3, 'erev' : -80e-3}
presyn1 = {'name': 'presyn1', 'rate' : 10, 'refractT' : 1e-3, 'delay' : 5e-3}
presyn2 = {'name': 'presyn2', 'rate' : 70, 'refractT' : 1e-3, 'delay' : 5e-3}

# Create HH Na and K channels with their corresponding conductances, and place them in a library
# so that they can be implemented into the compartments (if so desired)
cond_set = {'Na': 120e-3*1e4, 'K': 36e-3*1e4}
libraryName = '/library'
compType = 'Compartment'
u.createChanLib(libraryName, chan_set, rateParams, CaParams = None)

# Create copies of the prototype channels and connect them to each compartment under the neuron hierarchy
u.addChannelSet(cond_set, libraryName, neuron, compType)

# Create data tables to store the voltage for the soma and dendrite, and ignore current output
# since there is no pulse being provided to the soma
data = moose.Neutral('/data')
soma_Vm, _ = u.createDataTables(compname = soma, data_hierarchy = data)
dend_Vm, _ = u.createDataTables(compname = dend, data_hierarchy = data)

# Create the excitatory synchan and return the handler to be used in the creation of the RandSpike
glu_handler = u.addSynChan(dend, glu)
presyn1 = u.createRandSpike(presyn1, glu_handler)

# Create the inhibitory synchan and return the handler to be used in the creation of the RandSpike
GABA_handler = u.addSynChan(dend, GABA)
presyn2 = u.createRandSpike(presyn2, GABA_handler)

# Create tables to store the output time of the RandSpikes
spike_table = moose.Table('st')
moose.connect(presyn1, 'spikeOut', spike_table, 'spike')
spike_table2 = moose.Table('st2')
moose.connect(presyn2, 'spikeOut', spike_table2, 'spike')

# Run the experiment and plot the results
simtime = 1
simdt = 0.25e-5
plotdt = 0.25e-3
for i in range(10):
    moose.setClock(i, simdt)
moose.setClock(8, plotdt)
moose.reinit()
moose.start(simtime)

t=np.linspace(0,simtime, len(dend_Vm.vector))
plt.ion()
plt.plot(t,soma_Vm.vector*1e3, 'r', label = 'soma Vm (mV)')
plt.plot(t,dend_Vm.vector*1e3, 'b', label = 'dend Vm (mV)')
plt.legend()
