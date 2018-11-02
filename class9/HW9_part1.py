# Import all libraries necessary to run the file
import numpy as np
import matplotlib.pyplot as plt
import moose
import util as u
from chan_proto_part1 import chan_set, rateParams

neuron = moose.Neutral('/neuron')
spikegen = moose.SpikeGen('/neuron/spikegen')
spikegen.threshold = 0.0
spikegen.refractT = 1e-3
soma_l = 30e-6
soma_rad = 10e-6
dend_l = 100e-6
dend_rad = 1e-6
RM = 2
CM = 0.01
RA = 10
Em = -65e-3
pulse_dur = 100e-3
pulse_amp = 0.1e-9
pulse_delay = 50e-3
pulse_delay2 = 1e9
soma = u.createCompartment(neuron,'swagginSoma',soma_l,soma_rad,RM,CM,RA,Em,Em)
dend = u.createCompartment(neuron,'swagginDend',dend_l,dend_rad,RM,CM,RA,Em,Em)
moose.connect(soma,'axialOut',dend,'handleAxial')

glu = {'name': 'glu', 'Gbar' : 1e-8, 'tau1' : 2e-3, 'tau2' : 2e-3, 'erev' : -0.010}
presyn1 = {'name': 'presyn1', 'rate' : 10, 'refractT' : 1e-3, 'delay' : 5e-3}
presyn2 = {'name': 'presyn2', 'rate' : 0, 'refractT' : 1e-3, 'delay' : 5e-3}
cond_set = {'Na': 120e-3*1e4, 'K': 36e-3*1e4}
libraryName = '/library'
compType = 'Compartment'

# Create the channels necessary for the giant squid axon and place them in a library
u.createChanLib(libraryName, chan_set, rateParams, CaParams = None)

# Create copies of the prototype channels and connect them to the
# compartment, along with their conductances
u.addChannelSet(cond_set, libraryName, neuron, compType)

# Create data tables to store the voltage for the soma and each compartment
# making up the dendritic branch
data = moose.Neutral('/data')
soma_pulse = u.createPulse(soma, 'rollingWave', pulse_dur, pulse_amp, pulse_delay, pulse_delay2)
soma_Vm, soma_Iex = u.createDataTables(soma, data, soma_pulse)
dend_Vm, dend_Iex = u.createDataTables(dend, data, soma_pulse)

glu_handler = u.addSynChan(dend, glu)

# the commands above specify and connect the randSpike to a synapse on the 
# dendrite
#pre_syn = moose.RandSpike('presyn_input')
#pre_syn.rate = 10
#pre_syn.refractT = 1e-3
#glu_handler.synapse.num = 1
#glu_handler.synapse[0].delay = 5e-3
#moose.connect(pre_syn, 'spikeOut', glu_handler.synapse[0], 'addSpike')

presyn1 = u.createRandSpike(presyn1, glu_handler)
# commands above connect another spike from a pre-synaptic neuron to the 
# synapse
#pre_syn2 = moose.RandSpike('pre_syn2')
#pre_syn2.rate = 0
#pre_syn2.refractT = 1e-3
#glu_handler.synapse.num=2
#glu_handler.synapse[1].delay = 5e-3
#moose.connect(pre_syn2, 'spikeOut', glu_handler.synapse[1], 'addSpike')

presyn2 = u.createRandSpike(presyn2, glu_handler)

spike_table = moose.Table('st')
moose.connect(presyn1, 'spikeOut', spike_table, 'spike')
spike_table2 = moose.Table('st2')
moose.connect(presyn2, 'spikeOut', spike_table2, 'spike')
moose.reinit()
moose.start(1)
spike_table.vector
spike_table2.vector

t=np.linspace(0,1, len(dend_Vm.vector))
plt.ion()
for rate in np.arange(10):
    presyn1.rate = rate
    moose.reinit()
    moose.start(1)     
    plt.plot(t,dend_Vm.vector, label=str(rate))

plt.legend()
