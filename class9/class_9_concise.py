import moose
import util as u
import pylab as plt
import numpy as np

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

# Create data tables to store the voltage for the soma and each compartment
# making up the dendritic branch
data = moose.Neutral('/data')
soma_pulse = u.createPulse(soma, 'rollingWave', pulse_dur, pulse_amp, pulse_delay, pulse_delay2)
soma_Vm, soma_Iex = u.createDataTables(soma, data, soma_pulse)
dend_Vm, dend_Iex = u.createDataTables(dend, data, soma_pulse)
synchan = moose.SynChan('/neuron/swagginDend/glu')
synchan.Gbar = 1e-8
synchan.tau1 = 2e-3
synchan.tau2 = 2e-3
synchan.Ek = -0.010
msg = moose.connect(dend, 'channel', synchan, 'channel')
sh = moose.SimpleSynHandler('/neuron/swagginDend/glu/synhandler')
moose.connect(sh, 'activationOut', synchan, 'activation')
sh.synapse.num = 1
sh.synapse[0].delay = 1e-3

moose.showmsg(sh)
moose.showmsg(synchan)
moose.showfield(sh.synapse[0])

# the commands above specify and connect the randSpike to a synapse on the 
# dendrite
pre_syn = moose.RandSpike('presyn_input')
pre_syn.rate = 10
pre_syn.refractT = 1e-3
sh.synapse.num = 1
sh.synapse[0].delay = 5e-3
moose.connect(pre_syn, 'spikeOut', sh.synapse[0], 'addSpike')

# commands above connect another spike from a pre-synaptic neuron to the 
# synapse
pre_syn2 = moose.RandSpike('pre_syn2')
pre_syn2.rate = 0
pre_syn2.refractT = 1e-3
sh.synapse.num=2
sh.synapse[1].delay = 5e-3
moose.connect(pre_syn2, 'spikeOut', sh.synapse[1], 'addSpike')

moose.showfield(sh)
sh.synapse.num
spike_table = moose.Table('st')
moose.connect(pre_syn, 'spikeOut', spike_table, 'spike')
spike_table2 = moose.Table('st2')
moose.connect(pre_syn2, 'spikeOut', spike_table2, 'spike')
moose.reinit()
moose.start(1)
moose.showfield('st2')
spike_table.vector
spike_table2.vector

t=np.linspace(0,1, len(dend_Vm.vector))
plt.ion()
for rate in [1,10]:
    pre_syn.rate = rate
    moose.reinit()
    moose.start(1)     
    plt.plot(t,dend_Vm.vector, label=str(rate))

plt.legend()
