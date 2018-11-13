import moose
import numpy as np
from matplotlib import pyplot as plt
from pprint import pprint # allows for one to see nested dictionaries clearer
plt.ion()
nmdachan = moose.NMDAChan('/nmda')
moose.showfield(nmdachan)

nmdachan.KMg_A = 0.17
nmdachan.KMg_B = 0.012
nmdachan.CMg = 1.4
 
# Function to plot and see magnesium effect on NMDA channel
def magnesium_term(A,B,C):
     eta = 1.0/A
     gamma = 1.0/B
     v = np.linspace(-100e-3,60e-3,1000) # voltage vector
     mg_term = 1/(1+eta*C*np.exp(-gamma*v))
     plt.plot(v,mg_term)
 

synparams = {}
mgparams = {'A' : (1/6.0), 'B' : (1/80.0), 'conc' : 1.4}
synparams['ampa'] = {'Erev' : 5e-3, 'tau1' : 1.0e-3, 'tau2' : 5e-3, 
                     'Gbar' : 1e-9}
synparams['nmda'] = {'Erev' : 5e-3, 'tau1' : 1.1e-3, 'tau2' : 37.5e-3, 
                     'Gbar' : 2e-9, 'mgparams' : mgparams,
                     'temperature' : 303, 'extCa' : 2,
                     'condFraction' : 0.1}

pprint(synparams, width = 40)

moose.Neutral('cell2')
moose.Compartment('/cell2/dend')

# Looping through synparams dictionary and creates synaptic channels
# Moose uses GHK current equation to solve how calcium affects the NMDA
# channel -- temp, extCa, and condFraction required for solving this
for key, params in synparams.items():
    if key == 'nmda':
        chan = moose.NMDAChan('/cell2/dend/' + key)
        chan.KMg_A = params['mgparams']['A']
        chan.KMg_B = params['mgparams']['B']
        chan.CMg = params['mgparams']['conc']
        chan.temperature = params['temperature']
        chan.extCa = params['extCa']
        chan.condFraction = params['condFraction']
    else:
        chan = moose.SynChan('/cell2/dend/' + key)
    chan.Gbar = params['Gbar']
    chan.tau1 = params['tau1']
    chan.tau2 = params['tau2']
    chan.Ek = params['Erev']




# Below code is just an example -- in HW we will modify our existing CaPool
# code

capool = moose.CaConc('/cell2/dend/capool')
moose.showfield(capool)

nmdachan = moose.element('/cell2/dend/nmda')
moose.connect(nmdachan, 'ICaOut', capool, 'current')
moose.connect(capool, 'concOut', nmdachan, 'setIntCa')

moose.showmsg(nmdachan)


