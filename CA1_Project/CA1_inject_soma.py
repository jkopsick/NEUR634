# Import necessary libraries to perform the experiment and plot results
import moose
import pylab as plt
import numpy as np
import util as u
import plot_channel as pc
import random # to create random synapses to connect randSpikes to
import copy # to create n copies of the synapse dictionary
from chan_proto_part1 import chan_set, rateParams

plt.ion()

# Define the variables needed for the creation of the compartments and pulse
#EREST_ACT = -69e-3 #: Resting membrane potential
#RM_soma = 3.988 # for somatic compartments RM (uniform)
#RM = 3.988*1.54 # for non somatic compartments RM (uniform)
#CM_soma = 1.0153e-6*1e4 # for somatic compartments CM (uniform)
#CM = (1.0153e-6*1e4)/1.54 # for non somatic compartments CM (uniform)
#RA_soma = 2.61 # for non somatic compartments RA (uniform)
#RA = 2.61*1.54 # for non somatic compartments RA (uniform)
#Em = EREST_ACT + 10.613e-3
#initVm = EREST_ACT

EREST_ACT = -69e-3 #: Resting membrane potential
RM_soma = 6.29013 # for somatic compartments RM (uniform)
RM = 3.988*1.54 # for non somatic compartments RM (uniform)
RM_end = 3.1916 # for non somatic compartments RM (uniform)
RM_halfdist= 100.05
RM_slope= 50.48
CM_soma = 1.06e-6*1e4 # for somatic compartments CM (uniform)
CM = (1.06e-6*1e4)*1.54 # for non somatic compartments CM (uniform)
RA_soma = 2.18 # for non somatic compartments RA (uniform)
RA = 2.18 # for non somatic compartments RA (uniform)
Em = EREST_ACT + 10.613e-3
initVm = EREST_ACT

#cond_set = {'Na': 120e-3*1e4, 'K': 36e-3*1e4, 'HCN' : 0e-9*1e12}
#cond_set = {'Na': {(0, 30e-6): 120e-3*1e4, (30e-6, 1) : 0e-3*1e-4}, 
#	    'K': {(0, 30e-6): 0e-3*1e4, (30e-6, 1) : 0e-3*1e-4}, 
#	    'HCN' : {(0, 30e-6): 2e-9*1e12, (30e-6, 1) : 8e-9*1e12}}

# Set of conductances that are being placed non-uniformly but in a discrete fashion for the different
# compartment types in the CA1 morphology. Values have been multipled to reflect conversion from
# physiological to SI units
cond_set = {'Na' : {'soma' : 0e-3*1e4, 'dend' : 0e-3*1e4, 'apical' : 0e-3*1e4}, 
	    'K' : {'soma' : 0e-3*1e4, 'dend' : 0e-3*1e4, 'apical' : 0e-3*1e4}}

# Set the parameters for the pulse applied to the soma (in non-synapse experiments)
pulse_dur = 400e-3
pulse_amp = -30e-12
pulse_delay1 = 0e-3
pulse_delay2 = 1e9

# Define a dictionary for the excitatory synapse channel type (AMPA) based off of information 
# provided in Golding et al. Figure 6
glu = {'name': 'glu', 'Gbar' : 0.1e-9, 'tau1' : 0.5e-3, 'tau2' : 5e-3, 'erev' : 0}

# Create a loop that will create N many pre-synaptic neurons and places them into a dictionary
N = 40
dict = {}
for i in range(1,N+1):
    dict[i] = {'name': 'presyn_' + str(i), 'rate' : 1, 'refractT' : 1e-3, 'delay' : 5e-3}

# Load a multi-compartment model into MOOSE and give it channels
swcfile = 'ri04.CNG.swc'
container = 'CA1_cell'
libraryName = '/library'
compType = 'Compartment'
CA1_cell = u.createMultiCompCell(swcfile, container, libraryName, compType, 					 chan_set, cond_set, rateParams, CaParams = None,
				 CaPoolParams = None, cell_RM = RM, 
				 cell_CM = CM, cell_RA = RA, cell_initVm = initVm, cell_Em = initVm,
				 dist_dep = True)

# Define the variables needed to view the undelying curves for the channel kinetics
plot_powers = True
VMIN = -120e-3
VMAX = 50e-3
CAMIN = 0.01e-3
CAMAX = 40e-3
channelList = ('Na', 'K', 'HCN')

# Graph the curves
for chan in channelList:
        libchan=moose.element('/library/'+chan)
        pc.plot_gate_params(libchan,plot_powers, VMIN, VMAX, CAMIN, CAMAX)

# Re-create the python variable pointing to the gran_cell to limit results just to 
# type compartment (excludes the spines and allows for createDataTables function to
# work properly)
CA1_cell = moose.wildcardFind('/CA1_cell/#[TYPE=Compartment]')

# Acquire the distance and the name of each compartment relative to the origin, and place them into
# a list to be used in the distance-dependent conductance. Origin is close to the soma (soma is 4 um away)
nameList = []
distList = []

for comp in CA1_cell:
     dist, name = u.get_dist_name(comp)
     nameList.append(name)
     distList.append(dist)

distList2 = []
distList2.append(nameList)
distList2.append(distList)

# Work around for now to change the somatic compartments parameters to non spine adjusted values
somaNameListIndex = [i for i, n in enumerate(distList2[0]) if "soma" in n]
somaNameList = [nameList[x] for x in somaNameListIndex]
for soma_comp in somaNameList:
    comp = moose.element('/CA1_cell/' + soma_comp)
    u.setSpecificCompParameters(comp, RM_soma, CM_soma, RA_soma, initVm, initVm)

# Work around for now to change the dendritic compartment parameters to non-uniform values
dendNameListIndex = [i for i, n in enumerate(distList2[0]) if "soma" not in n]
dendNameList = [nameList[x] for x in dendNameListIndex]
for dend_comp in dendNameList:
    comp = moose.element('/CA1_cell/' + dend_comp)
    u.setSpecificCompParametersNonUniform(comp = comp, RM_soma = RM_soma, RM_end = RM_end, 
					  RM_halfdist = RM_halfdist, RM_slope = RM_slope,
					  CM = CM, RA = RA, initVm = initVm, E_leak = initVm)


cond_set_test = {'HCN' : 0e-12*1e12}
minq=2.1582  		# units are pS/um2
maxq=1.5969  		# units are pS/um2
qhalfdis=98.753
qsteep=50.07
# Work around for now to implement non uniform conductance
for comp in nameList:
    comp = moose.element('/CA1_cell/' + comp)
    dist, _ = u.get_dist_name(comp)
    for chan_name, cond in cond_set_test.items():
        SA = np.pi*comp.length*comp.diameter
        proto = moose.element(libraryName + '/' + chan_name)
        chan = moose.copy(proto, comp, chan_name)[0]
	dist = dist*1e6
	hpoint=minq+(maxq-minq)/(1+np.exp(-(dist-qhalfdis)/qsteep))
        chan.Gbar = hpoint*SA*0
        m = moose.connect(chan, 'channel', comp, 'channel')

# Selecting 40 random AMPA synapses to implement that are along the apical dendrite

# First create a list of the indexes that contain apical dendrites
apicalindexList = [i for i, n in enumerate(distList2[0]) if "apical" in n]
# create a new list that uses the indices of the apicalindexList to get a list of their values
apicaldistList =[distList[x] for x in apicalindexList]
# find all apical dendrites that are between 242 and 390 um away from the soma
apicalpotentialsynList = [x for x in apicaldistList if (x >= 242e-6) & (x <=390e-6)]
# get the index for these potential apical dendrites
apicalpotentialSynListIndex = [i for i, x in enumerate(distList2[1]) if x in apicalpotentialsynList]
# get the names for these potential apical dendrites
apicalpotentialnameList = [nameList[x] for x in apicalpotentialSynListIndex]

# pick 40 random values from this list
randsynList = random.sample(apicalpotentialnameList, N)

# Create the pulse and apply it to the granule cell's soma
CA1_soma_pulse = u.createPulse(CA1_cell[0], 'rollingWave', pulse_dur, pulse_amp, 
                           pulse_delay1, pulse_delay2)

# Create a neutral object to store the data in
CA1_data = moose.Neutral('/CA1_data')
CA1_tables = []
for comp in CA1_cell:
    CA1_tables.append(u.createDataTables(comp,CA1_data,CA1_soma_pulse))

# Choose 3 of the compartments along the apical dendrite to monitor their voltage changes
synPlotListName = random.sample(randsynList,3)
indexforSynPlots = [x for x, n in enumerate(distList2[0]) if n in synPlotListName]


# Create a variable to store the voltage table generated for the soma to be used in the plot
CA1_soma_Vm = CA1_tables[0][0]

# Choose voltage tables for apical dendrites you are interested in recording from
synPlotTables = []
for i in indexforSynPlots:
    table = CA1_tables[i][0] 
    synPlotTables.append(table)

# Retrieve the names and put them in a format that allows for them to be labeled neatly when plotting
synPlotTableNames = []
for i in indexforSynPlots:
    path = CA1_tables[i][0].path
    path = path.split('/')[-1]
    path = path.strip(']')
    path = path.replace('[','')
    path = path[:-1] # remove trailing zero from the name of the voltage table
    synPlotTableNames.append(path)

# Plot the simulation
simTimenonSyn = 600e-3 # simulation time for non-synapse simulations
simdt = 25e-6
plotdt = 0.2e-3
for i in range(10):
    moose.setClock(i, simdt)

moose.setClock(8, plotdt)

# Use the hsolve method for the experiment
u.hsolve(CA1_cell[0], simdt)
moose.reinit()
moose.start(simTimenonSyn)
plt.figure()
t = np.linspace(0,simTimenonSyn,len(CA1_soma_Vm.vector))
plt.plot(t,CA1_soma_Vm.vector * 1e3, 'r',label = 'CA1_soma_Vm (mV)')
plt.plot(t,synPlotTables[0].vector * 1e3, 'k',label = synPlotTableNames[0] + ' (mV)')
plt.plot(t,synPlotTables[1].vector * 1e3, 'b',label = synPlotTableNames[1] + ' (mV)')
plt.plot(t,synPlotTables[2].vector * 1e3, 'g',label = synPlotTableNames[2] + ' (mV)')
plt.xlabel('time (ms)')
plt.legend()
