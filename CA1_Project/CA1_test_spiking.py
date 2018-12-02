# Import necessary libraries to perform the experiment and plot results
import moose
import pylab as plt
import numpy as np
import util as u
import plot_channel as pc
import random # to create random synapses to connect randSpikes to
import copy # to create n copies of the synapse dictionary
from itertools import chain # allow for checking of multiple lists at a time
from chan_proto_CA1 import chan_set, rateParams
from moose.neuroml.NeuroML import NeuroML
import spine_scale_list as ssl

# Allow for each figure created in the simulation to be shown without a plt.show
plt.ion()

# Define some variables that will be utilized in the simulation
EREST_ACT = -69e-3 #: Resting membrane potential
RM_soma = 6.29013*1.25 # for somatic compartments RM (uniform)

RM_end = 3.1916*1.25 # for non somatic compartments RM (uniform)
RM_halfdist= 100.05
RM_slope= 50.48
CM = 1.0595e-6*1e4 # for somatic compartments CM (uniform)
RA = 2.18 # for non somatic compartments RA (uniform)
initVm = EREST_ACT
E_leak = -69e-3

libraryName = '/library'
cell_path = '/library/CA1'

# Read in the CA1 cell morphology
filename = 'CA1_nrn_morph_nobiophys.xml'
neuromlR = NeuroML()
neuromlR.readNeuroMLFromFile(filename)

# Re-create the python variable pointing to the CA1_cell to limit results just to 
# type compartment (excludes the spines and allows for createDataTables function to
# work properly)
CA1_cell = moose.wildcardFind(cell_path + '/' + '#[TYPE=Compartment]')

# Define the origin to be used in the simulation as the center of the soma, as the soma is defined by
# multiple compartments
soma_xloc = CA1_cell[0].x0 + (CA1_cell[0].x - CA1_cell[0].x0)/2 # x for center of compartment
soma_yloc = CA1_cell[0].y0 + (CA1_cell[0].y - CA1_cell[0].y0)/2 # x for center of compartment
soma_zloc = CA1_cell[0].z0 + (CA1_cell[0].z - CA1_cell[0].z0)/2 # x for center of compartment
soma_center = [soma_xloc, soma_yloc, soma_zloc]

# Create library of prototype channels to be copied to compartments now that the cell morphology has been
# loaded in
u.createChanLib(libraryName,chan_set,rateParams,CaParams=None)

# Declare maximal conductances for the classic HH Na and K Channels which will be placed in the soma
cond_set = {'Na' : 140*20, 'K' : 70*20}

cond_set_test = {'HCN' : 0e-12*1e12}
minq=2.1582  		# units are pS/um2
maxq=1.5969  		# units are pS/um2
qhalfdist=98.753
qslope=50.07

# Define a dictionary for the excitatory synapse channel type (AMPA) based off of information 
# provided in Golding et al. Figure 6
AMPA = {'name': 'AMPA', 'Gbar' : 2e-9, 'tau1' : 1.1e-3, 'tau2' : 5.75e-3, 'erev' : 5e-3}

# Create a loop that will create N many pre-synaptic neurons and places them into a dictionary
N = 100
dict = {}
for i in range(1,N+1):
    dict[i] = {'syn_num' : i, 'name': 'presyn_' + str(i), 'rate' : 20, 'refractT' : 1e-3, 'delay' : 20e-3}

# Define the variables needed to view the undelying curves for the channel kinetics
plot_powers = True
VMIN = -120e-3
VMAX = 50e-3
CAMIN = 0.01e-3
CAMAX = 40e-3
channelList = ('Na', 'K', 'HCN')

# Graph the curves
for chan in channelList:
        libchan=moose.element(libraryName + '/' + chan)
        pc.plot_gate_params(libchan,plot_powers, VMIN, VMAX, CAMIN, CAMAX)

# Acquire the distance and the name of each compartment relative to the origin, and place them into
# a list to be used in the distance-dependent conductance. Origin is close to the soma (soma is 4 um away)
nameList = []
distList = []

for comp in CA1_cell:
     dist, name = u.get_dist_name_from_soma(comp,soma_center)
     nameList.append(name)
     distList.append(dist)

distList2 = []
distList2.append(nameList)
distList2.append(distList)

# Implement non-uniform membrane resistances and conductances for each compartment
for comp in distList2[0]:
    comp = moose.element(cell_path + '/' + comp)
    u.setSpecificCompParametersNonUniform(comp = comp, origin = soma_center,
					  RM_soma = RM_soma, RM_end = RM_end, 
					  RM_halfdist = RM_halfdist, RM_slope = RM_slope,
					  CM = CM, RA = RA, initVm = initVm, E_leak = E_leak)

# Implement non-uniform membrane conductances for each compartment
u.setNonUniformConductance(comp_list = distList2[0], cell_path = cell_path, libraryName = libraryName,
			   condSet = cond_set_test, minq = minq, maxq = maxq, qhalfdist = qhalfdist,
			   qslope = qslope, origin = soma_center)

# Create a list of lists of all compartments that are a part of the primary apical dendrite
prim_ap_name_list = []
comp_names = distList2[0]

for i in reversed(ssl.prim_apical):
    name_index = [j for j, n in enumerate(comp_names) if i in n]
    names = [comp_names[x] for x in name_index]
    prim_ap_name_list.append(names)
    comp_names = [x for x in comp_names if x not in names]

# Add HH Na and K channels to somatic compartments
soma_list = prim_ap_name_list[-1]
for comp in soma_list:
    comp = moose.element(cell_path + '/' + comp)
    u.addChanList(libraryName, cond_set, comp)

# Create a list of list of all the names that have no spine scaling, and then scale all of the compartments
# Rm, Cm, and Gbar by the spine scaling factor associated with them
no_spine_name_list = []
comp_names = distList2[0]

for i in reversed(ssl.no_spine_scale):
    name_index = [j for j, n in enumerate(comp_names) if i in n]
    names = [comp_names[x] for x in name_index]
    no_spine_name_list.append(names)
    comp_names = [x for x in comp_names if x not in names]

for comp in list(chain(*no_spine_name_list)):
    comp = moose.element(cell_path + '/' + comp)
    u.scaleCompRmCm(comp = comp, scaling_factor = ssl.scale_factors['no_spine_scale'])
    for chan in cond_set_test:
	chan = moose.element(comp.path + '/' + chan)
        u.scaleChannelCond(chan = chan, scaling_factor = ssl.scale_factors['no_spine_scale'])

# Create a list of list of all the names that have basal spine scaling, and then scale all of the compartments
# Rm, Cm, and Gbar by the spine scaling factor associated with them
basal_name_list = []
comp_names = distList2[0]

for i in reversed(ssl.basal_scale):
    name_index = [j for j, n in enumerate(comp_names) if i in n]
    names = [comp_names[x] for x in name_index]
    basal_name_list.append(names)
    comp_names = [x for x in comp_names if x not in names]

for comp in list(chain(*basal_name_list)):
    comp = moose.element(cell_path + '/' + comp)
    u.scaleCompRmCm(comp = comp, scaling_factor = ssl.scale_factors['basal_scale'])
    for chan in cond_set_test:
	chan = moose.element(comp.path + '/' + chan)
        u.scaleChannelCond(chan = chan, scaling_factor = ssl.scale_factors['basal_scale'])


# Create a list of list of all the names that have med spine rad scaling, and then scale all of the compartments
# Rm, Cm, and Gbar by the spine scaling factor associated with them
med_spine_rad_name_list = []
comp_names = distList2[0]

for i in reversed(ssl.med_spine_rad_scale):
    name_index = [j for j, n in enumerate(comp_names) if i in n]
    names = [comp_names[x] for x in name_index]
    med_spine_rad_name_list.append(names)
    comp_names = [x for x in comp_names if x not in names]

for comp in list(chain(*med_spine_rad_name_list)):
    comp = moose.element(cell_path + '/' + comp)
    u.scaleCompRmCm(comp = comp, scaling_factor = ssl.scale_factors['med_spine_rad_scale'])
    for chan in cond_set_test:
	chan = moose.element(comp.path + '/' + chan)
        u.scaleChannelCond(chan = chan, scaling_factor = ssl.scale_factors['med_spine_rad_scale'])

# Create a list of list of all the names that have med spine LM scaling, and then scale all of the compartments
# Rm, Cm, and Gbar by the spine scaling factor associated with them
med_spine_LM_name_list = []
comp_names = distList2[0]

for i in reversed(ssl.med_spine_LM_scale):
    name_index = [j for j, n in enumerate(comp_names) if i in n]
    names = [comp_names[x] for x in name_index]
    med_spine_LM_name_list.append(names)
    comp_names = [x for x in comp_names if x not in names]

for comp in list(chain(*med_spine_LM_name_list)):
    comp = moose.element(cell_path + '/' + comp)
    u.scaleCompRmCm(comp = comp, scaling_factor = ssl.scale_factors['med_spine_LM_scale'])
    for chan in cond_set_test:
	chan = moose.element(comp.path + '/' + chan)
        u.scaleChannelCond(chan = chan, scaling_factor = ssl.scale_factors['med_spine_LM_scale'])

# Create a list of list of all the names that have max spine rad scaling, and then scale all of the compartments
# Rm, Cm, and Gbar by the spine scaling factor associated with them
max_spine_rad_name_list = []
comp_names = distList2[0]

for i in reversed(ssl.max_spine_rad_scale):
    name_index = [j for j, n in enumerate(comp_names) if i in n]
    names = [comp_names[x] for x in name_index]
    max_spine_rad_name_list.append(names)
    comp_names = [x for x in comp_names if x not in names]

for comp in list(chain(*max_spine_rad_name_list)):
    comp = moose.element(cell_path + '/' + comp)
    u.scaleCompRmCm(comp = comp, scaling_factor = ssl.scale_factors['max_spine_rad_scale'])
    for chan in cond_set_test:
	chan = moose.element(comp.path + '/' + chan)
        u.scaleChannelCond(chan = chan, scaling_factor = ssl.scale_factors['max_spine_rad_scale'])

# Create a list of list of all the names that have thin rad spine scaling, and then scale all of the compartments
# Rm, Cm, and Gbar by the spine scaling factor associated with them
thin_rad_name_list = []
comp_names = distList2[0]

for i in reversed(ssl.thin_rad_spine_scale):
    name_index = [j for j, n in enumerate(comp_names) if i in n]
    names = [comp_names[x] for x in name_index]
    thin_rad_name_list.append(names)
    comp_names = [x for x in comp_names if x not in names]

for comp in list(chain(*thin_rad_name_list)):
    comp = moose.element(cell_path + '/' + comp)
    u.scaleCompRmCm(comp = comp, scaling_factor = ssl.scale_factors['thin_rad_spine_scale'])
    for chan in cond_set_test:
	chan = moose.element(comp.path + '/' + chan)
        u.scaleChannelCond(chan = chan, scaling_factor = ssl.scale_factors['thin_rad_spine_scale'])

# Create a list of list of all the names that have thin LM spine scaling, and then scale all of the compartments
# Rm, Cm, and Gbar by the spine scaling factor associated with them
thin_LM_name_list = []
comp_names = distList2[0]

for i in reversed(ssl.thin_LM_spine_scale):
    name_index = [j for j, n in enumerate(comp_names) if i in n]
    names = [comp_names[x] for x in name_index]
    thin_LM_name_list.append(names)
    comp_names = [x for x in comp_names if x not in names]

for comp in list(chain(*thin_LM_name_list)):
    comp = moose.element(cell_path + '/' + comp)
    u.scaleCompRmCm(comp = comp, scaling_factor = ssl.scale_factors['thin_LM_spine_scale'])
    for chan in cond_set_test:
	chan = moose.element(comp.path + '/' + chan)
        u.scaleChannelCond(chan = chan, scaling_factor = ssl.scale_factors['thin_LM_spine_scale'])

# Create a list of distances all the SLM compartments to be used in selection of the synapse sites
SLM_names = list(chain(list(chain(*med_spine_LM_name_list)), list(chain(*thin_LM_name_list))))
SLM_names_index = [j for j, n in enumerate(distList2[0]) if n in SLM_names]
SLM_dist_list = [distList[x] for x in SLM_names_index]

# find all SLM dendritic compartments that are between 242 and 390 um away from the soma
SLM_potential_syn_list = [x for x in SLM_dist_list if (x >= 242e-6) & (x <=390e-6)]
# get the index for these potential SLM dendritic compartments
SLM_potential_syn_list_index = [i for i, x in enumerate(distList2[1]) if x in SLM_potential_syn_list]
# get the names for these potential SLM dendritic compartments
SLM_potential_syn_name_list = [nameList[x] for x in SLM_potential_syn_list_index]

# pick 40 random values from this list
randsynList = random.sample(SLM_potential_syn_name_list, N)

# Add AMPA synapses to each of these random SLM compartments
synapseHandler = []
for comp in randsynList:
    slm = moose.element(cell_path + '/' + comp)
    sh = u.addSynChan(slm, AMPA)
    synapseHandler.append(sh)

# Connect the presynaptic neuron to each synapse
for i in dict.values():
	shnum = i['syn_num'] - 1 # get the correct index to be used for the synapseHandler
        presyncell = u.createRandSpike(i, synapseHandler[shnum])

# Create a neutral object to store the data in
CA1_data = moose.Neutral('/CA1_data')
CA1_tables = []
for comp in CA1_cell:
    CA1_tables.append(u.createDataTables(compname = comp, data_hierarchy = CA1_data))

# Create a variable to store the voltage table generated for the soma to be used in the plot
CA1_soma_Vm = CA1_tables[0][0]

# Get all the distances associated with the primary apical dendritic tree to be used later for setting up
# synapses and measuring voltage attenuation
prim_ap_names = list(chain(*prim_ap_name_list)) # unordered list of prim apical dendrite compartments
prim_ap_names_index = [j for j, n in enumerate(distList2[0]) if n in prim_ap_names]
prim_ap_dist_list = [distList[x] for x in prim_ap_names_index]
prim_ap_names = [nameList[x] for x in prim_ap_names_index] # ordered list of prim apical dend compartments

# Choose 3 of the compartments along the primary apical dendrite to monitor voltage attenuation. Selective
# dendrites were chosen at roughly 100 um, 225 um, and 350 um (furthest compartment on primary apical dend)
primApicalPlotListName = [prim_ap_names[67], prim_ap_names[133], prim_ap_names[215]]
indexforprimApicalPlots = [x for x, n in enumerate(distList2[0]) if n in primApicalPlotListName]

# Choose voltage tables for apical dendrites you are interested in recording from
primApicalPlotTables = []
for i in indexforprimApicalPlots:
    table = CA1_tables[i][0] 
    primApicalPlotTables.append(table)

# Retrieve the names and put them in a format that allows for them to be labeled neatly when plotting
primApicalPlotTableNames = []
for i in indexforprimApicalPlots:
    path = CA1_tables[i][0].path
    path = path.split('/')[-1]
    path = path.strip(']')
    path = path.replace('[','')
    path = path[:-1] # remove trailing zero from the name of the voltage table
    primApicalPlotTableNames.append(path)

# Plot the simulation
simTime = 1000e-3 # simulation time for non-synapse simulations
simdt = 10e-6
plotdt = 0.1e-3
for i in range(10):
    moose.setClock(i, simdt)

moose.setClock(8, plotdt)

# Use the hsolve method for the experiment
u.hsolve(CA1_cell[0], simdt)
moose.reinit()
moose.start(simTime)
plt.figure()
t = np.linspace(0,simTime,len(CA1_soma_Vm.vector))
plt.plot(t,CA1_soma_Vm.vector * 1e3, 'r',label = 'CA1_soma_Vm (mV)')
plt.plot(t,primApicalPlotTables[0].vector * 1e3, 'k',label = primApicalPlotTableNames[0] + ' (mV)')
plt.plot(t,primApicalPlotTables[1].vector * 1e3, 'b',label = primApicalPlotTableNames[1] + ' (mV)')
plt.plot(t,primApicalPlotTables[2].vector * 1e3, 'g',label = primApicalPlotTableNames[2] + ' (mV)')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.legend()
