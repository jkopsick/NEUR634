# This file contains code that will create ball and stick neurons in a ring that have
# excitatory syns and NetCons, along with an AlphaSynapse to allow for input current
# for the simulation

# Import the modules needed to run the simulation in NEURON
from neuron import h, gui
import ballandstick as bs
import numpy as np
import util as u

# Define the properties for the structure of the ball and stick neurons
soma_d = 12.6157
soma_l = 12.6157
dend_d = 1
dend_l = 100
dend_nseg = 101

soma_loc = 0.5 # recording location for the soma

# Define the biophysics for the ball and stick neurons
Ra = 100
Cm = 1
gna_max = 0.12
gk_max = 0.036
gl_max = 0.0003
e_leak = -54.3
g_pas = 0.001 # passive conductance for the dendrite
e_pas = -65 # passive e_rev for the dendrite

# Define the properties to be used by the ball and stick neuron
cell_geometry = u.cell_geometry(soma_d = soma_d, soma_l = soma_l, dend_d = dend_d, 
				dend_l = dend_l, dend_nseg = dend_nseg)

# Define the number of Ball and Stick Neurons for the simulation
N = 3

# Define and set the variance to be used in the simulation
var = 0.2
variance = np.random.standard_normal(N)*var

# Defining the structure and biophysics with varied conductances utilizing a list of dictionaries. This
# list will then be used to create the ring of Ball and Stick neurons 
bs_properties = []
for i in range(N):
    properties = {}
    properties['cell_geometry'] = cell_geometry
    # The biophysics must be defined in this way to allow for variance in the conductances
    properties['cell_biophysics'] = u.cell_biophysics(Ra = Ra, Cm = Cm, 
						      gna_max = gna_max + gna_max*variance[i], 
						      gk_max = gk_max + gk_max*variance[i], 
						      gl_max = gl_max + gl_max*variance[i], 
						      e_leak = e_leak, g_pas = g_pas, e_pas = e_pas)
    bs_properties.append(properties)

# Creation of ball and stick neurons using the ball and stick function using list comprehension
# (credit: NEURON developers)
list_of_bs = [bs.BallAndStick(properties) for properties in bs_properties]

# Creation of the excitatory Esyns and NetCons to be utilized in the simulation
Esyn = []
ENetCon = []
for i in range(N):
    source = list_of_bs[i]
    target = list_of_bs[(i+1)%N]
    esyn = u.createEsyn(synapse_location = target.dend(0.5), erev = 0, tau = 0.1)
    Enetcon = u.connectNetStim(stim = source.soma(0.5)._ref_v, synapse = esyn, delay = 5,
			       weight = 0.02, section = source.soma)
    Esyn.append(esyn)
    ENetCon.append(Enetcon)


# Creation of an alpha synapse so that the small network has an input
asyn = u.createAsyn(synapse_location = list_of_bs[0].dend(0.5), erev = 0, gmax = 0.01, 
		    onset = 20, tau = 1)

# Set the recording location of the soma to be used in the simulation
soma1_loc = list_of_bs[0].soma(soma_loc)
soma2_loc = list_of_bs[1].soma(soma_loc)
soma3_loc = list_of_bs[2].soma(soma_loc)
