# Import the modules needed to run the simulation in NEURON
from neuron import h, gui
import ballandstick as bs
import numpy as np
from collections import namedtuple
import util as u

# Define the properties to be used by the ball and stick neuron
cell_geometry = u.cell_geometry(soma_d = 12.6157, soma_l = 12.6157, dend_d = 1, dend_l = 100, dend_nseg = 101)

# Define the number of Ball and Stick Neurons for the simulation
N = 3
var = 0.01
variance = np.random.standard_normal(N)*var
bs_properties = []
for i in range(N):
    properties = {}
    properties['cell_geometry'] = cell_geometry
    # The biophysics must be defined in this way to allow for variance in the conductances
    properties['cell_biophysics'] = u.cell_biophysics(Ra = 100, Cm = 1, gna_max = 0.12 + 0.12*variance[i], gk_max = 0.036 + 0.036*variance[i], gl_max = 0.0003 + 0.0003*variance[i], e_leak = -54.3, g_pas = 0.001, e_pas = -65)
    bs_properties.append(properties)

# Creation of a ball and stick neuron using the ball and stick function (credit: NEURON developers)
list_of_bs = [bs.BallAndStick(props) for props in bs_properties]

ENetCon = []
INetCon = []
Esyn = []
Isyn = []
for i in range(N):
    source = list_of_bs[i]
    target = list_of_bs[(i+1)%N]
    esyn = u.createEsyn(target.dend(0.5), 0, 0.1)
    Esyn.append(esyn)
    Enetcon = u.connectNetStim(source.soma(0.5)._ref_v, esyn, 0.02, 5,  source.soma)
    ENetCon.append(Enetcon)
    
    isyn = u.createEsyn(target.dend(0.5), -80, 0.1)
    Isyn.append(isyn)
    Inetcon = u.connectNetStim(source.soma(0.5)._ref_v, isyn, 0.08, 5, source.soma)
    INetCon.append(Inetcon)


# Creation of an alpha synapse so that the small network has an input
asyn = u.createAsyn(list_of_bs[0].dend(0.5), 0, 0.01, 20, 1)





# Set the recording location of the soma to be used in the simulation
soma1_loc = list_of_bs[0].soma(0.5)
soma2_loc = list_of_bs[1].soma(0.5)
soma3_loc = list_of_bs[2].soma(0.5)
