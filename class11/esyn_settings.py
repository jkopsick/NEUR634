# This file contains settings needed to create and set Esyn parameters on the proximal, middle,
# and distal dendrite.

# Minimum synaptic weight to produce an AP for distal: 0.020971 (40 ms simulation)
# Minimum synaptic weight to produce an AP for middle: 0.016135 (40 ms simulation)
# Minimum synaptic weight to produce an AP for proximal: 0.012684 (40 ms simulation)

# Do I see five spikes with interval of 10 ms with same minimum weight to ellicit spiking?
# proximal = no, just one; decreasing to interval of 1 ms I still see one spike
# middle = no, just one; decreasing to interval of 1 ms I still see one spike
# distal = no, just one; decreasing to interval of 1 ms I still see one spike

# Import the modules needed to run the simulation in NEURON
from neuron import h, gui
import ballandstick as bs

# Creation of a ball and stick neuron using the ball and stick function (credit: NEURON developers)
cell = bs.BallAndStick()

# Set the recording location of the soma to be used in the simulation
soma_loc = cell.soma(0.5)

# Set the Esyn parameters to be utilized for the simulation
esyn_prox_loc = cell.dend(0.1)
esyn_prox_erev = 0
esyn_prox_tau = 0.1

esyn_mid_loc = cell.dend(0.5)
esyn_mid_erev = 0
esyn_mid_tau = 0.1

esyn_dist_loc = cell.dend(0.9)
esyn_dist_erev = 0
esyn_dist_tau = 0.1

# Set the NetStim parameters to be utilized for the simulation
prox_stim_num = 1
prox_stim_start = 5
prox_stim_length = 10

mid_stim_num = 1
mid_stim_start = 5
mid_stim_length = 10

dist_stim_num = 1
dist_stim_start = 5
dist_stim_length = 10

# Set the NetCon Parameters to be utilized for the simulation
prox_netcon_delay = 5
prox_netcon_weight = 0.012684

mid_netcon_delay = 5
mid_netcon_weight = 0

dist_netcon_delay = 5
dist_netcon_weight = 0
