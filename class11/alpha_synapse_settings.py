# This file contains settings needed to create and set alpha synapse parameters on the proximal, middle,
# and distal dendrite. It is no longer needed for HW 11 after the assignment was updated, but this code
# may be useful in future HW assignments involving NEURON

# Import the modules needed to run the simulation in NEURON
from neuron import h, gui
import ballandstick as bs

# Creation of a ball and stick neuron using the ball and stick function (credit: NEURON developers)
cell = bs.BallAndStick()

# Set the recording location of the soma to be used in the simulation
soma_loc = cell.soma(0.5)

# Set the recording location and parameters for the proximal, middle, and distal parts of the dendrite
asyn_prox_loc = cell.dend(0.1)
asyn_prox_erev = 0
ayn_prox_gmax = 7.89075e-4
asyn_prox_onset = 20
asyn_prox_tau = 1

asyn_mid_loc = cell.dend(0.5)
asyn_mid_erev = 0
ayn_mid_gmax = 0
asyn_mid_onset = 20
asyn_mid_tau = 1

asyn_dist_loc = cell.dend(0.9)
asyn_dist_erev = 0
ayn_dist_gmax = 0
asyn_dist_onset = 20
asyn_dist_tau = 1
