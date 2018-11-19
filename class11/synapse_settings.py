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

# NetStim parameters set to ones that were utilized in class and stated in HW
prox_stim_num = 1
prox_stim_start = 5
prox_stim_length = 1

mid_stim_num = 1
mid_stim_start = 5
mid_stim_length = 1

dist_stim_num = 1
dist_stim_start = 5
dist_stim_length = 1

# NetCon parameters set to ones that were utilized in class
prox_netcon_delay = 5
prox_netcon_weight = 0.04

mid_netcon_delay = 5
mid_netcon_weight = 0.04

dist_netcon_delay = 5
dist_netcon_weight = 0.04
