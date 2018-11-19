# NEURON code that will use Alpha Synapses on the proximal, middle, and distal dendrites to see how
# location, timing, and conductance of excitatory synapses affects firing rate. Currently using a spherical
# soma and dendrite with nseg = 5, and these can be changed in ballandstick.py

# Import the modules needed to run the simulation in NEURON
from neuron import h, gui
import ballandstick as bs
import util as u
import pylab as plt

plt.ion() # allows for immediate viewing of the figure without declaring a plt.show()

# Import parameters needed to run the simulation now that the ball and stick neuron has been created
import synapse_settings as syn_set 


# Adding AlphaSynapses to the proximal, middle, and distal dendrites
asyn_prox = u.createAsyn(syn_set.asyn_prox_loc, syn_set.asyn_prox_erev, syn_set.ayn_prox_gmax,
			 syn_set.asyn_prox_onset, syn_set.asyn_prox_tau) # synapse at proximal dendrite

asyn_mid =  u.createAsyn(syn_set.asyn_mid_loc, syn_set.asyn_mid_erev, syn_set.ayn_mid_gmax,
			 syn_set.asyn_mid_onset, syn_set.asyn_mid_tau) # synapse at middle dendrite

asyn_dist = u.createAsyn(syn_set.asyn_dist_loc, syn_set.asyn_dist_erev, syn_set.ayn_dist_gmax,
			 syn_set.asyn_dist_onset, syn_set.asyn_dist_tau) # synapse at distal dendrite

# Set up recording tables for the experiment for the soma and each synapse
soma_v_vec = u.recordCompNeuron(syn_set.soma_loc)
prox_dend_v_vec = u.recordCompNeuron(syn_set.asyn_prox_loc)
mid_dend_v_vec = u.recordCompNeuron(syn_set.asyn_mid_loc)
dist_dend_v_vec = u.recordCompNeuron(syn_set.asyn_dist_loc)
t_vec = u.recordTimeNeuron()

# Run the simulation
u.simulate(tstop=40)

# Plot results from the simulation
plt.figure()
plt.plot(t_vec, soma_v_vec, 'r', label = 'soma Vm (mv)')
plt.plot(t_vec, prox_dend_v_vec, 'b', label = 'prox dend Vm (mv)')
plt.plot(t_vec, mid_dend_v_vec, 'g', label = 'mid dend Vm (mv)')
plt.plot(t_vec, dist_dend_v_vec, 'k', label = 'dist dend Vm (mv)')
plt.legend()




