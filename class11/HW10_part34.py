# Code that will be utilized in parts 3 and 4 of HW, where we will be using NetStim and NetCon objects
# which I believe will involve the use of Esyns and not Asyns (as it appears Asyns are not amenable
# to NetStim and NetCon objects, as they are activated manually and not by a previously declared
# pre-synaptic neuron)

# Will most likely need to create additional settings for the Esyn in synapse_settings


# Import the modules needed to run the simulation in NEURON
from neuron import h, gui
import ballandstick as bs
import util as u
import pylab as plt

plt.ion() # allows for immediate viewing of the figure without declaring a plt.show()

# Import parameters needed to run the simulation now that the ball and stick neuron has been created
import synapse_settings as syn_set 

# Adding ESyns to the proximal, middle, and distal dendrites
esyn_prox = u.createEsyn(syn_set.esyn_prox_loc, syn_set.esyn_prox_erev, syn_set.esyn_prox_tau)
esyn_mid = u.createEsyn(syn_set.esyn_mid_loc, syn_set.esyn_mid_erev, syn_set.esyn_mid_tau)
esyn_dist = u.createEsyn(syn_set.esyn_dist_loc, syn_set.esyn_dist_erev, syn_set.esyn_dist_tau)

# Create NetStim objects so that ESyns can be connected to them and NetCons
nc_stim_prox = u.createNetStim(syn_set.prox_stim_num, syn_set.prox_stim_start, syn_set.prox_stim_length)
nc_stim_mid = u.createNetStim(syn_set.mid_stim_num, syn_set.mid_stim_start, syn_set.mid_stim_length)
nc_stim_dist = u.createNetStim(syn_set.dist_stim_num, syn_set.dist_stim_start, syn_set.dist_stim_length)

# Create NetCon objects for each Esyn so that they can be connected to the "pre-synaptic" neuron
netcon_prox = u.connectNetStim(nc_stim_prox, esyn_prox, syn_set.prox_netcon_delay,
			       syn_set.prox_netcon_weight)

netcon_mid = u.connectNetStim(nc_stim_mid, esyn_mid, syn_set.mid_netcon_delay,
			       syn_set.mid_netcon_weight)

netcon_dist = u.connectNetStim(nc_stim_dist, esyn_dist, syn_set.dist_netcon_delay,
			       syn_set.dist_netcon_weight)

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
