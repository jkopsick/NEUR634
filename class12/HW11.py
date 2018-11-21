# NEURON code that uses Esyns on the proximal, middle, and distal dendrites to see how
# location, timing, and conductance of excitatory synapses affects firing rate. 
# Currently using a spherical soma and dendrite with nseg = 101, and these can be changed 
# in ballandstick.py (only values that differ from original BallAndStick class in NEURON is
# for nseg, which was changed from 5 to 101, and for the length of the dendrite, which was changed
# from 200 to 100 to reflect length in previous HWs). Pre-synaptic input comes from NetStim and 
# NetCon objects -- the parameters for these objects, along with those for the Esyns, 
# can be found in the module esyn_settings. Finally, locations for the proximal, middle, and distal
# were as follows: {'proximal' : 0.1, 'middle' : 0.5, 'distal' : 0.9}. You can see these in the 
# esyn_settings file as well.


# Import the modules needed to run the simulation in NEURON
from neuron import h, gui
import ballandstick as bs
import util as u
import pylab as plt

# Allow for immediate viewing of the figure without declaring a plt.show()
plt.ion() 

# Import parameters needed to run the simulation -- Ball and Stick Neuron is created in here as well
import esyn_settings as syn_set

# Set up recording tables for the experiment for the soma and each synapse
soma1_v_vec = u.recordCompNeuron(syn_set.soma1_loc)
soma2_v_vec = u.recordCompNeuron(syn_set.soma2_loc)
soma3_v_vec = u.recordCompNeuron(syn_set.soma3_loc)
t_vec = u.recordTimeNeuron()

# Run the simulation
u.simulate(tstop=40)

# Plot results from the simulation
plt.figure()
plt.plot(t_vec, soma1_v_vec, 'r', label = 'soma1 Vm (mv)')
plt.plot(t_vec, soma2_v_vec, 'b', label = 'soma2 Vm (mv)')
plt.plot(t_vec, soma3_v_vec, 'g', label = 'soma3 Vm (mv)')
plt.legend()
