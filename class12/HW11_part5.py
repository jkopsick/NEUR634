# Main script file for part 5 of HW11


# Import the modules needed to run the simulation in NEURON
from neuron import h, gui
import ballandstick as bs
import util as u
import pylab as plt

# Allow for immediate viewing of the figure without declaring a plt.show()
plt.ion() 

# Import parameters needed to run the simulation -- Ball and Stick neurons are created in here as well
import esyn_settings_part5 as syn_set

# Set up recording tables for the experiment for each soma in the ball and stick neuron ring
soma1_v_vec = u.recordCompNeuron(syn_set.soma1_loc)
soma2_v_vec = u.recordCompNeuron(syn_set.soma2_loc)
soma3_v_vec = u.recordCompNeuron(syn_set.soma3_loc)
t_vec = u.recordTimeNeuron()

# Run the simulation
u.simulate(tstop=100)

# Plot results from the simulation
plt.figure()
plt.plot(t_vec, soma1_v_vec, 'r', label = 'soma1 Vm (mv)')
plt.plot(t_vec, soma2_v_vec, 'b', label = 'soma2 Vm (mv)')
plt.plot(t_vec, soma3_v_vec, 'g', label = 'soma3 Vm (mv)')
plt.legend()
