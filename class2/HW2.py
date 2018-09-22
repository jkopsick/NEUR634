# Import necessary libraries to perform the experiment and plot results
import moose
import numpy as np
import pylab

# Create an object to store the simple hierarchy, and the compartment object soma -- must use characters around compartment name
def createCompartment(compname,length,diameter):
	# Create a neutral directory for neuron and create a 		  compartment with the name defined by the user
	neuron = moose.Neutral('/neuron')
	compname = moose.Compartment('/neuron/%s' % compname)
	# Define the length, diameter, total resistivity and 		  compacitance and surface area
	len = length
	dia = diameter
	RM = 20000 * (1/100.)**2 
	CM = 1.0 * (1e-6*100**2)
	#Xarea = np.pi*rad*rad where rad is radius passed in
	SA = np.pi*dia*len
	# Setting the resistance and capacitance for the compartment
	compname.Rm = RM/SA
	compname.Cm = CM*SA
	# Set values for other necessary values for the experiment
	compname.initVm = -65e-3
	compname.Em = -65e-3
	return compname

# Create the pulse that will be connected to the compartment
def createPulse(compname,pulsename,duration,amplitude,delay):
	pulsename = moose.PulseGen('%s' % pulsename)
	# Define pulse duration, amplitude, and delay (if desired)
	pulsename.width[0] = duration
	pulsename.level[0] = amplitude
	pulsename.delay[0] = delay
	# Connect the pulse with the compartment defined by user
	moose.connect(pulsename, 'output', compname, 'injectMsg')
	return pulsename


# Main code
def main(args=None):
	l = createCompartment('swagginWagon',30e-6,20e-6)
	n = createPulse(l, 'rollingWave', 100e-3, 0.1e-9, 50e-3)
	out = moose.Neutral('/out')
	vmtab = moose.Table('/out/somaVm')
	moose.connect(vmtab, 'requestOut', l, 'getVm')
	# Set the current injection to zero to see sole effect of pulse
	l.inject = 0
	# Set the clock so that all components of the simulation run 		  appropriately
	#l.tick = 2
	#moose.setClock(2, 0.25e-3)
	# Run the experiment
	moose.reinit()
	moose.start(0.3)
	# Compute delta VM at the steady state
	delta_Vm = l.Rm*n.level[0]*(1-np.exp(-0.3/(l.Rm*l.Cm)))
	print delta_Vm
	# Compute actual delta Vm at the steady state
	actual_delta_Vm = np.max(vmtab.vector) - np.min(vmtab.vector)
	print actual_delta_Vm
	# Plot the Voltage for the compartment for length of the 		  experiment
	t = pylab.linspace(0, 300e-3, len(vmtab.vector))
	pylab.plot(t, vmtab.vector)
	pylab.show()

if __name__ == '__main__':
	main()
