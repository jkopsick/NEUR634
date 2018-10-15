import moose
import numpy as np
from collections import namedtuple
#from neuron import h

# Create an object to store the simple hierarchy, and the compartment object soma -- must use characters around compartment name
def createCompartment(directory,compname,length,radius,RM,CM,RA,initVm,Em):
    # Create a neutral directory for neuron and create a           compartment with the name defined by the user
    compname = moose.Compartment(directory.path + '/' + compname)
    # Define the x-sectional area, diameter, and surface area
    Xarea = np.pi*radius*radius
    dia = 2*radius
    SA = np.pi*dia*length
    # Set the length and diameter of the compartment
    compname.length = length
    compname.diameter = dia
    # Set the resistance, capacitance, and axial resistance of compartment
    # using the total resistivity and capacticance, and x-sectional area
    compname.Rm = RM/SA
    compname.Cm = CM*SA
    compname.Ra = RA*length/Xarea
    #compname.Ra = RA*compname.length/Xarea
    #RM = 20000 * (1/100.)**2 this would be just 2 if in SI already
    #CM = 1.0 * (1e-6*100**2)
    # Setting the resistance and capacitance for the compartment
    # Set values for other necessary values for the experiment
    compname.initVm = initVm
    compname.Em = Em
    return compname

# Create the pulse that will be connected to the compartment
def createPulse(compname,pulsename,duration,amplitude,delay1,delay2):
    pulsename = moose.PulseGen('%s' % pulsename)
    # Define pulse duration, amplitude, and delay (if desired)
    pulsename.width[0] = duration
    pulsename.level[0] = amplitude
    pulsename.delay[0] = delay1
    pulsename.delay[1] = delay2
    # Connect the pulse with the compartment defined by user
    moose.connect(pulsename, 'output', compname, 'injectMsg')
    return pulsename

# Create data tables that will store the current and voltage associated with a compartment
def createDataTables(compname,data_hierarchy,pulsename):
    # Create a new path using string manipulation to be used in the creation
    # of a unique data table for each compartment
    comp_path = compname.path.split('/')[-1]
    comp_path = comp_path.strip(']')
    comp_path = comp_path.replace('[','')
    # Create the unique membrane potential table and connect this to the compartment
    Vmtab = moose.Table(data_hierarchy.path + '/' + comp_path + '_Vm')
    moose.connect(Vmtab, 'requestOut', compname, 'getVm')
    # Create the unique external current table and connect this to the compartment
    current_tab = moose.Table(data_hierarchy.path + '/' + comp_path + '_Iex')
    moose.connect(current_tab, 'requestOut', pulsename, 'getOutputValue')
    return Vmtab, current_tab

# Function that will discretize a compartment into N many segments
def discretize(modelLoc,numComps,length,radius,RM,CM,RA,Em):
    # Create an array of n compartments as designated by user
    compArray = moose.vec('%s/comp' % (modelLoc.path), n=numComps,
                          dtype='Compartment')
    # Define the x-sectional area, diameter, and surface area
    Xarea = np.pi*radius*radius
    dia = 2*radius
    SA = np.pi*dia*length
    # Loop through created compartments and assign resistance, capacitance
    # and axial resistance
    for comp in compArray:
        comp.Rm = RM/SA
        comp.Cm = CM*SA
        comp.Ra = RA*length/Xarea
        comp.initVm = Em
        comp.Em = Em
    # Connect the compartments
    for i in range(len(compArray[0:-2])):
        moose.connect(compArray[i],'axialOut',compArray[i+1],'handleAxial')
    return compArray

# Function that will set the parameters for all compartments in a model neuron
def setCompParameters(compvector,comptype,RM,CM,RA,initVm,E_leak):
    # Loop through all the compartments in the vector
    for comp in moose.wildcardFind(compvector.path + '/' + '#[TYPE=' + comptype + ']'):
    	SA = np.pi*comp.length*comp.diameter
	X_area = np.pi*comp.diameter*comp.diameter/4.0
	comp.Rm = RM/SA
	comp.Cm = CM*SA
	comp.Ra = RA*comp.length/X_area
	comp.initVm = initVm
	comp.Em = E_leak

# Function that will create a pulse in NEURON
def createNeuronPulse(compname,pulsename,duration,amplitude,delay):
    pulsename = h.IClamp(compname)
    pulsename.dur = duration
    pulsename.amp = amplitude
    pulsename.delay = delay
    return pulsename

# Function that will record the voltage for a compartment in NEURON
def recordCompNeuron(location):
    # Store the time and voltage in a vector
    v_vec = h.Vector()
    v_vec.record(location._ref_v)
    return v_vec

# Function that will record the time in a NEURON simulation
def recordTimeNeuron():
    t_vec = h.Vector()
    t_vec.record(h._ref_t)
    return t_vec

# These named tuples provide a template for creating the biophysics associated with a given
# channel type
AlphaBetaChanParams = namedtuple('AlphaBetaChannelParams',
                         'A_rate A_B A_C A_vhalf A_vslope B_rate B_B B_C B_vhalf B_vslope')
ChannelSettings = namedtuple('ChannelSettings', 'Xpow Ypow Erev name Xparam Yparam')

# Function that will create from the biophysics defined for a channel a MOOSE channel
def createChanProto(libraryName, channelParams, rateParams):
    # Create a library to store the channel prototypes
    if not moose.exists(libraryName):
        lib = moose.Neutral(libraryName)
    else:
        lib = moose.element(libraryName)
    
    # Create the channel and set the powers and reversal potential
    channel = moose.HHChannel(lib.path + '/' + channelParams.name)
    channel.Ek = channelParams.Erev
    channel.Xpower = channelParams.Xpow
    channel.Ypower = channelParams.Ypow
    
    # Define the activation gating kinetics if they exist
    if channel.Xpower > 0:
        xGate = moose.HHGate(channel.path + '/' + 'gateX')
        xGate.setupAlpha(channelParams.Xparam + rateParams)
        
    # Define the inactivation gating kinetics if they exist    
    if channel.Ypower > 0:
        yGate = moose.HHGate(channel.path + '/' + 'gateY')
        yGate.setupAlpha(channelParams.Yparam + rateParams)
    
    # Set the tick for the channel so that it is set appropriately when it
    # is copied to a channel from the library    
    channel.tick = -1
    
    return channel

# Function that will create channels using the createChanProto function and place them
# in a MOOSE library, which will be utilized in adding channels to a particular compartment
def createChanLib(libraryName, channelSet, rateParams):
    # Create a library to store the channel prototypes
    if not moose.exists(libraryName):
        lib = moose.Neutral(libraryName)
    else:
        lib = moose.element(libraryName)

    # Add all the channels to the MOOSE library
    for params in channelSet.values():
        chan = createChanProto(libraryName, params, rateParams)

# Function that will create a multi-compartment model in MOOSE from a .p or .swc file
def createMultiCompCell(file_name, container_name, library_name, comp_type, channelSet, condSet,
                        rateParams, cell_RM = None, cell_CM = None, cell_RA = None, 
                        cell_initVm = None, cell_Em = None):
    # Create the channel types and store them in a library to be used by each compartment
    # in the model
    createChanLib(library_name, channelSet, rateParams)
    # Load in the model in question
    if file_name.endswith('.p'):
        cell = moose.loadModel(file_name, container_name)
        for comp in moose.wildcardFind(cell.path + '/' + '#[TYPE=' + comp_type + ']'):
            for chan_name, cond in condSet.items():
                SA = np.pi*comp.length*comp.diameter
                proto = moose.element(library_name + '/' + chan_name)
                chan = moose.copy(proto, comp, chan_name)[0]
                chan.Gbar = gbar*SA
                m = moose.connect(chan, 'channel', comp, 'channel')
    else:
        cell = moose.loadModel(file_name, container_name)
        setCompParameters(cell, comp_type, cell_RM, cell_CM, cell_RA, cell_initVm, cell_Em)
        for comp in moose.wildcardFind(cell.path + '/' + '#[TYPE=' + comp_type + ']'):
            for chan_name, cond in condSet.items():
                SA = np.pi*comp.length*comp.diameter
                proto = moose.element(library_name + '/' + chan_name)
                chan = moose.copy(proto, comp, chan_name)[0]
                chan.Gbar = cond*SA
                m = moose.connect(chan, 'channel', comp, 'channel')
    return cell
