import moose
import numpy as np
import pylab as plt
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
def createDataTables(compname,data_hierarchy,pulsename=None):
    # Create a new path using string manipulation to be used in the creation
    # of a unique data table for each compartment
    comp_path = compname.path.split('/')[-1]
    comp_path = comp_path.strip(']')
    comp_path = comp_path.replace('[','')
    comp_path = comp_path[:-1]
    # Create the unique membrane potential table and connect this to the compartment
    Vmtab = moose.Table(data_hierarchy.path + '/' + comp_path + '_Vm')
    moose.connect(Vmtab, 'requestOut', compname, 'getVm')
    # Create the unique external current table and connect this to the compartment
    if pulsename is not None:
        current_tab = moose.Table(data_hierarchy.path + '/' + comp_path + '_Iex')
        moose.connect(current_tab, 'requestOut', pulsename, 'getOutputValue')
    else:
        current_tab = None
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

# Function that will allow for acquiring the distance of a compartment from the soma (credit: moose_nerp)
def get_dist_name(comp):
    name = comp.name
    xloc = comp.x
    yloc = comp.y
    zloc = comp.z
    dist = np.sqrt(xloc*xloc+yloc*yloc+zloc*zloc)
    return dist,name

# Function that will acquire the distance of a compartment from an origin selected by the user
def get_dist_name_from_soma(comp, origin):
    name = comp.name
    xloc = comp.x0 + (comp.x - comp.x0)/2 # x for center of compartment
    yloc = comp.y0 + (comp.y - comp.y0)/2 # y for center of compartment
    zloc = comp.z0 + (comp.z - comp.z0)/2 # z for center of compartment
    dist = np.sqrt((xloc-origin[0])**2+(yloc-origin[1])**2+(zloc-origin[2])**2)
    return dist,name

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

# Function that will set the parameters for specific compartments in a model neuron
def setSpecificCompParameters(comp,RM,CM,RA,initVm,E_leak):
    	SA = np.pi*comp.length*comp.diameter
	X_area = np.pi*comp.diameter*comp.diameter/4.0
	comp.Rm = RM/SA
	comp.Cm = CM*SA
	comp.Ra = RA*comp.length/X_area
	comp.initVm = initVm
	comp.Em = E_leak

# Function that will set up compartment with non uniform membrane resistance
def setSpecificCompParametersNonUniform(comp,origin,RM_soma,RM_end,RM_halfdist,RM_slope,CM,RA,initVm,E_leak):
    	SA = np.pi*comp.length*comp.diameter
	X_area = np.pi*comp.diameter*comp.diameter/4.0
	dist, _ = get_dist_name_from_soma(comp, origin)
	dist = dist*1e6 # convert length so that it is amenable to non-uniform equation
	RM_point=RM_end+(RM_soma-RM_end)/(1+np.exp((dist-RM_halfdist)/RM_slope))
	comp.Rm = RM_point/SA
	comp.Cm = CM*SA
	comp.Ra = RA*comp.length/X_area
	comp.initVm = initVm
	comp.Em = E_leak

# Function that will set up compartments with non uniform conductance
def setNonUniformConductance(comp_list, cell_path, libraryName, condSet, minq, maxq, qhalfdist,
			     qslope, origin, sag):
    Hpoints = []
    Gbars = []
    for comp in comp_list:
        comp = moose.element(cell_path + '/' + comp)
        dist, _ = get_dist_name_from_soma(comp, origin)
        for chan_name, cond in condSet.items():
            SA = np.pi*comp.length*comp.diameter
            proto = moose.element(libraryName + '/' + chan_name)
            chan = moose.copy(proto, comp, chan_name)[0]
	    dist = dist*1e6
	    hpoint=minq+(maxq-minq)/(1+np.exp(-(dist-qhalfdist)/qslope))
            chan.Gbar = hpoint*SA*sag
            m = moose.connect(chan, 'channel', comp, 'channel')
	    Hpoints.append(hpoint)
	    Gbars.append(chan.Gbar)
    return Hpoints, Gbars

# Function that will scale Rm and Cm for a compartment by some scaling factor
def scaleCompRmCm(comp, scaling_factor):
    comp.Rm = comp.Rm/scaling_factor
    comp.Cm = comp.Cm*scaling_factor


# Function that will scale the conductance of a channel by some scaling factor
def scaleChannelCond(chan, scaling_factor):
    chan.Gbar = chan.Gbar*scaling_factor

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
ChannelSettings = namedtuple('ChannelSettings', 'Xpow Ypow Zpow Erev name Xparam Yparam Zparam chan_type')

CadepParamSettings = namedtuple('CadepParams', 'Kd power tau')

HCNParamSettings = namedtuple('HCNParams', 'alpha_0 a z Vhalf exp_temp sim_temp q10')

CaPoolSettings = namedtuple('CaPoolSettings', 'CaBasal CaThick CaTau BufCapacity caName')

# Function that will add a set of channels and their corresponding maximal conductances 
# to each compartment in a model
def addChannelSet(condSet, library_name, cell, comp_type):
    for comp in moose.wildcardFind(cell.path + '/' + '#[TYPE=' + comp_type + ']'):
        for chan_name, cond in condSet.items():
            SA = np.pi*comp.length*comp.diameter
            proto = moose.element(library_name + '/' + chan_name)
            chan = moose.copy(proto, comp, chan_name)[0]
            chan.Gbar = cond*SA
            m = moose.connect(chan, 'channel', comp, 'channel')
    return

# Function to connect calcium to the channel
def connect_cal2chan(chan_name, chan_type, cellname, calname, comp_type):
    for comp in moose.wildcardFind(cellname.path + '/' + '#[TYPE=' + comp_type + ']'):
        capool = moose.element(comp.path + '/' + calname)
        chan = moose.element(comp.path + '/' + chan_name)
        if chan_type == 'ca_permeable':
            m = moose.connect(chan, 'IkOut', capool, 'current')
        elif chan_type == 'ca_dependent':
            m = moose.connect(capool, 'concOut', chan, 'concen')
        else:
            print ('unknown calcium connection type')
    return

# Function to create the prototype for the Calcium Pool
def CaPoolProto(libraryName, CaPoolParams):
    if not moose.exists(libraryName):
        lib = moose.Neutral(libraryName)
    else:
        lib = moose.element(libraryName)
    poolproto = moose.CaConc(lib.path + '/' + CaPoolParams.caName)
    poolproto.CaBasal = CaPoolParams.CaBasal
    poolproto.ceiling = 1
    poolproto.floor = 0
    poolproto.thick = CaPoolParams.CaThick
    poolproto.tau = CaPoolParams.CaTau
    return poolproto

# Function to copy the Calcium Pool from the library to the compartments
def add_calcium(libraryName, cellname, CaPoolParams, comp_type):
    Fday = 96485.33289
    caproto = CaPoolProto(libraryName, CaPoolParams)
    for comp in moose.wildcardFind(cellname.path + '/' + '#[TYPE=' + comp_type + ']'):
        capool = moose.copy(caproto, comp, CaPoolParams.caName,1)
        capool.length = comp.length
        capool.diameter = comp.diameter
        SA = np.pi*capool.length*capool.diameter
        vol = SA*capool.thick
        capool.B = 1/(Fday*vol*2)/CaPoolParams.BufCapacity
    return
# Function that will create from the biophysics defined for a channel a MOOSE channel
def createChanProto(libraryName, channelParams, rateParams, CaParams = None):
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
    channel.Zpower = channelParams.Zpow
    
    # Define the activation gating kinetics if they exist
    if channel.Xpower > 0 and 'HCN' not in channelParams.name :
        xGate = moose.HHGate(channel.path + '/' + 'gateX')
        xGate.setupAlpha(channelParams.Xparam + rateParams)
	print channelParams.Xpow

    # Define custom activation gating kinetics if they exist (Hyperpolarized channel kinetics)
    if channel.Xpower > 0 and 'HCN' in channelParams.name :
	Fday = 9.648e4
	R = 8.315
	FbyRT= Fday/ \
               (R*(273.16 + channelParams.Xparam.sim_temp))
	channel.Xpower = channelParams.Xpow
        xGate = moose.HHGate(channel.path + '/' + 'gateX')
	v_array = np.linspace(rateParams[1], rateParams[2], rateParams[0]+1)
	xGate.min = rateParams[1]
	xGate.max = rateParams[2]
	# custom equation for HCN gating kinetics
	alpha = (channelParams.Xparam.alpha_0*np.exp(-channelParams.Xparam.a*channelParams.Xparam.z \
		*(v_array-channelParams.Xparam.Vhalf)*FbyRT))
	beta = (channelParams.Xparam.alpha_0*np.exp((1-channelParams.Xparam.a) \
	       *channelParams.Xparam.z*(v_array-channelParams.Xparam.Vhalf)*FbyRT))
	q10 = channelParams.Xparam.q10**((channelParams.Xparam.sim_temp-channelParams.Xparam.exp_temp)/10)
	inf_x = alpha/(alpha+beta)
	tau_x = 1/(q10*(alpha+beta))
	#tau_x = [tau_x if tau_x > 2e-3 else 2e-3 for tau_x in tau_x]
	tau_x = np.array(tau_x)

	xGate.tableA = inf_x /tau_x
	xGate.tableB = 1 / tau_x
	print xGate.tableA
	print xGate.tableB
	print channelParams.Xpow
    
    # Define the inactivation gating kinetics if they exist    
    if channel.Ypower > 0:
        yGate = moose.HHGate(channel.path + '/' + 'gateY')
        yGate.setupAlpha(channelParams.Yparam + rateParams)

    if channel.Zpower > 0:
        channel.Zpower = channelParams.Zpow
        ca_array = np.linspace(CaParams[0], CaParams[1], CaParams[2])
        zGate = moose.HHGate(channel.path + '/' + 'gateZ')
        zGate.min = CaParams[0]
        zGate.max = CaParams[1]
        # custom equation for the Ca dynamics
        caterm = (ca_array/channelParams.Zparam.Kd)
        caterm = caterm**channelParams.Zparam.power
        inf_z = caterm/(1+caterm)
        tau_z = channelParams.Zparam.tau*np.ones(len(ca_array))
        # How to fill the A and B tables with custom gating
        zGate.tableA = inf_z / tau_z #tableA is equivalent to alpha
        zGate.tableB = 1 / tau_z # tableB is equivalent to alpha+beta
        # Set this for z Gate channel to use concentration
        channel.useConcentration = True
    
    # Set the tick for the channel so that it is set appropriately when it
    # is copied to a channel from the library    
    channel.tick = -1
    
    return channel

# Function that will create channels using the createChanProto function and place them
# in a MOOSE library, which will be utilized in adding channels to a particular compartment
def createChanLib(libraryName, channelSet, rateParams, CaParams):
    # Create a library to store the channel prototypes
    if not moose.exists(libraryName):
        lib = moose.Neutral(libraryName)
    else:
        lib = moose.element(libraryName)

    # Add all the channels to the MOOSE library
    for params in channelSet.values():
        chan = createChanProto(libraryName, params, rateParams, CaParams)

# Function that will add one channel to a compartment, which can be called in a loop for all
# channels provided
def addOneChan(library_name, channelName, conductance, compName):
    SA = np.pi*compName.length*compName.diameter
    proto = moose.element(library_name + '/' + channelName)
    chan = moose.copy(proto, compName, channelName)[0]
    chan.Gbar = conductance*SA
    m = moose.connect(chan, 'channel', compName, 'channel')

# Function that will add a channel set to multiple compartments 
def addChanList(library_name, condSet, comp):
    for chan_name, cond in condSet.items():
        SA = np.pi*comp.length*comp.diameter
        proto = moose.element(library_name + '/' + chan_name)
        chan = moose.copy(proto, comp, chan_name)[0]
        chan.Gbar = cond*SA
        m = moose.connect(chan, 'channel', comp, 'channel')

# Function that will create a multi-compartment model in MOOSE from a .p or .swc file
def createMultiCompCell(file_name, container_name, library_name, comp_type, channelSet, condSet,
                        rateParams, CaParams = None, CaPoolParams = None,
			cell_RM = None, cell_CM = None, cell_RA = None, cell_initVm = None, 
			cell_Em = None, dist_dep = False):
    # Create the channel types and store them in a library to be used by each compartment
    # in the model
    createChanLib(library_name, channelSet, rateParams, CaParams)

    # Load in the model in question
    if file_name.endswith('.p'):
        cell = moose.loadModel(file_name, container_name)
        for comp in moose.wildcardFind(cell.path + '/' + '#[TYPE=' + comp_type + ']'):
            for chan_name, cond in condSet.items():
                SA = np.pi*comp.length*comp.diameter
                proto = moose.element(library_name + '/' + chan_name)
                chan = moose.copy(proto, comp, chan_name)[0]
                chan.Gbar = cond*SA
                m = moose.connect(chan, 'channel', comp, 'channel')
        # Add the calcium pool to each compartment in the cell if it has been specified
        if (CaPoolParams != None):
	    add_calcium(library_name, cell, CaPoolParams, comp_type)
	    for key in channelSet.keys():
	        if ("Ca" in key):
		    connect_cal2chan(channelSet[key].name, channelSet[key].chan_type, cell,
                    		     CaPoolParams.caName, comp_type)
    else:
        cell = moose.loadModel(file_name, container_name)
        setCompParameters(cell, comp_type, cell_RM, cell_CM, cell_RA, cell_initVm, cell_Em)
	if (dist_dep == False):
            for comp in moose.wildcardFind(cell.path + '/' + '#[TYPE=' + comp_type + ']'):
                for chan_name, cond in condSet.items():
                    SA = np.pi*comp.length*comp.diameter
                    proto = moose.element(library_name + '/' + chan_name)
                    chan = moose.copy(proto, comp, chan_name)[0]
                    chan.Gbar = cond*SA
                    m = moose.connect(chan, 'channel', comp, 'channel')
            # Add the calcium pool to each compartment in the cell if it has been specified
            if (CaPoolParams != None):
	        add_calcium(library_name, cell, CaPoolParams, comp_type)
	        for key in channelSet.keys():
	            if ("Ca" in key):
		        connect_cal2chan(channelSet[key].name, channelSet[key].chan_type, cell,
                    		         CaPoolParams.caName, comp_type)
	# Work in progress code for distance dependent conductance
	if (dist_dep == True):
            for comp in moose.wildcardFind(cell.path + '/' + '#[TYPE=' + comp_type + ']'):
                for chan_name, cond in condSet.items():
		    if 'soma' in comp.name:
		        conductance = cond['soma']
		    elif 'apical' in comp.name:
		        conductance = cond['apical']
		    else:
			conductance = cond['dend'] 
                    SA = np.pi*comp.length*comp.diameter
                    proto = moose.element(library_name + '/' + chan_name)
                    chan = moose.copy(proto, comp, chan_name)[0]
                    chan.Gbar = conductance*SA
                    m = moose.connect(chan, 'channel', comp, 'channel')
            # Add the calcium pool to each compartment in the cell if it has been specified
            if (CaPoolParams != None):
	        add_calcium(library_name, cell, CaPoolParams, comp_type)
	        for key in channelSet.keys():
	            if ("Ca" in key):
		        connect_cal2chan(channelSet[key].name, channelSet[key].chan_type, cell,
                    		         CaPoolParams.caName, comp_type)	
    return cell

# Function that will add excitatory or inhibitory synchans (with a SimpleSynHandler) to a single or
# set of compartments
def addSynChan(compname, synparams):
    synchan = moose.SynChan(compname.path + '/' + synparams['name'])
    synchan.Gbar = synparams['Gbar']
    synchan.tau1 = synparams['tau1']
    synchan.tau2 = synparams['tau2']
    synchan.Ek = synparams['erev']
    msg = moose.connect(compname, 'channel', synchan, 'channel')
    sh = moose.SimpleSynHandler(synchan.path + '/' + synparams['name'])
    moose.connect(sh, 'activationOut', synchan, 'activation')
    return sh

# Function that will create a SpikeGen
def createSpikeGen(spikeParams):
    spikegen = moose.SpikeGen(spikeParams['name'])
    spikegen.threshold = spikeParams['threshold']
    spikegen.refractT = spikeParams['refractT']
    return spikegen


# Function that will produce a RandSpike, which can be used as a substitute for a pre-synaptic neuron,
# and connect it to a channel on the post-synaptic compartment 
def createRandSpike(spikeParams, synHandler):
    pre_syn = moose.RandSpike(spikeParams['name'])
    pre_syn.rate = spikeParams['rate']
    pre_syn.refractT = spikeParams['refractT']
    index = synHandler.synapse.num
    synHandler.synapse.num = synHandler.synapse.num + 1
    synHandler.synapse[index].delay = spikeParams['delay']
    moose.connect(pre_syn, 'spikeOut', synHandler.synapse[index], 'addSpike')
    return pre_syn


# Function that will allow for the use of the hsolve implicit numerical method (credit: moose_nerp)
def hsolve(soma, simdt):
    hsolve = moose.HSolve(soma.parent.path + '/hsolve')
    hsolve.dt = simdt
    hsolve.target = soma.path


# some neuron utilities for plotting synapses and networks
def attach_current_clamp(cell, delay=5, dur=1, amp=.1, loc=1):
    """Attach a current Clamp to a cell.

    :param cell: Cell object to attach the current clamp.
    :param delay: Onset of the injected current.
    :param dur: Duration of the stimulus.
    :param amp: Magnitude of the current.
    :param loc: Location on the dendrite where the stimulus is placed.
    """
    stim = h.IClamp(cell.dend(loc))
    stim.delay = delay
    stim.dur = dur
    stim.amp = amp
    return stim

def set_recording_vectors(cell):
    """Set soma, dendrite, and time recording vectors on the cell.

    :param cell: Cell to record from.
    :return: the soma, dendrite, and time vectors as a tuple.
    """
    soma_v_vec = h.Vector()   # Membrane potential vector at soma
    dend_v_vec = h.Vector()   # Membrane potential vector at dendrite
    t_vec = h.Vector()        # Time stamp vector
    soma_v_vec.record(cell.soma(0.5)._ref_v)
    dend_v_vec.record(cell.dend(0.5)._ref_v)
    t_vec.record(h._ref_t)
    return soma_v_vec, dend_v_vec, t_vec

def simulate(tstop=25):
    """Initialize and run a simulation.

    :param tstop: Duration of the simulation.
    """
    h.tstop = tstop
    h.run()

def show_output(soma_v_vec, dend_v_vec, t_vec, new_fig=True):
    """Draw the output.

    :param soma_v_vec: Membrane potential vector at the soma.
    :param dend_v_vec: Membrane potential vector at the dendrite.
    :param t_vec: Timestamp vector.
    :param new_fig: Flag to create a new figure (and not draw on top
            of previous results)
    """
    if new_fig:
        plt.figure(figsize=(8,4)) # Default figsize is (8,6)
    soma_plot = plt.plot(t_vec, soma_v_vec, color='black')
    dend_plot = plt.plot(t_vec, dend_v_vec, color='red')
    plt.legend(soma_plot + dend_plot, ['soma', 'dend(0.5)'])
    plt.xlabel('time (ms)')
    plt.ylabel('mV')
