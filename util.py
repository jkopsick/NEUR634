import moose
import numpy as np

# Create an object to store the simple hierarchy, and the compartment object soma -- must use characters around compartment name
def createCompartment(directory,compname,length,radius,RM,CM,RA,Em):
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
    compname.initVm = Em
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

def createDataTables(compname,data_hierarchy):
    # Create a new path using string manipulation to be used in the creation
    # of a unique data table for each compartment
    comp_path = compname.path.split('/')[-1]
    comp_path = comp_path.strip(']')
    comp_path = comp_path.replace('[','')
    # Create the unique table and connect this to the compartment
    Vmtab = moose.Table(data_hierarchy.path + '/' + comp_path + '_Vm')
    moose.connect(Vmtab, 'requestOut', compname, 'getVm')
    return Vmtab

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
    
def setCompParameters(compvector,comptype,RM,CM,RA,E_leak):
    # Loop through all the compartments in the vector
    for comp in moose.wildcardFind(compvector.path + '/' + '#[TYPE=' + comptype + ']'):
    	SA = comp.length*comp.diameter
	X_area = np.pi*comp.diameter*comp.diameter/4.0
	comp.Rm = RM/SA
	comp.Cm = CM*SA
	comp.Ra = RA*comp.length/X_area
	comp.initVm = E_leak
	comp.Em = E_leak    
    
