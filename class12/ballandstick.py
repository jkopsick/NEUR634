from neuron import h

class BallAndStick(object):
    """Two-section cell: A soma with active channels and
    a dendrite with passive properties."""
    def __init__(self, bs_properties):
        self.create_sections()
        self.build_topology()
        self.build_subsets()
	geometry = bs_properties.get('cell_geometry')
	soma_d, soma_l, dend_d, dend_l, dend_nseg = geometry
        self.define_geometry(soma_d, soma_l, dend_d, dend_l, dend_nseg)
	biophysics = bs_properties.get('cell_biophysics')
	Ra, Cm, gna_max, gk_max, gl_max, e_leak, g_pas, e_pas = biophysics
        self.define_biophysics(Ra, Cm, gna_max, gk_max, gl_max, e_leak,
			       g_pas, e_pas)
    #
    def create_sections(self):
        """Create the sections of the cell."""
        # NOTE: cell=self is required to tell NEURON of this object.
        self.soma = h.Section(name='soma', cell=self)
        self.dend = h.Section(name='dend', cell=self)
    #
    def build_topology(self):
        """Connect the sections of the cell to build a tree."""
        self.dend.connect(self.soma(1))
    #
    def define_geometry(self, soma_d, soma_l, dend_d, dend_l, dend_nseg):
        """Set the 3D geometry of the cell."""
        self.soma.L = soma_l # microns
	self.soma.diam = soma_d # microns
        self.dend.L = dend_l                      # microns
        self.dend.diam = dend_d                     # microns
        self.dend.nseg = dend_nseg # number of divisions for dendrite
        h.define_shape() # Translate into 3D points.
    #
    def define_biophysics(self, Ra, Cm, gna_max, gk_max, gl_max, e_leak,
			  g_pas, e_pas):
        """Assign the membrane properties across the cell."""
        for sec in self.all: # 'all' defined in build_subsets
            sec.Ra = Ra    # Axial resistance in Ohm * cm
            sec.cm = Cm      # Membrane capacitance in micro Farads / cm^2
        # Insert active Hodgkin-Huxley current in the soma
        self.soma.insert('hh')
        for seg in self.soma:
            seg.hh.gnabar = gna_max  # Sodium conductance in S/cm2
            seg.hh.gkbar = gk_max  # Potassium conductance in S/cm2
            seg.hh.gl = gl_max    # Leak conductance in S/cm2
            seg.hh.el = e_leak     # Reversal potential in mV
        # Insert passive current in the dendrite
        self.dend.insert('pas')
        for seg in self.dend:
            seg.pas.g = g_pas  # Passive conductance in S/cm2
            seg.pas.e = e_pas    # Leak reversal potential mV
    #
    def build_subsets(self):
        """Build subset lists. For now we define 'all'."""
        self.all = h.SectionList()
        self.all.wholetree(sec=self.soma)
