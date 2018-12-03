from neuron import h,gui
h.load_file('ri06_figs6-7_run.hoc')
import numpy as np
for sec in h.allsec():
    gbar = np.pi*sec.gbar_h*(sec.L*sec.diam/2*sec.nseg)*1e-12
    print sec.name(), gbar
