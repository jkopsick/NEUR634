from neuron import h,gui
h.load_file('ri04_figs6-7_run.hoc')
import numpy as np
for sec in h.allsec():
    gbar = np.pi*sec.gbar_h*sec.L*sec.diam*1e-12
    print sec.name(), gbar
