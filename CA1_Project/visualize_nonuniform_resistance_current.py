# These gating kinetics look much better now, but will need to make sure that these are converted
# properly into SI units, or since the mod file is defined in physiological units, would need
# to convert the rest of the model into physiological units (not what I really want to do though)

import numpy as np
import pylab as plt
plt.ion()

rmsoma2=6.29013
rmend2=3.1916
rmhalfdis2= 100.05
rmsteep2= 50.48
dis2 = np.linspace(0,533,5001)
rmpoint2=rmend2+(rmsoma2-rmend2)/(1+np.exp((dis2-rmhalfdis2)/rmsteep2))

plt.figure()
plt.plot(dis2,rmpoint2)
plt.xlabel('Distance from Soma (um)')
plt.ylabel('Membrane Resistance (ohm/m^2)')

minq=2.1582  		# units are pS/um2
maxq=1.5969  		# units are pS/um2
qhalfdis=98.753
qsteep=50.07
hpoint=minq+(maxq-minq)/(1+np.exp(-(dis2-qhalfdis)/qsteep))

plt.figure()
plt.plot(dis2,hpoint)
plt.xlabel('Distance from Soma (um)')
plt.ylabel('Ih Conductance (S/m^2)')

rmsoma3= 7.4697
rmend3=1.0216
rmhalfdis3= 115.07
rmsteep3= 20
dis3 = np.linspace(0,611,5001)
rmpoint3=rmend3+(rmsoma3-rmend3)/(1+np.exp((dis3-rmhalfdis3)/rmsteep3))

plt.figure()
plt.plot(dis3,rmpoint3)
plt.xlabel('Distance from Soma (um)')
plt.ylabel('Membrane Resistance (ohm/m^2)')

minq=0.1002  		# units are pS/um2
maxq=14.349  		# units are pS/um2
qhalfdis=216.65
qsteep=79.4
hpoint=minq+(maxq-minq)/(1+np.exp(-(dis3-qhalfdis)/qsteep))

plt.figure()
plt.plot(dis3,hpoint)
plt.xlabel('Distance from Soma (um)')
plt.ylabel('Ih Conductance (S/m^2)')
