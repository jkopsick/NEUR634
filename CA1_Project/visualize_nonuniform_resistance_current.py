# These gating kinetics look much better now, but will need to make sure that these are converted
# properly into SI units, or since the mod file is defined in physiological units, would need
# to convert the rest of the model into physiological units (not what I really want to do though)

import numpy as np
import pylab as plt
plt.ion()

rmsoma=0.629013*100000
rmend=0.31916*100000
rmhalfdis=0.10005*1000
rmsteep=0.5048*100
dis = np.linspace(0,533,5001)
rmpoint=rmend+(rmsoma-rmend)/(1+np.exp((dis-rmhalfdis)/rmsteep))
#{ g_pas=1/(rmpoint/scale_spines) cm=Cm*scale_spines Ra=global_ra }

plt.figure()
plt.plot(dis,rmpoint)


rmsoma2=6.29013
rmend2=3.1916
rmhalfdis2= 100.05
rmsteep2= 50.48
dis2 = np.linspace(0,533,5001)
rmpoint2=rmend2+(rmsoma2-rmend2)/(1+np.exp((dis2-rmhalfdis2)/rmsteep2))

plt.figure()
plt.plot(dis2,rmpoint2)

minq=2.1582  		# units are pS/um2
maxq=1.5969  		# units are pS/um2
qhalfdis=98.753
qsteep=50.07
hpoint=minq+(maxq-minq)/(1+np.exp(-(dis-qhalfdis)/qsteep))

plt.figure()
plt.plot(dis,hpoint)