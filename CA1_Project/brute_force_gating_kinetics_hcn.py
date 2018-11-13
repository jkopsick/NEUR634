# These gating kinetics look much better now, but will need to make sure that these are converted
# properly into SI units, or since the mod file is defined in physiological units, would need
# to convert the rest of the model into physiological units (not what I really want to do though)

import numpy as np
import pylab as plt
plt.ion()
a0 = 0.00057
ab = 0.4
vhalf = -81
zeta = 7
temp = 33
gas = 8.315
farad = 9.648e4
celsius = 35
qten = 4.5

varray = np.linspace(-100,50,5001)
alpha = a0*np.exp((0.001)*(-ab)*zeta*(varray-vhalf)*farad/(gas*(273.16+celsius)))
beta = a0*np.exp((0.001)*(1-ab)*zeta*(varray-vhalf)*farad/(gas*(273.16+celsius)))
q10 = qten**((celsius-temp)/10)
a = q10*alpha
b = q10*beta
inf = a/(a+b)
tau = 1/(a+b)

plt.figure()
plt.plot(varray,inf)
plt.figure()
plt.plot(varray,tau)
