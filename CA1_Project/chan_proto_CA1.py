# This python file contains the necessary code to create all the channels necessary for
# a compartment and to load these channels into a library in MOOSE

import moose
import util as u
from collections import namedtuple

#: Create tuples for each kinetic type for the channels used and then set
# the tuples to the kinetics desired 

# naF params taken from CA1 moose_nerp
Na_m_params= u.AlphaBetaChanParams(
    A_rate=-12.0052e-3*2.14e6,
    A_B=-0.4*2.14e6,
    A_C=-1,
    A_vhalf=30.013e-3,
    A_vslope=-7.2e-3,
    B_rate=3.722e-3*2.14e6,
    B_B=0.124*2.14e6,
    B_C=-1,
    B_vhalf=30.013e-3,
    B_vslope=7.2e-3)



Na_h_params= u.AlphaBetaChanParams(
    A_rate=(45.013e-3+15.0e-3)*0.03*2.14e6,
    A_B=0.03*2.14e6,
    A_C=-1,
    A_vhalf=45.013e-3+15e-3,
    A_vslope=3.5e-3,
    B_rate=-(45.013e-3+15e-3)*0.01*2.14e6,
    B_B=-0.01*2.14e6,
    B_C=-1,
    B_vhalf=45.013e-3+15e-3,
    B_vslope=-3.5e-3)


Na_param = u.ChannelSettings(Xpow = 3, Ypow = 1, Zpow = 0, 
			     Erev = 50e-3, 
			     name = 'Na', Xparam = Na_m_params, 
                             Yparam = Na_h_params, Zparam = [],
                             chan_type = [])

# Delayed rectifier K params taken from CA1 moose_nerp
K_n_params= u.AlphaBetaChanParams(
    A_rate=6.5,
    A_B=0,
    A_C=0,
    A_vhalf=0,
    A_vslope=-12.5e-3,
    B_rate=24,
    B_B=0,
    B_C=0,
    B_vhalf=0,
    B_vslope=33.5e-3)

K_param = u.ChannelSettings(Xpow = 1, Ypow = 0, Zpow = 0, Erev = -90e-3, name = 'K',
			      Xparam = K_n_params, Yparam = [], Zparam = [], 
			      chan_type = [])

# HCN params taken from Golding et. al 2005 
HCN_n_params = u.HCNParamSettings(alpha_0=0.00057*1e3, a=0.4, z=7, Vhalf=-81e-3,
				  exp_temp = 33, sim_temp = 35, q10 = 4.5)

HCN_param = u.ChannelSettings(Xpow = 1, Ypow = 0, Zpow = 0, Erev = -25e-3, name = 'HCN',
                              Xparam = HCN_n_params, Yparam = [], Zparam = [],
                              chan_type = [])


#: We define the rate parameters, which are functions of Vm as
#: interpolation tables looked up by membrane potential.
#: Minimum x-value for the interpolation table
VMIN = -120e-3
#: Maximum x-value for the interpolation table
VMAX = 50e-3
#: Number of divisions in the interpolation table
VDIVS = 3400

# Creation of the channels and placing them into a specified MOOSE library
chan_set = {'Na': Na_param, 'K': K_param, 'HCN' : HCN_param}

# rateParams defined this way as this allows for the list of parameters to be loaded in appropriately
# via setupAlpha MOOSE function
rateParams = (VDIVS, VMIN, VMAX)
