# This python file contains the necessary code to create all the channels necessary for
# a compartment (in this case, a patch of the HH giant squid axon) and to load these channels
# into a library in MOOSE

import moose
import util as u
from collections import namedtuple

# Set the resting membrane potential to be used in creating the channels
# and defining their biophysics
EREST_ACT = -70e-3 #: Resting membrane potential


#: Create tuples for each kinetic type for the channels used and then set
# the tuples to the kinetics desired 

Na_m_params = u.AlphaBetaChanParams(1e5 * (25e-3 + EREST_ACT),   # 'A_A':
                -1e5,                       # 'A_B':
                -1.0,                       # 'A_C':
                -25e-3 - EREST_ACT,         # 'A_D':
               -10e-3,                      # 'A_F':
               4e3,                     # 'B_A':
                0.0,                        # 'B_B':
                0.0,                        # 'B_C':
                0.0 - EREST_ACT,            # 'B_D':
                18e-3                       # 'B_F':
                )

Na_h_params = u.AlphaBetaChanParams(70.0,                        # 'A_A':
                0.0,                       # 'A_B':
                0.0,                       # 'A_C':
                0.0 - EREST_ACT,           # 'A_D':
                0.02,                     # 'A_F':
                1000.0,                       # 'B_A':
                0.0,                       # 'B_B':
                1.0,                       # 'B_C':
                -30e-3 - EREST_ACT,        # 'B_D':
                -0.01                    # 'B_F':
                )

K_n_params = u.AlphaBetaChanParams(1e4 * (10e-3 + EREST_ACT),   #  'A_A':
               -1e4,                      #  'A_B':
               -1.0,                       #  'A_C':
               -10e-3 - EREST_ACT,         #  'A_D':
               -10e-3,                     #  'A_F':
               0.125e3,                   #  'B_A':
               0.0,                        #  'B_B':
               0.0,                        #  'B_C':
               0.0 - EREST_ACT,            #  'B_D':
               80e-3                       #  'B_F':  
               )

HCN_n_params = u.HCNParamSettings(0.00057, 0.4, 7, -81e-3, 9.648e4, 8.315, 35)


Na_param = u.ChannelSettings(Xpow = 3, Ypow = 1, Zpow = 0, 
			     Erev = 115e-3 + EREST_ACT, 
			     name = 'Na', Xparam = Na_m_params, 
                             Yparam = Na_h_params, Zparam = [],
                             chan_type = [])

K_param = u.ChannelSettings(Xpow = 4, Ypow = 0, Zpow = 0, Erev = -12e-3 + EREST_ACT, name = 'K',
                            Xparam = K_n_params, Yparam = [], Zparam = [],
                            chan_type = [])

HCN_param = u.ChannelSettings(Xpow = 1, Ypow = 0, Zpow = 0, Erev = -25e-3, name = 'HCN',
                            Xparam = HCN_n_params, Yparam = [], Zparam = [],
                            chan_type = [])


#: We define the rate parameters, which are functions of Vm as
#: interpolation tables looked up by membrane potential.
#: Minimum x-value for the interpolation table
VMIN = -30e-3 + EREST_ACT
#: Maximum x-value for the interpolation table
VMAX = 120e-3 + EREST_ACT
#: Number of divisions in the interpolation table
VDIVS = 3400

# q10 factor for HCN channel
q10 = 4.5

# Creation of the channels and placing them into a specified MOOSE library
chan_set = {'Na': Na_param, 'K': K_param, 'HCN' : HCN_param}
rateParams = (VDIVS, VMIN, VMAX)
HCNParams = (-100e-3, 50e-3, VDIVS+1, q10)
