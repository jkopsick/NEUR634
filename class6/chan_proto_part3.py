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

# activation constants for alphas and betas (obtained by
# matching m2 to Tkatch et al., 2000 Figs 2c, and mtau to fig 2b)
qfactKaF = 2
KaF_X_params = u.AlphaBetaChanParams(A_rate = 1.8e3*qfactKaF,
                                      A_B = 0,
                                      A_C = 1.0,
                                      A_vhalf = 18e-3,
                                      A_vslope = -13.0e-3,
                                      B_rate = 0.45e3*qfactKaF,
                                      B_B = 0.0,
                                      B_C = 1.0,
                                      B_vhalf = -2.0e-3,
                                      B_vslope = 11.0e-3)

#inactivation consts for alphas and betas obtained by matching Tkatch et al., 2000 Fig 3b,
#and tau voltage dependence consistent with their value for V=0 in fig 3c.
#slowing down inact improves spike shape tremendously
KaF_Y_params = u.AlphaBetaChanParams(A_rate = 0.105e3*qfactKaF,
                                      A_B = 0,
                                      A_C = 1.0,
                                      A_vhalf = 121e-3,
                                      A_vslope = 22.0e-3,
                                      B_rate = 0.065e3*qfactKaF,
                                      B_B = 0.0,
                                      B_C = 1.0,
                                      B_vhalf = 55.0e-3,
                                      B_vslope = -11.0e-3)

# Create the channel settings associated with each channel type to be used in the model
Na_param = u.ChannelSettings(Xpow = 3, Ypow = 1, Erev = 115e-3 + EREST_ACT, name = 'Na',
                           Xparam = Na_m_params, Yparam = Na_h_params)

K_param = u.ChannelSettings(Xpow = 4, Ypow = 0, Erev = -12e-3 + EREST_ACT, name = 'K',
                          Xparam = K_n_params, Yparam = [])

KaF_param = u.ChannelSettings(Xpow = 2, Ypow = 1, Erev = -20e-3 + EREST_ACT, name = 'KaF',
                              Xparam = KaF_X_params, Yparam = KaF_Y_params)


#: We define the rate parameters, which are functions of Vm as
#: interpolation tables looked up by membrane potential.
#: Minimum x-value for the interpolation table
VMIN = -30e-3 + EREST_ACT
#: Maximum x-value for the interpolation table
VMAX = 120e-3 + EREST_ACT
#: Number of divisions in the interpolation table
VDIVS = 3000

# Creation of the channels and placing them into a specified MOOSE library
chan_set = {'Na': Na_param, 'K': K_param, 'KaF': KaF_param}
rateParams = (VDIVS, VMIN, VMAX)
