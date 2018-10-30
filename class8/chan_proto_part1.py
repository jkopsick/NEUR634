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

CaL_X_params = u.AlphaBetaChanParams(A_rate = -880, A_B = -220e3, A_C = -1.0,
                                     A_vhalf = 4.0003e-3, A_vslope = -7.5e-3,
                                     B_rate = -284, B_B = 71e3, B_C = -1.0,
                                     B_vhalf = -4.0003e-3, B_vslope = 5e-3)

SK_Z_params = u.CadepParamSettings(Kd = 0.57e-3, power = 5.2, tau = 4.9e-3)

Ca_pool_params = u.CaPoolSettings(CaBasal = 50e-6, CaThick = 1e-6, 
                                  CaTau = 20e-3, BufCapacity = 20,
                                  caName = 'CaPool')

Na_param = u.ChannelSettings(Xpow = 3, Ypow = 1, Zpow = 0, 
			     Erev = 115e-3 + EREST_ACT, 
			     name = 'Na', Xparam = Na_m_params, 
                             Yparam = Na_h_params, Zparam = [],
                             chan_type = [])

K_param = u.ChannelSettings(Xpow = 4, Ypow = 0, Zpow = 0, Erev = -12e-3 + EREST_ACT, name = 'K',
                            Xparam = K_n_params, Yparam = [], Zparam = [],
                            chan_type = [])

KaF_param = u.ChannelSettings(Xpow = 2, Ypow = 1, Zpow = 0, Erev = -20e-3 + EREST_ACT, name = 'KaF',
                              Xparam = KaF_X_params, Yparam = KaF_Y_params, Zparam = [],
			      chan_type = [])

sk_params = u.ChannelSettings(Xpow = 0, Ypow = 0, Zpow = 1, Erev = -87e-3, name = 'SKCa', Xparam = [], 
			      Yparam = [], Zparam = SK_Z_params, chan_type = 'ca_dependent')

CaLparam = u.ChannelSettings(Xpow = 1, Ypow = 0, Zpow = 0, Erev = 130e-3,
                             Xparam = CaL_X_params, Yparam = [], Zparam = [],
                             name = 'CaL', chan_type = 'ca_permeable')


#: We define the rate parameters, which are functions of Vm as
#: interpolation tables looked up by membrane potential.
#: Minimum x-value for the interpolation table
VMIN = -30e-3 + EREST_ACT
#: Maximum x-value for the interpolation table
VMAX = 120e-3 + EREST_ACT
#: Number of divisions in the interpolation table
VDIVS = 3000

CaMIN = 0

CaMAX = 1

CaDIVS = 10000

# Creation of the channels and placing them into a specified MOOSE library
chan_set = {'Na': Na_param, 'K': K_param, 'KaF': KaF_param, 'SKCa' : sk_params, 'CaL' : CaLparam}
rateParams = (VDIVS, VMIN, VMAX)
CaParams = (CaMIN, CaMAX, CaDIVS)
