from pyrene.transient_absorption import GlobalAnalysis

# suppose we want to model as branching scheme like this one
#                       A
#                [k0]|     | [k1]
#                    v     v
#                    B     C
#                [k2]|     | [k3]
#                    v     v

# the rate equations for the species are 
# dA/dt = -(k0 + k1) * A    + 0 * B         + 0 * C
# dB/dt = +k0 * A           -k2 * B         + 0 * C
# dB/dt = +k1 * A           + 0 * B         - k3 * C

# which can be written in matrix form
#       [A]    [-(k0 + k1),     0,      0] [A]
# d/dt  [B] =  [+k0,          -k2,      0] [B]
#       [C]    [+k1,            0,    -k3] [C]

# for the target analysis we need to provide this K-matrix 
# as a numpy array to the class
import numpy as np
K = np.array([['-k0 -k1', '0', '0'],
              ['+k0', '-k2', '0'],
              ['+k1', '0', '-k3']])

# as well as the starting concentrations at t=0
# assuming everything is in A at t=0
C0 = [1, 0, 0]

g = GlobalAnalysis(files=['data/dA.npy'], 
                   error=['data/Error.npy'],                    # error matrix (optional but encouraged)
                   model='target',
                   C0=[C0],
                   K=[K],
                   wavelength=[(np.arange(350, 600, 40))],      # wavelengths to plot fitted kinetic at
                   scale=[(-30, 30)],
                   p0=[(1, 10, 20, 100)],                           # guess parameters of the fit (tau's)
                   y_cuts=[(0.3, 1800)], x_cuts=[(330, 720)],
                   scatter=[(580-30, 580+40)],
                   yscale='log')
g.fit()