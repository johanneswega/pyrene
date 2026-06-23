from pyrene.transient_absorption import GlobalAnalysis
import numpy as np

g = GlobalAnalysis(files=['data/dA.npy'], 
                   error=['data/Error.npy'],                    # error matrix (optional but encouraged)
                   wavelength=[(350, 400, 450, 500, 550)],      # wavelengths to plot fitted kinetic at
                   scale=[(-30, 30)],
                   p0=[(1, 10, 100)],                           # guess parameters of the fit (tau's)
                   y_cuts=[(0.3, 1800)], x_cuts=[(330, 720)],
                   scatter=[(580-30, 580+40)],
                   yscale='log')
g.fit()