from pyrene.transient_absorption import GlobalAnalysis
import numpy as np

g = GlobalAnalysis(files=['r.npy'],                   # error matrix (optional but encouraged)
                   wavelength=[(550, 560, 570, 580, 590)],      # wavelengths to plot fitted kinetic at
                   scale=[(-0.6, 0.6)],
                   p0=[[1, 20]],                           # guess parameters of the fit (tau's)
                   y_cuts=[(0.3, 1800)], x_cuts=[(550, 600)],
                   yscale='log')
g.fit()