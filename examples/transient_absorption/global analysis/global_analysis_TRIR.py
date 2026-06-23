from pyrene.transient_absorption import GlobalAnalysis
import numpy as np

g = GlobalAnalysis(files=['data/data.pdat'], 
                   IR=True,
                   wavelength=[(np.arange(2080, 2280, 40))],      # wavelengths to plot fitted kinetic at
                   scale=[(-10, 10)],
                   p0=[(1, 10, 100)],                           # guess parameters of the fit (tau's)
                   y_cuts=[(0.3, 1800)],
                   yscale='log')
g.fit()