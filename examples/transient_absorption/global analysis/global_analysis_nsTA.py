from pyrene.transient_absorption import GlobalAnalysis
import numpy as np

g = GlobalAnalysis(files=['data/nano.npy'], 
                   experiment='nano',
                   wavelength=[(np.arange(350, 600, 40))],
                   scale=[(-30, 30)],
                   p0=[(25, 1e3, 10e3)], 
                   y_cuts=[(2, 500e3)], x_cuts=[(380, 720)],
                   yscale='log')
g.fit()