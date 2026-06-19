from pyrene.transient_absorption import Kinetics

k = Kinetics(files=['data/dA1.npy'], 
             wavelength=[350],
             marker=['o'],
             markersize=[3],
             alphas=[0.3], 
             x_cuts=[(0.3, 1800)],
             xscale='log')

# to fit function --> choose fit function from fit_functions module 
# (src/pyrene/fitting/fit_functions.py)
from pyrene.fitting.fit_functions import mono_exp, mono_exp_with_bg, bi_exp_with_bg
# then load fit function and file index
k.fit(file_index=0, model=bi_exp_with_bg, p0=[1, 1, 10, 100, 3])
k.show()