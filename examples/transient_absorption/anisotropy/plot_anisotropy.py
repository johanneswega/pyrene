from pyrene.transient_absorption import Kinetics, Spectra

delays = [0.5, 1, 2, 5, 10, 20, 50, 100, 200]
s = Spectra(files=['r.npy'],
          x_cuts=[(400, 780)],  
          delay=delays,
          alphas=[0.5 for _ in delays],   
          overview=True,
          outside=True,
          ylim=[-0.6, 0.6],       
          steady_state=[['abs.txt', (400, 780), -0.5, 'b', 'Abs.'], 
                        ['exported/MA/spectra/1_ps.txt', (400, 780), +0.5, 'gray', r'$\Delta A_\text{MA}$ at 1 ps']])
s.ylabel = r'$r$'
s.show()

k = Kinetics(files=['r.npy'], 
             wavelength=[570],
             figsize=[6, 3.5],
             outside=True,
             marker=['o'],
             markersize=[3],
             alphas=[0.3], 
             savefig='kinetics_fit.png',
             x_cuts=[(0.3, 500)],
             xscale='log')
k.ylabel = r'$r$'
# to fit function --> choose fit function from fit_functions module 
# (src/pyrene/fitting/fit_functions.py)
from pyrene.fitting.fit_functions import mono_exp, bi_exp, mono_exp_with_bg
# then load fit function and file index
k.fit(file_index=0, model=bi_exp, p0=[0.1, 1, 0.4, 100])
#k.fit(file_index=0, model=mono_exp_with_bg, p0=[0.4, 100, 0.1])
k.show()