from pyrene.electrochemistry import CV
from pyrene.standard.misc import rainbow
import numpy as np 
import os

### get files ###
# list of scan rates in mV/s
v = np.array([10, 63, 161, 305, 494, 729, 1010, 1330, 1710, 2120, 2580, 3090])
files = np.sort(['data/scan_rate_dependence_1/' + i for i in os.listdir('data/scan_rate_dependence_1') if '.csv' in i])
colors = rainbow(files)
labels = ['%s mV'%i for i in v]
Fc = [0.60508 for _ in files]

c = CV(files=files, colors=colors, labels=labels, Fc=Fc, outside=True, figsize=[6, 4], savefig='scanrates.png')

# peak separation as function of scan rate
#c.peak_separation(scan_rates=v/1000, E_range=[-1.2, -0.8])

# perform a cotroll fit (according to Macedo et al.: https://pubs.acs.org/doi/10.1021/acs.analchem.3c04181)
#c.cottrell_fit(scan_rates=v/1000, E_range=[-3, 2], plot=True)

### Randles-Ševčík fit ###
# WE electrode glassy carbon d = 3 mm = 0.3 cm | CH Instruments 
d = 0.3
# electrode area in cm2
A = np.pi*(d/2)**2
# concentration of the substrate in mM
conc = 0.2
c.randles_sevcik(scan_rates=v/1000, E_range=[-3, 2], conc=conc, A=A)

c.show()