from pyrene.electrochemistry import CV
from pyrene.standard.misc import rainbow
import numpy as np 
import os

### get files ###
# list of scan rates in mV/s
v = np.array([100, 305, 494, 729, 1000, 2580, 3090, 4890])
files = np.sort(['data/scan_rate_dependence_2/' + i for i in os.listdir('data/scan_rate_dependence_2') if '.csv' in i])
colors = rainbow(files)
labels = ['%s mV'%i for i in v]
Fc = [0.38 for _ in files]

c = CV(files=files, colors=colors, labels=labels, Fc=Fc, outside=True, figsize=[6, 4])
# get slope of peak irreversible forward peak potential vs. log(scan_rate)
# this can give mechanistic information on the irreversible chemical reaction following the reduction/oxidation
# if slope around 30 mV/s --> just reaction e.g. R- --> X
# if slope around 20 mV/s --> dimerization e.g. R- + R- --> R2-
c.EC_pot_shift(scan_rates=v/100, E_range=[-3, 2])
c.show()