from pyrene.steady_state import Absorption
from pyrene.standard.misc import rainbow
import os 
import numpy as np

files = ['Data/conc/' + i for i in np.sort(os.listdir('Data/conc')) if not '.D' in i]
labels = ["%s"%i for i in range(len(files))]
colors = rainbow(files)

a = Absorption(files=files,
               labels=labels,
               colors=colors,
                x_cuts=[(300, 500) for _ in files],
                wn=False)

a.epsilon_regression(M=178.23, m=np.array([0.26, 0.48, 0.67, 0.95]),
                     Vstock=1, Vcuv=4, Vadd=np.array([100, 50, 50, 50]), l=1, wl=375)
a.show()