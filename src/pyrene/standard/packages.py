import numpy as np
np.seterr(divide='ignore')
import os
from matplotlib import pyplot as plt
plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'style.mplstyle'))
from matplotlib.animation import FuncAnimation, writers
from scipy.optimize import curve_fit, least_squares
import scipy.constants as sc