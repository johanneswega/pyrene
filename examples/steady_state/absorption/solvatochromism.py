from pyrene.steady_state import Absorption
from pyrene.standard.misc import rainbow, sort_with_respect
from pyrene.standard.solvents import solvent_dic
import os 

# load files
files = [i for i in os.listdir('Data/solvchrom') if not '.DS' in i]
# get labels 
lables = [s[:s.find('.csv')] for s in files]
# sort with respect to dielectric constnat
eps_r = [solvent_dic[s][1] for s in lables]
eps_r, files, lables = sort_with_respect(eps_r, files, lables)
files = ['Data/solvchrom/' + s for s in files]

# get colors
colors = rainbow(files)

# standard figure
a = Absorption(
    files=files,
    x_cuts=[(550, 800) for _ in files],
    baseline_at=[800 for _ in files],
    outside=True,
    colors=colors,
    norm=[True for _ in files],
    labels=lables)
a.show()

# waterfall figure
a = Absorption(
    files=files,
    figsize=(4, 8),
    waterfall=2,
    fill=[True for _ in files],
    yticks=False,
    x_cuts=[(550, 800) for _ in files],
    baseline_at=[800 for _ in files],
    outside=True,
    colors=colors,
    norm=[True for _ in files],
    labels=lables)
a.solvchrom()
a.show()