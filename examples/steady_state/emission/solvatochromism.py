from pyrene.steady_state import Absorption, Emission
from pyrene.standard.misc import rainbow, sort_with_respect
from pyrene.standard.solvents import solvent_dic
import os 

### load absorption files ###
files_abs = [i for i in os.listdir('Data/solvchrom_abs') if not '.DS' in i]
# get labels/solvents 
lables = [s[:s.find('.csv')] for s in files_abs]
# sort with respect to dielectric constnat
eps_r = [solvent_dic[s][1] for s in lables]
eps_r, files_abs, lables = sort_with_respect(eps_r, files_abs, lables)
files_abs = ['Data/solvchrom_abs/' + s for s in files_abs]

### load emission files ###
files_em = [i for i in os.listdir('Data/solvchrom_em') if not '.DS' in i]
# get labels/solvents 
lables = [s[:s.find('.dat')] for s in files_em]
# sort with respect to dielectric constnat
eps_r = [solvent_dic[s][1] for s in lables]
eps_r, files_em, lables = sort_with_respect(eps_r, files_em, lables)
files_em = ['Data/solvchrom_em/' + s for s in files_em]

# get colors
colors = rainbow(files_em)

### absorption figure = slave ###
a = Absorption(files=files_abs, slave=True, waterfall=2,
    marker=['--' for _ in files_abs], yticks=False, x_cuts=[(550, 800) for _ in files_abs],
    baseline_at=[800 for _ in files_abs], outside=True, colors=colors, norm=[True for _ in files_abs], labels=lables)

### emission figure = master ###
e = Emission(files=files_em, figsize=(6, 9), waterfall=2, yticks=False, x_cuts=[(550, 1200) for _ in files_em], corr=False,
              outside=True, colors=colors, norm=[True for _ in files_abs], labels=lables)

e.plot_absorption(AbsClass=a)
e.solv_chrom(AbsClass=a)
e.show()