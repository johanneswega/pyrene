from pyrene.steady_state import Absorption

a = Absorption(
    files=['Data/abs_file1.csv'],
    x_cuts=[(350, 800)], 
    eps=[(7.72e3, 524)], # --> provide extinction coefficient (to get it see get_extinction_coefficient.py)
    l=[1],
    wn=False,
    baseline_at=[700],
    colors=['r'],
    labels=['your molecule'])

a.get_concentration(0)
a.show()