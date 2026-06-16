from pyrene.steady_state import Absorption

a = Absorption(
    files=['Data/abs_file1.csv'],
    x_cuts=[(350, 800)],
    wn=True,
    colors=['r'],
    labels=['your molecule'])

a.show()