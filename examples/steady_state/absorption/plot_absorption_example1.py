from pyrene.steady_state import Absorption

a = Absorption(
    files=['abs_file1.csv'],
    x_cuts=[(350, 800)],
    wn=True,
    colors=['r'],
    norm=[True],
    labels=['your molecule'])

a.show()