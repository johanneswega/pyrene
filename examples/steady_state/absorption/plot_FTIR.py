from pyrene.steady_state import Absorption

a = Absorption(
    files=['Data/FTIR_2.0'],
    x_cuts=[(1125, 1275)],
    baseline_at=[(1275)],
    IR=True,
    colors=['r'],
    labels=['your molecule'])

a.show()