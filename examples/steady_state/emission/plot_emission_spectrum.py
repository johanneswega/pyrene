from pyrene.steady_state import Emission

e = Emission(
    files=['Data/em_spectrum_1.dat', 'Data/em_spectrum_2.dat'],
    x_cuts=[(350, 600)], norm=[True], TDM=True,
    wn=True)
e.show()