from pyrene.steady_state import Emission, Absorption

# make emission class
e = Emission(files=['Data/em_spectrum_1.dat'], x_cuts=[(350, 700)], norm=[True], labels=['emission'], yticks=False)

# make absorption class as slave to plot inside ax of em class
a = Absorption(files=['Data/abs.csv'], x_cuts=[(300, 500)], norm=[True], labels=['absorption'], fill=[True], marker=[None], slave=True)

e.plot_absorption(a)
e.show()