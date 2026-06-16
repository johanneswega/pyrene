from pyrene.steady_state import Emission, Absorption

# make emission class
em = Emission(files=['Data/em_spectrum_1.dat'], x_cuts=[(350, 700)], norm=[True], labels=['emission'])

# make class for excitation spectra as slave
# same as emission but no correction required
ex = Emission(files=['Data/ex_spectrum.dat'], slave=True, corr=False, marker=['--'], colors=['k'],
              x_cuts=[(350, 500)], norm=[True], labels=['excitation'])

# make absorption class as slave to plot inside ax of em class
a = Absorption(files=['Data/abs.csv'], x_cuts=[(350, 500)], norm=[True], labels=['absorption'], colors=['b'], slave=True)

# plot absorption and excitation spectra inside em figure
em.plot_absorption(a, ExClass=ex)
em.show()