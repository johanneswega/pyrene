from pyrene.transient_absorption import Kinetics

k = Kinetics(files=['data/dA1.npy', 'data/dA2.npy'], wavelength=[350, 350],
             outside=True, figsize=[6, 3.5], savefig='compare_kinetics.png',
             x_cuts=[(0.2, 1800)], yscale='log')
k.show()
k.export()