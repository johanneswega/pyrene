from pyrene.transient_absorption import Spectra

s = Spectra(files=['data/dA1.npy', 'data/dA2.npy'], delay=[1, 1], norm=[False],
            x_cuts=[(300, 500)])
s.export()
s.show()