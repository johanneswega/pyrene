from pyrene.transient_absorption import Spectra

s = Spectra(files=['dA1.npy', 'dA2.npy'], delay=[1, 1])
s.show()