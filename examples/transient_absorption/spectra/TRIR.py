from pyrene.transient_absorption import Spectra

s = Spectra(files=['data/data.pdat'], overview=True, IR=True)
s.show()