from pyrene.transient_absorption import Spectra

# to normalize use "devide" argument manually
s = Spectra(files=['dA1.npy'], x_cuts=[(350, 750)], overview=True, scatter=[(520, 540)], outside=True,
            steady_state=[['abs.txt', (360, 700), -10.0, 'b', 'Abs.']])
s.show()