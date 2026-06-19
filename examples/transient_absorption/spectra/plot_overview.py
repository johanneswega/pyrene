from pyrene.transient_absorption import Spectra

# to normalize use "devide" argument manually
s = Spectra(files=['data/dA_nano.npy'], x_cuts=[(370, 750)], experiment='nano',
            overview=True, outside=True, export_overview=True,
            steady_state=[['data/abs.txt', (360, 700), -10.0, 'b', 'Abs.']])
s.show()