from pyrene.transient_absorption import Spectra

s = Spectra(files=['data/dA1.npy', 'data/dA2.npy'], x_cuts=[(400, 750)], delay=[1, 1], wn=False,
            scatter=[None, (510, 540)], 
            norm=[True], norm_at=[417, 417], 
            steady_state=[['data/abs.txt', (360, 700), -5.0, 'b', 'Abs.']])

# steady state spectra (optional) provided as lists [filename, cuts, scale factor, label]
# absorption/emission spectra can be exported from respectrive classes with export method
s.show()