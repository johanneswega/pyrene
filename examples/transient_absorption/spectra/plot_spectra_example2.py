from pyrene.transient_absorption import Spectra

s = Spectra(files=['dA1.npy', 'dA2.npy'], x_cuts=[(400, 750)], delay=[1, 1], wn=False,
            scatter=[None, (510, 540)], ma=[True, True],
            norm=[True, True], norm_at=[417, 417], min_norm=[True, True], 
            steady_state=[['abs.txt', (360, 700), -5.0, 'b', 'Abs.']])

# steady state spectra (optional) provided as lists [filename, cuts, scale factor, label]
# absorption/emission spectra can be exported from respectrive classes with export method
s.show()