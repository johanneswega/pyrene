from pyrene.steady_state import Absorption

a = Absorption(files=['Data/abs_file1.csv', 'Data/abs_file2.csv'],      # file names 
             x_cuts=[(350, 800), (250, 800)],                           # wavelength cuts (optional)
             ma=[False, True],                                          # whether to apply moving average filter (optional)
             ma_npoints = [5, 5],                                       # window size for moving avergae (optional)
             baseline_at=[700, None],                                   # wavelength to correct baseline at (optional)
             colors=['r', 'b'],                                         # colors for the files 
             labels=['molecule 1', 'molecule 2'],                       # labels for the files
             norm=[True, True])                                         # whether to normalize (optional)
a.show()
