from pyrene.steady_state import Absorption

a = Absorption(files=['abs_file1.csv', 'abs_file2.csv'],    # file names 
             x_cuts=[(300, 800), (300, 800)],               # wavelength cuts (optional)
             baseline_at=[700, 700],                       # wavelength to correct baseline at (optional)
             colors=['r', 'b'],                             # colors for the files 
             labels=['molecule 1', 'molecule 2']           # labels for the files
             )                             
a.plot_diff(0, 1)
a.show()
