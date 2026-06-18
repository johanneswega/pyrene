from pyrene.transient_absorption import Movie

m = Movie(files=['data/dA1.npy'],
          scatter=[(520, 560)],     # scatter region to exclude (optional)
          experiment='femto',       # experiment type (optional)
          x_cuts=[(330, 780)],      # wavelength cuts in nm (optinal)
          y_cuts=[(0.3, 1800)],     # time cuts in ps if experiment = 'femto' else in ns (optional)
          colors=['crimson'],   
          labels=['your molecule'],
          movname='movie.mp4',      # file name to save (optional)
          ylim=[-12.5, 25.5],       # imits of y-axis (optional)
          before=True,              # decide to plot traces before (optional)
          steady_state=[['data/abs.txt', (330, 700), -9.5, 'b', 'Abs.']])

# steady state spectra (optional) provided as lists [filename, cuts, scale factor, label]
# absorption/emission spectra can be exported from respectrive classes with class.export() 

m.render()