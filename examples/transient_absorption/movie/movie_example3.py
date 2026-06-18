from pyrene.transient_absorption import Movie

m = Movie(files=['data/dA1.npy', 'data/dA2.npy'],
          scatter=[(520, 560), (520, 560)],
          colors=['darkorange', 'crimson'],
          labels=[r'Mol. 1', r'Mol. 2'],
          movname='comp.mp4',
          x_cuts=[(330, 750)],
          y_cuts=[(0.3, 1800)],
          ylim=[-0.65, 1.2],
          before=False,
          norm=[True], # nomalize spectra (optional)
          norm_at=[(380, 0.3)], # normalize spectra at a specific (wavelength, delay)
          yticks=False) # remove y-ticks of the plot

m.render()