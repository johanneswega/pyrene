from pyrene.transient_absorption import Movie

m = Movie(files=['data/dA1.npy'], 
          before=True,          # whether to plot slices before (optional)
          y_cuts=[(200, 500)],  # time cuts
          ylim=[-5.2, 6.2], 
          scatter=[(500, 600)])

m.render()