from pyrene.transient_absorption import Movie

m = Movie(files=['data/nano.npy'], 
          experiment='nano',
          slicing=[5],          # don't take all delays but only every 5th delay
          normall=True,          
          y_cuts=[(0.5, 500e3)],  # time cuts
          ylim=[-1.5, 1.5], 
          x_cuts=[(370, 720)])

m.render()