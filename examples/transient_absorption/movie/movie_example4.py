from pyrene.transient_absorption import Movie

m = Movie(files=['data/dA1.npy', 'data/dA2.npy'], 
          y_cuts=[(0.5, 500)], ylim=[-1.2, 1.2], 
          scatter=[(500, 600)], norm=[True], norm_at=[(490, 0.5)],
          ma=[True],
          steady_state=[['data/abs.txt', (360, 700), -1.0, 'b', 'Abs.']])


m.render()