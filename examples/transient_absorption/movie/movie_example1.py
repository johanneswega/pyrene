from pyrene.transient_absorption import Movie

m = Movie(files=['data/dA1.npy'], y_cuts=[(200, 500)], ylim=[-5.2, 5.2], scatter=[(500, 600)],
          ma=[True], steady_state=[['data/abs.txt', (360, 700), -5.0, 'b', 'Abs.']])

""" m.animate(0)
from matplotlib import pyplot as plt
plt.show() """

m.render()