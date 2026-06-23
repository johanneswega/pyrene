from pyrene.transient_absorption import Contour

c = Contour(files=['data/data.pdat'], 
            IR = True,
            y_cuts=[(0.5, 1000)],
            lines=[True],
            cmap=['seismic'],
            nlevels=[41],
            scale=[(-10, 10)],
            yscale='log')
c.show()