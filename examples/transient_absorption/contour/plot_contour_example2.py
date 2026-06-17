from pyrene.transient_absorption import Contour

c = Contour(files=['data/dA1.npy'], 
            scatter=[(520, 550)], 
            y_cuts=[(-0.5, 1000)], x_cuts=[(300, 700)], 
            scale=[(-30, 30)],
            nlevels=[51],
            lines=[True],
            white=[2],
            yscale='symlog')
c.show()