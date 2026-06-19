from pyrene.transient_absorption import Contour

c = Contour(files=['data/dA1.npy'], 
            scatter=[(520, 550)], 
            y_cuts=[(-0.5, 1000)], x_cuts=[(300, 500)], 
            scale=[(-30, 30)],
            lines=[True],
            yscale='symlog')
#c.export()
c.show()