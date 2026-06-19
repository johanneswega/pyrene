from pyrene.transient_absorption import Contour

c = Contour(files=['data/dA1.npy', 'data/dA2.npy'], 
            figsize=[12, 5],
            scatter=[(520, 550)], 
            y_cuts=[(-0.5, 1000)], x_cuts=[(300, 500)], 
            scale=[(-1.2, 1.2)],
            norm=[True],
            norm_at=[(350, 0.5)], # normat takes a tuple here 
            extend=['both'],
            savefig='test.svg',
            yscale='symlog')
c.export()
c.show()