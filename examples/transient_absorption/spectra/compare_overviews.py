from pyrene.transient_absorption import Spectra

# to normalize use "devide" argument manually
s = Spectra(files=['data/dA1.npy', 'data/dA2.npy'], 
            titles=['molecule A', 'molecule B'],
            figsize=(14, 5),
            x_cuts=[(350, 750)], 
            scatter=[(520, 550)],
            devide=[8.5, 8.9],
            export_overview=True,
            ma=[True],
            steady_state_ax=[[['data/abs.txt', (360, 700), -1.0, 'b', 'Abs.'], ['data/abs.txt', (360, 700), +1.0, 'g', 'juhu']],
                              [['data/em.txt', (360, 700), -1.0, 'dodgerblue', 'SE']]],
            overview=True, outside=True)
s.show()