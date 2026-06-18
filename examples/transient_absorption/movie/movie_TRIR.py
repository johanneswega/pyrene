from pyrene.transient_absorption import Movie

m = Movie(files=['data/data.pdat'],
        y_cuts=[(0.3, 500)],
        IR=True, # for TRIR
        before=True,  
        figsize=(5, 3.5), # size of the figure (optional)
        labels=[r'TCNQ$^{\bullet –}$'],
        colors=['r'], 
        steady_state=[['data/FTIR.txt', (1e3, 3e3), -8, 'b', 'FTIR']],
        movname='TCNQ.mp4', 
        ylim=[-10, 4]) # limit of the y-axis plot

m.render() 