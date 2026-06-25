from pyrene.TCSPC import TCSPC

t = TCSPC(files=['data/dilute.asc', 'data/IRF.asc'], 
          colors=['r', 'gray'], fill=[False, True], marker=['-', None], yscale='log', ylim=[2, 5e4])

# tail = time at which tail fit should start
t.fit(file_index=0, p0=[1e4, 2], tail=2.0)
t.show()