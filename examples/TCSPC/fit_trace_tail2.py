from pyrene.TCSPC import TCSPC

t = TCSPC(files=['data/MNP.asc'], yscale='log', ylim=[1, 2e3], x_cuts=[(-5, 250)])

t.fit(file_index=0, tail=1, p0=[100, 2, 1000, 50])
t.show()