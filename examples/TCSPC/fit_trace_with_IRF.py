from pyrene.TCSPC import TCSPC

t = TCSPC(files=['data/dilute.asc', 'data/IRF.asc'], 
          colors=['r', 'gray'], fill=[False, True], marker=['-', None], yscale='log', ylim=[2, 5e4])

# gaussian fit of the IRF
t.fit_IRF()

# by convention the IRF should always be the last file in files
# just give the number of parameters to fit : [A1, tau1, A2, tau2, ..., An, taun, shift]
# shift = color shift of the IRF
t.fit(file_index=0, p0=[100, 2, 1], IRF=True, IRF_cleanup=4)

t.show()