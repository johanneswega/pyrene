from pyrene.steady_state import Absorption, Emission

### get reference and sample absorbances at excitation wavelengths ###
a = Absorption(files=['abs/Py2Eth_CHX_N2.csv', 'abs/R6G_EthOH2.csv'], wn=False, figsize=(5,4),
               x_cuts=[(280, 700)], colors=['r', 'b'], 
               baseline_at=[500, 700],
               labels=['sample', 'reference (Rh6G in EtOH)'])
# sample excited at 330 nm
A_S = a.find(0, 330)
# reference excited at 480 nm
A_R = a.find(1, 480)
a.show()

### plot emission spectra and calculate fluoresence quantum yield ###
e = Emission(files=['em/Py2Eth_CHX_under_N2_ex_330.dat', 'em/RG6_EthOH_ex_480.dat'], wn=False, figsize=(5,4), 
             colors=['r', 'b'], corr=True, labels=['sample', 'reference (Rh6G in EtOH)'])

e.calc_FluQY(file_ind_S=0, file_ind_R=1, A_S=A_S, A_R=A_R, phi_R=0.916, solvents=['CHX', 'EtOH'])
e.show()