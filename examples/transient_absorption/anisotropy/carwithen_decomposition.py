from pyrene.standard.misc import load, save

# decompostion of TA spectrum into features with TDMs parallel to pumped TDM = dA_par_mol 
# and TDMs perpendicular to pumped TDM = dA_perp_mol 

# according to method of Carwithen et al. JACS (2026) : https://pubs.acs.org/doi/10.1021/jacs.5c22517

# anisotropy lifetime 
tau_rot = 350

if not tau_rot:
    # in the stationary case where the molecule is considered frozen:
    t, wl, dApar = load('par.npy')
    t, wl, dAperp = load('perp.npy')
    dApar_mol = 2*dApar - dAperp
    dAperp_mol = 3*dAperp - dApar
    save('par_mol.npy', t, wl, dApar_mol)
    save('perp_mol.npy', t, wl, dAperp_mol)
else:
    # more general case where molecules are considered to be able to rotate with tau_rot
    import numpy as np 
    t, wl, dApar = load('par.npy')
    t, wl, dAperp = load('perp.npy')

    dApar_mol = np.zeros((len(t), len(wl)))
    dAperp_mol = np.zeros((len(t), len(wl)))

    for i in range(len(t)):
        F = np.array([[4/15, -2/15], [-2/15, 1/15]]) * np.exp(-t[i]/tau_rot) + np.array([[1/3, 1/3], [1/3, 1/3]])
        F_inv = np.linalg.inv(F)
        dApar_mol[i, :] = F_inv[0,0] * dApar[i, :] + F_inv[0, 1] * dAperp[i, :]
        dAperp_mol[i, :] = F_inv[1, 0] * dApar[i, :] + F_inv[1, 1] * dAperp[i, :]
        
    save('par_mol.npy', t, wl, dApar_mol)
    save('perp_mol.npy', t, wl, dAperp_mol)

# make movie
from pyrene.transient_absorption import Movie
m = Movie(files=['par_mol.npy', 'perp_mol.npy'],
          x_cuts=[(400, 780)],     
          y_cuts=[(0.3, 1800)],    
          outside=[True],
          labels=[r'$\Delta A_{\parallel}^{\text{mol.}}$', r'$\Delta A_{\perp}^{\text{mol.}}$'],   
          colors=['r', 'b'],
          movname='anisotropy_decomp_carwithen.mp4',    
          ylim=[-12.5, 12.5],       
          steady_state=[['abs.txt', (400, 800), -8.5, 'b', 'Abs.'], ['exported/MA/spectra/1_ps.txt', (400, 780), +8.5, 'gray', r'$\Delta A_\text{MA}$ at 1 ps']])
m.render()
