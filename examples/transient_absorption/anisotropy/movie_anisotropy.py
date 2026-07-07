from pyrene.transient_absorption import Movie

# calculate magic angle from perpendicular and parallel 
from pyrene.standard.misc import load, save
t, wl, dApar = load('par.npy')
t, wl, dAperp = load('perp.npy')
dA_iso = (1/3)*(dApar + 2*dAperp)
save('MA_calc.npy', t, wl, dA_iso)

m = Movie(files=['par.npy', 'perp.npy', 'MA.npy', 'MA_calc.npy'],
          x_cuts=[(400, 780)],     
          y_cuts=[(0.3, 1800)],    
          outside=[True],
          labels=[r'$\Delta A_{\parallel}$', r'$\Delta A_{\perp}$', r'$\Delta A_{\text{MA}}$', 
                  r'$\frac{1}{3} \cdot (\Delta A_{\parallel} + 2 \Delta A_{\perp})$'],   
          colors=['r', 'b', 'g', 'orange'],
          movname='anisotropy.mp4',    
          ylim=[-7.5, 7.5],       
          steady_state=[['abs.txt', (400, 800), -5.5, 'b', 'Abs.']])
m.render()