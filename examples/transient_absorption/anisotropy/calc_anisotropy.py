from pyrene.transient_absorption import Contour, Spectra, Kinetics

### compare spectra ###
s = Spectra(files=['par.npy', 'perp.npy', 'MA.npy'],
          figsize=[20, 5],
          x_cuts=[(400, 780)],     
          overview=True,
          outside=True,
          titles=[r'$\Delta A_{\parallel}$', r'$\Delta A_{\perp}$', r'$\Delta A_{\text{MA}}$'],   
          ylim=[-7.5, 7.5],       
          steady_state=[['abs.txt', (400, 800), -6.5, 'b', 'Abs.']])
s.export()
s.show()

### compare contours ###
c = Contour(files=['par.npy', 'perp.npy', 'MA.npy'],
          y_cuts=[(0.3, 1800)],    
          figsize=[20, 5],
          x_cuts=[(400, 780)],   
          lines=[True], 
          yscale='log',  
          titles=[r'$\Delta A_{\parallel}$', r'$\Delta A_{\perp}$', r'$\Delta A_{\text{MA}}$'],   
          scale=[(-8, 8)])
c.show()

### calculate constructed magic angle signal ###
from pyrene.standard.misc import load, save
t, wl, dApar = load('par.npy')
t, wl, dAperp = load('perp.npy')
dA_iso = (1/3)*(dApar + 2*dAperp)
save('MA_calc.npy', t, wl, dA_iso)

### plot kinetics at wavelength to see if MA set correctly ###
wl_look = 400
k = Kinetics(files=['MA.npy', 'MA_calc.npy'], x_cuts=[(0.3, 1800)], xscale='log', wavelength=[wl_look, wl_look],
             labels=[r'$\Delta A_{\text{MA}}$', r'$(1/3) \cdot (\Delta A_{\parallel} + 2 \Delta A_{\perp})$'], figsize=[5, 3.5])
k.show()

### calculate anisotropy ###
r = (dApar - dAperp)/(dApar + 2*dAperp)
save('r.npy', t, wl, r)