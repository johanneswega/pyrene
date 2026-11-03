from pyrene.refractive_index import *

# script to calculate pulse front tilt angle 
# after passage through a prism (assuming symmetric passage; minimum deviation angle)

# let's use SF10 glass 
# and a center wavelength of 1300 nm
wl = 1300
# calculate refractive index
n = Sellmeier1(wl, *material["SF10"])
print(f"n = {n :.5g}")
# calculate group velocity
dndwl = Sellmeier1_first_derivative(wl, *material["SF10"])
vg = Group_velocity(wl, n, dndwl)
print(f"vg = {vg*1e3 :.3g} mm/ps")
# calculate phase velocity 
vph = sc.c*1e-12 / n
print(f"vph = {vph*1e3 :.3g} mm/ps")
# prism apex angle
gamma_deg = 60.6
# calculate pulse-front tilt
Phi = np.arctan(2*sc.c*1e-12 * (1/vg - 1/vph) * np.sin(np.deg2rad(gamma_deg/2)) / (1 - n**2 * (np.sin(np.deg2rad(gamma_deg/2)))**2)**(0.5))
Phi_deg = np.rad2deg(Phi)
print(f"Pulse front tilt Phi = {Phi_deg :.3g}°")