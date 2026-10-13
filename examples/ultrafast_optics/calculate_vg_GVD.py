from pyrene.refractive_index import *

# look into src/pyrene/refractive_index/refractive_index_data.py which 
# Sellmeier type equation is required for which material

# wavelength in nm as input for the functions
wl = 400
n = Sellmeier6(wl, *material['BBO(o)'])

# first and second derivatives of the refractive index at specified wavelength
# in units of µm-1 and µm-2
dndl = Sellmeier6_first_derivative(wl, *material['BBO(o)'])
dn2dl = Sellmeier6_second_derivative(wl, *material['BBO(o)'])

# group velocity in mm/fs
vg = Group_velocity(wl, n, dndl)
# group index c/vg
ng = (sc.c*1e-12)/vg

# group velocity dispersion in fs^2/mm
gvd = GVD(wl, dn2dl)

print("")
print(f"l = {wl:.5g} nm")
print(f"n = {n:.5g} ")
print(f"dn/dl = {dndl:.5g} µm-1")
print(f"d2n/dl2 = {dn2dl:.5g} µm-2")
print(f"vg = {vg*1000:.5g} mm/ps")
print(f"ng = {ng:.5g}")
print(f"GVD = {gvd:.5g} fs^2/mm")
print("")