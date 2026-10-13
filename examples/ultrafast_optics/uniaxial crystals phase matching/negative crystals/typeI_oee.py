from pyrene.standard.packages import *
from pyrene.refractive_index import *

# example to find phase matching angle theta in a negative 
# uniaxial crystal = here BBO for the case of SHG
# om3 = om1 + om2
# with Type I (ooe) phase matching

wl3 = 400
wl1 = 800

om3, om1 = lam_to_om(wl3, wl1)
om2 = om3 - om1
wl2 = om_to_lam(om2)

# look into src/pyrene/refractive_index/refractive_index_data.py which 
# Sellmeier type equation is required for which material

# calculate refractive indices
n3o = Sellmeier6(wl3, *material['BBO(o)'])
n3e = Sellmeier6(wl3, *material['BBO(e)'])

n1o = Sellmeier6(wl1, *material['BBO(o)'])
n2o = Sellmeier6(wl2, *material['BBO(o)'])

# find phase matching angle theta

# in case of non-colinear beams:
# alpha is the angle between om3 and om1 (i.e. signal-pump angle for OPA)
# beta is the angle between om2 and om1 (i.e. signal-idler angler for OPA)

theta_deg, beta_deg = find_phasematching_angle_typeI_ooe(wl1, n1o, wl2, n2o, wl3, n3o, n3e, alpha_deg=0)

print(f'theta = {theta_deg : .2f}°')