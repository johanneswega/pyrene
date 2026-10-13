from pyrene.standard.packages import *
from pyrene.refractive_index import *

# example to find phase matching angle theta in a negative 
# uniaxial crystal = here BBO for the case of SFG
# om3 = om1 + om2
# with Type II (oee) phase matching

wl1 = 400
wl2 = 800

om2, om1 = lam_to_om(wl2, wl1)
om3 = om1 + om2
wl3 = om_to_lam(om3)

# look into src/pyrene/refractive_index/refractive_index_data.py which 
# Sellmeier type equation is required for which material

# calculate refractive indices
n1o = Sellmeier6(wl1, *material['BBO(o)'])
n1e = Sellmeier6(wl1, *material['BBO(e)'])

n2o = Sellmeier6(wl2, *material['BBO(o)'])

n3o = Sellmeier6(wl3, *material['BBO(o)'])
n3e = Sellmeier6(wl3, *material['BBO(e)'])

# find phase matching angle theta

# in case of non-colinear beams:
# alpha is the angle between om3 and om1 (i.e. signal-pump angle for OPA)
# beta is the angle between om2 and om1 (i.e. signal-idler angler for OPA)

theta_deg, beta_deg = find_phasematching_angle_typeII_eoe(wl1, n1o, n1e, wl2, n2o, wl3, n3o, n3e, alpha_deg=0.0)
print(f'theta = {theta_deg : .2f}°')