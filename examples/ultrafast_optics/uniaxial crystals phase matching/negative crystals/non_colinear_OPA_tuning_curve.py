from pyrene.standard.packages import *
from pyrene.refractive_index import *

# in this example we pump an OPA consisting of BBO 
# with the SH of an Yb-Laser at wl3 = 515 nm 
# the goal is to find the signal/pump angle for which 
# phase matching is achieved over the broadest wavelength range
# we consider type I (ooe) phasematching

# specify angle between signal and pump
alpha_deg = 2.3

wl3 = 515
wl1 = np.linspace(600, 1200, 100)

theta = np.zeros(len(wl1))
beta = np.zeros(len(wl1))
idler = np.zeros(len(wl1))
for i in range(len(wl1)):
    om3, om1 = lam_to_om(wl3, wl1[i])
    om2 = om3 - om1
    wl2 = om_to_lam(om2)

    # look into src/pyrene/refractive_index/refractive_index_data.py which 
    # Sellmeier type equation is required for which material

    # calculate refractive indices
    n3o = Sellmeier6(wl3, *material['BBO(o)'])
    n3e = Sellmeier6(wl3, *material['BBO(e)'])

    n1o = Sellmeier6(wl1[i], *material['BBO(o)'])
    n2o = Sellmeier6(wl2, *material['BBO(o)'])

    # find phase matching angle theta

    # in case of non-colinear beams:
    # alpha is the angle between om3 and om1 (i.e. signal-pump angle for OPA)
    # beta is the angle between om2 and om1 (i.e. signal-idler angler for OPA)
    theta_deg, beta_deg = find_phasematching_angle_typeI_ooe(wl1[i], n1o, wl2, n2o, wl3, n3o, n3e, alpha_deg)
    theta[i] = theta_deg
    beta[i] = beta_deg
    idler[i] = wl2

fig, ax = plt.subplots(2,1, figsize=(8,8), gridspec_kw={'height_ratios': [2, 1]})
ax[0].set_title(r'BBO (Type I), $\lambda_\text{pump} = %i \,\text{nm}$'%(int(round(wl3))))
ax[0].plot(wl1, theta, '-g', label=r'$\alpha = %.3g^\circ$'%(alpha_deg))
ax[0].set_ylabel(r'phase matching angle $\theta / ^\circ$')
ax[0].set_xlabel(r'signal wavelength / nm')

ax[1].plot(idler, beta, '-r')
ax[1].set_ylabel(r'pump-idler angle $\beta / ^\circ$')
ax[1].set_xlabel(r'idler wavelength / nm')

ax[1].grid(alpha=0.2)
ax[0].grid(alpha=0.2)
ax[0].legend()
fig.tight_layout()
plt.show()