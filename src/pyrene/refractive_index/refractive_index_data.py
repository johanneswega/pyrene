import numpy as np
import scipy.constants as sc
from scipy.optimize import brentq

# sellmeier parameters for matrials
# look which sellmeier eqn to use
# input wavelength in nm into functions
# dn/dl will be returned in µm-1
# d2n/dl2 will be returned in µm-2
# vg will be returned in mm/fs when putting in dn/dl in µm-1
# GVD will be returned in fs^2/mm when putting in d2n/dl in µm-2

# Quartz: Sellmeier1
# ACN: Sellmeier2
# CHF (chloroform): Sellmeier3 
# CHX: Sellmeier4 
# THF: Sellmeier3
# EtOH: Sellmeier3
# OctOH: Sellmeier3
# Hex: Sellmeier5
# BZN: Sellmeier5
# TOL: Sellmeier5
# BK7: Sellmeier1
# SF10: Sellmeier1
# LAK21: Sellmeier1
# CaF2: Sellmeier1
# BBO: Sellmeier6

# types of Sellmeier equations: 
# 1.) n = [A*x**2/(x**2 - B) + C*x**2/(x**2 - D) + E*x**2/(x**2 - F) + 1]**(0.5)
# 2.) n = (A - B*wl**2 + C*wl**(-2) + D*wl**(-4) + E*wl**(-6))**0.5
# 3.) n = ((A*x**2/(x**2 - B)) + (C*x**2/(x**2 - D)) + 1)**(0.5)
# 4.) n = A + B*x**(-2) + C*x**(-4)
# 5.) n = ((A*x**2 / (x**2 - B)) + 1)**(0.5)
# 6.) n = ((A + B/(x**2 - C) - D*x**2))**(0.5)

material = {'quartz': [0.6961663, 0.0684043**2, 0.4079426, 0.1162414**2, 0.8974794, 9.896161**2],
            'ACN': [1.78207648, 0.001147808, 0.003252916, 0.000655275, 0.000001251],
            'CHF': [1.04647, 0.01048, 0.00345, 0.15207], 
            'CHX': [1.41545, 0.00369, 0.00004],
            'THF': [0.9544, 0.0935**2, 0.0147, 3.450**2],
            'EtOH': [0.8312, 0.0964**2, 0.0192**2, 2.97**2],
            'OctOH': [1.013, 0.0965**2, 0.0112, 2.90**2],
            'Hex': [0.8660, 0.0988**2],
            'BZN': [1.180, 0.138**2],
            'TOL': [1.171, 0.134**2],
            'BK7': [1.03961212, 0.00600069867, 0.231792344, 0.0200179144, 1.01046945, 103.560653],
            'SF10': [1.62153902, 0.0122241457, 0.256287842, 0.0595736775, 1.64447552, 147.468793],
            'LAK21': [1.22718116, 0.00602075682, 0.420783743, 0.0196862889, 1.01284843, 88.4370099],
            'CaF2': [0.5675888, 0.050263605**2, 0.4710914, 0.1003909**2, 3.8484723, 34.649040**2],
            'BBO(o)': [2.7405, 0.0184, 0.0179, 0.0155],
            'BBO(e)': [2.3730, 0.0128, 0.0156, 0.0044]}

# Sellmeier type 6
def Sellmeier6(wl, A, B, C, D):
    wl = wl / 1000
    x = wl
    n = np.sqrt(A + B / (x**2 - C) - D * x**2)
    return n

def Sellmeier6_first_derivative(wl, A, B, C, D):
    wl = wl / 1000
    x = wl
    n = np.sqrt(A + B / (x**2 - C) - D * x**2)
    dndwl = ( -B * x / (x**2 - C)**2 - D * x) / n
    return dndwl

def Sellmeier6_second_derivative(wl, A, B, C, D):
    wl = wl / 1000
    x = wl
    n_squared = A + B / (x**2 - C) - D * x**2
    n = np.sqrt(n_squared)
    f_prime = -2 * B * x / (x**2 - C)**2 - 2 * D * x
    f_second = (2 * B * (3 * x**2 + C) / (x**2 - C)**3 - 2 * D)
    dn2dwl2 = (f_second / (2 * n)- f_prime**2 / (4 * n**3))
    return dn2dwl2

# sellmeier for Hex
def Sellmeier5(wl, A, B):
    wl = wl/1000
    x = wl 
    n = ((A*x**2 / (x**2 - B)) + 1)**(0.5)
    return n

def Sellmeier5_first_derivative(wl, A, B):
    wl = wl/1000
    x = wl
    dndwl = -(A * B * x) / ((x**2 - B)**2 * np.sqrt((A * x**2) / (x**2 - B) + 1))
    return dndwl

def Sellmeier5_second_derivative(wl, A, B):
    wl = wl/1000
    x = wl
    dn2dwl2 = (A * B * ((3 * A + 3) * x**4 - 2 * B * x**2 - B**2)) / ((x**2 - B)**4 * ((A * x**2) / (x**2 - B) + 1)**(3 / 2))
    return dn2dwl2

# sellmeier for CHX
def Sellmeier4(wl, A, B, C):
    wl = wl/1000
    x = wl
    n = A + B*x**(-2) + C*x**(-4)
    return n

def Sellmeier4_first_derivative(wl, A, B, C):
    wl = wl/1000
    x = wl
    dndwl = -(2 * B * x**2 + 4 * C) / x**5
    return dndwl

def Sellmeier4_second_derivative(wl, A, B, C):
    wl = wl/1000
    x = wl
    dn2dwl2 = (6 * B * x**2 + 20 * C) / x**6
    return dn2dwl2

# type of Sellmeier of chloroform 
def Sellmeier3(wl, A, B, C, D):
    wl = wl/1000
    x = wl
    n = ((A*x**2/(x**2 - B)) + (C*x**2/(x**2 - D)) + 1)**(0.5)
    return n

def Sellmeier3_first_derivative(wl, A, B, C, D):
    wl = wl/1000
    x = wl
    dndwl = ((2 * C * x) / (x**2 - D) - (2 * C * x**3) / (x**2 - D)**2 + (2 * A * x) / (x**2 - B) - (2 * A * x**3) / (x**2 - B)**2) / (2 * np.sqrt((C * x**2) / (x**2 - D) + (A * x**2) / (x**2 - B) + 1))
    return dndwl

def Sellmeier3_second_derivative(wl, A, B, C, D):
    wl = wl/1000
    x = wl
    dn2dwl2 = ((2 * C) / (x**2 - D) - (10 * C * x**2) / (x**2 - D)**2 + (8 * C * x**4) / (x**2 - D)**3 + (2 * A) / (x**2 - B) - (10 * A * x**2) / (x**2 - B)**2 + (8 * A * x**4) / (x**2 - B)**3) / (2 * np.sqrt((C * x**2) / (x**2 - D) + (A * x**2) / (x**2 - B) + 1)) - ((2 * C * x) / (x**2 - D) - (2 * C * x**3) / (x**2 - D)**2 + (2 * A * x) / (x**2 - B) - (2 * A * x**3) / (x**2 - B)**2)**2 / (4 * ((C * x**2) / (x**2 - D) + (A * x**2) / (x**2 - B) + 1)**(3 / 2))
    return dn2dwl2

# type of Sellmeier used for ACN
def Sellmeier2(wl, A, B, C, D, E):
    wl = wl/1000
    n = (A - B*wl**2 + C*wl**(-2) + D*wl**(-4) + E*wl**(-6))**0.5
    return n

def Sellmeier2_first_derivative(wl, A, B, C, D, E):
    wl = wl / 1000  # Convert to microns
    x = wl
    dn_dwl = -(B * x**8 + C * x**4 + 2 * D * x**2 + 3 * E) / (x**7 * np.sqrt(-B * x**2 + C / x**2 + D / x**4 + E / x**6 + A))
    return dn_dwl

def Sellmeier2_second_derivative(wl, A, B, C, D, E):
    wl = wl / 1000  # Convert to microns
    x = wl
    d2n_dwl2 = -(A * B * x**14 + 6 * B * C * x**12 + (15 * B * D - 3 * A * C) * x**10 + (28 * B * E - 10 * A * D - 2 * C**2) * x**8 + (-21 * A * E - 9 * C * D) * x**6 + (-18 * C * E - 6 * D**2) * x**4 - 19 * D * E * x**2 - 12 * E**2) / (x**14 * (-B * x**2 + C / x**2 + D / x**4 + E / x**6 + A)**(3 / 2))
    return d2n_dwl2

# type of Sellmeier equation for quartz
def Sellmeier1(wl, A, B, C, D, E, F):
    # wavelength is in um
    wl = wl/1000
    x = wl
    n = (A*x**2/(x**2 - B) + C*x**2/(x**2 - D) + E*x**2/(x**2 - F) + 1)**(0.5)
    return n

# in um^-1
def Sellmeier1_first_derivative(wl, A, B, C, D, E, F):
    wl = wl/1000
    x = wl
    dndwl = ((2 * E * x) / (x**2 - F) - (2 * E * x**3) / (x**2 - F)**2 + (2 * C * x) / (x**2 - D) - (2 * C * x**3) / (x**2 - D)**2 + (2 * A * x) / (x**2 - B) - (2 * A * x**3) / (x**2 - B)**2) / (2 * np.sqrt((E * x**2) / (x**2 - F) + (C * x**2) / (x**2 - D) + (A * x**2) / (x**2 - B) + 1))
    return dndwl

# in um^-2
def Sellmeier1_second_derivative(wl, A, B, C, D, E, F):
    wl = wl/1000
    x = wl
    dn2dwl2 = ((2 * E) / (x**2 - F) - (10 * E * x**2) / (x**2 - F)**2 + (8 * E * x**4) / (x**2 - F)**3 + (2 * C) / (x**2 - D) - (10 * C * x**2) / (x**2 - D)**2 + (8 * C * x**4) / (x**2 - D)**3 + (2 * A) / (x**2 - B) - (10 * A * x**2) / (x**2 - B)**2 + (8 * A * x**4) / (x**2 - B)**3) / (2 * np.sqrt((E * x**2) / (x**2 - F) + (C * x**2) / (x**2 - D) + (A * x**2) / (x**2 - B) + 1)) - ((2 * E * x) / (x**2 - F) - (2 * E * x**3) / (x**2 - F)**2 + (2 * C * x) / (x**2 - D) - (2 * C * x**3) / (x**2 - D)**2 + (2 * A * x) / (x**2 - B) - (2 * A * x**3) / (x**2 - B)**2)**2 / (4 * ((E * x**2) / (x**2 - F) + (C * x**2) / (x**2 - D) + (A * x**2) / (x**2 - B) + 1)**(3 / 2))
    return dn2dwl2

def Group_velocity(wl, n, dndwl):
    wl = wl/1000
    # convert c into fs/mm
    A = sc.c*1e-12/n
    B = 1 - (wl/n)*dndwl
    vg = A/B
    # return vg in fs/mm
    return vg

def GVD(wl, dn2dwl2):
    wl = wl/1000
    # speed of light in um/fs
    cfs = sc.c*(10**6 / (10**(15)))
    # GVD in fs^2 / um
    gvd = (wl**3/(2*np.pi*cfs**2))*dn2dwl2
    # in fs^2 / mm
    gvd = gvd*1000
    return gvd

# function to convert wavelength in nm to angular freq in s^-1
def lam_to_om(*args):
    oms = []
    for lami in args:
        oms.append(2*np.pi*sc.c/(lami*1e-9))
    if len(oms)==1:
        return oms[0]
    else:
        return oms

# function to convert angular freq in s^-1 to nm
def om_to_lam(*args):
    lams = []
    for omi in args:
        lams.append(2*np.pi*sc.c/(omi*1e-9))
    if len(lams)==1:
        return lams[0]
    else:
        return lams

# find phase matching angle for SFG/OPA of a uniaxial crystal
# negative crytal typeI ooe
# alpha_deg = angle between pump and signal for non-colinear phase matching
def find_phasematching_angle_typeI_ooe(wl1, n1o, wl2, n2o, wl3, n3o, n3e, alpha_deg=0):
    om1, om2, om3 = lam_to_om(wl1, wl2, wl3)

    # calculate angle between pump and idler
    alpha = np.deg2rad(alpha_deg)
    beta = np.arcsin(((n1o*om1)/(n2o*om2))*np.sin(alpha))
    beta_deg = np.rad2deg(beta)

    # brentq function finds roots of function between two boundraries
    # the function f(x) of which we want to find x = 0 is: 
    f = lambda theta: (np.cos(theta)**2/n3o**2 + np.sin(theta)**2/n3e**2)**(-1/2) * om3 - n1o * om1 * np.cos(alpha) - n2o * om2 * np.cos(beta)
    try:
        root = brentq(f, 0, np.pi/2)
        theta_deg = np.rad2deg(root)

        return theta_deg, beta_deg
    except ValueError:
        print('No phase matching angle found.')
        return np.nan, np.nan

# find phase matching angle for SFG/OPA of a uniaxial crystal
# negative crytal typeII oee
# alpha_deg = angle between pump and signal for non-colinear phase matching
def find_phasematching_angle_typeII_oee(wl1, n1o, wl2, n2o, n2e, wl3, n3o, n3e, alpha_deg=0):
    om1, om2, om3 = lam_to_om(wl1, wl2, wl3)
    alpha = np.deg2rad(alpha_deg)

    # brentq function finds roots of function between two boundraries
    # the function f(x) of which we want to find x = 0 is: 
    f = lambda theta: (np.cos(theta)**2/n3o**2 + np.sin(theta)**2/n3e**2)**(-1/2) * om3 - n1o * om1 * np.cos(alpha) - (np.cos(theta)**2/n2o**2 + np.sin(theta)**2/n2e**2)**(-1/2) * om2 * np.cos(np.arcsin(((n1o*om1)/(((np.cos(theta)**2/n2o**2 + np.sin(theta)**2/n2e**2)**(-1/2))*om2))*np.sin(alpha)))
    try:
        theta = brentq(f, 0, np.pi/2)
        theta_deg = np.rad2deg(theta)

        # calculate angle between pump and idler
        beta = np.arcsin(((n1o*om1)/(((np.cos(theta)**2/n2o**2 + np.sin(theta)**2/n2e**2)**(-1/2))*om2))*np.sin(alpha))
        beta_deg = np.rad2deg(beta)

        return theta_deg, beta_deg
    except ValueError:
        print('No phase matching angle found.')
        return np.nan, np.nan

# find phase matching angle for SFG/OPA of a uniaxial crystal
# negative crytal typeII eoe
# alpha_deg = angle between pump and signal for non-colinear phase matching
def find_phasematching_angle_typeII_eoe(wl1, n1o, n1e, wl2, n2o, wl3, n3o, n3e, alpha_deg=0):
    om1, om2, om3 = lam_to_om(wl1, wl2, wl3)

    # calculate angle between pump and idler
    alpha = np.deg2rad(alpha_deg)

    # brentq function finds roots of function between two boundraries
    # the function f(x) of which we want to find x = 0 is: 
    f = lambda theta: (np.cos(theta)**2/n3o**2 + np.sin(theta)**2/n3e**2)**(-1/2) * om3 - (np.cos(theta)**2/n1o**2 + np.sin(theta)**2/n1e**2)**(-1/2) * om1 * np.cos(alpha) - n2o * om2 * np.cos(np.arcsin(((((np.cos(theta)**2/n1o**2 + np.sin(theta)**2/n1e**2)**(-1/2))*om1)/(n2o*om2))*np.sin(alpha)))
    try:
        theta = brentq(f, 0, np.pi/2)
        theta_deg = np.rad2deg(theta)

        beta = np.arcsin(((((np.cos(theta)**2/n1o**2 + np.sin(theta)**2/n1e**2)**(-1/2))*om1)/(n2o*om2))*np.sin(alpha))
        beta_deg = np.rad2deg(beta)

        return theta_deg, beta_deg
    except ValueError:
        print('No phase matching angle found.')
        return np.nan, np.nan

# calculate, ne(theta), ne'(theta), ne''(theta), vg(theta), gvd(theta)
def angle_dependent_ne(wl_um, theta_deg, no, dno, d2no, ne, dne, d2ne):
    c0 = sc.c
    wl_m = wl_um * 1e-6

    ct2 = np.cos(np.deg2rad(theta_deg))**2
    st2 = 1 - ct2

    Q = ct2/no**2 + st2/ne**2
    dQ = -2 * (ct2*dno/no**3 + st2*dne/ne**3)
    d2Q = (ct2*(6*dno**2/no**4 - 2*d2no/no**3)+ st2*(6*dne**2/ne**4 - 2*d2ne/ne**3))

    ne_theta = Q**-0.5
    dne_theta = -0.5 * Q**-1.5 * dQ
    d2ne_theta = 0.75 * Q**-2.5 * dQ**2 - 0.5 * Q**-1.5 * d2Q

    ng = ne_theta - wl_um * dne_theta
    vg = c0 / ng
    vg_fs_mm = vg*1e-12

    # Convert d²n/dλ² from µm⁻² to m⁻²
    beta2 = wl_m**3 / (2*np.pi*c0**2) * d2ne_theta * 1e12
    gvd_fs2_mm = beta2 * 1e27

    return ne_theta, dne_theta, d2ne_theta, vg_fs_mm, gvd_fs2_mm
