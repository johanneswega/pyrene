import numpy as np 
import scipy.constants as sc 

def E_to_P(E, RP):
    # energy per pulse in J 
    # rep-rate in Hz
    P = E * RP
    print("P = %.3g W = %.3g mW = %.3g µW"%(P, P*1e3, P*1e6))
    
def I_to_E(I):
    # irradiance in W/cm^2 
    # convert to W/m^2 
    I = I*1e4
    # calculate E in V/m 
    E = np.sqrt(2*I/(sc.epsilon_0*sc.c))
    print("I = %.3g W/cm^2 --> E = %.3g V/m"%(I*1e-4, E))

def I_from_laser(En, t, d):
    # energy per pulse in J
    # focus diameter (assume circle) in µm
    # pulse width in s

    # calculate d in cm 
    d = d*1e-4
    # calculate area in cm^2
    A = np.pi * (d/2)**2
    # calculate peak irradiance in W/cm^2 
    I = En/(t*A)
    I_to_E(I)

#I_to_E(1e6)
I_from_laser(1e-6, 100e-15, 100)
E_to_P(1e-6, 1e3)