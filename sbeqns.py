from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as units
import astropy.constants as const

def k_from_m1m2(m1,m2, P, i, e=0):
    """
    Parameters
    ----------
    i: float
        inclination with units
        
    FIXME: Add in non-zero eccentricity
    """
    f1 = m2**3*np.sin(np.radians(i))**3/(m1+m2)**2
    f2 = m1**3*np.sin(np.radians(i))**3/(m1+m2)**2
    K1 = ((f1*2*np.pi*const.G/P)**(1/3.)).si
    K2 = ((f2*2*np.pi*const.G/P)**(1/3.)).si
    return K1, K2
    
def m1m2_from_k(K1,K2,P,i, e=0):
    M1 = ((P/2/np.pi/const.G/np.sin(i)**3*(K1+K2)**2*K2)/const.M_sun).si.value
    M2 = ((P/2/np.pi/const.G/np.sin(i)**3*(K1+K2)**2*K1)/const.M_sun).si.value
    return M1,M2 
    
if __name__=="__main__":
    P = 2.608223*units.d
    inc = 75*units.deg
    K1mn,K2mn = k_from_m1m2(2.216*const.M_sun, 0.715*const.M_sun,P, inc)
    print("K1: {:5.1f}km/s".format((K1mn/(units.km/units.s)).si.value))
    print("K2: {:5.1f}km/s".format((K2mn/(units.km/units.s)).si.value))
    M1s = []
    M2s = []
    ntry = 1000
    for i in range(ntry):
        K1 = np.random.normal(K1mn.value, 2.7e3)*units.m/units.s
        K2 = np.random.normal(K2mn.value, 4.7e3)*units.m/units.s
        M1,M2 = m1m2_from_k(K1,K2,P,inc)
        M1s.append(M1)
        M2s.append(M2)
    M1s = np.array(M1s)
    M2s = np.array(M2s)
    print("M1: {:5.2f} +/- {:5.2f}".format(np.mean(M1s), np.std(M1s)))
    print("M2: {:5.2f} +/- {:5.2f}".format(np.mean(M2s), np.std(M2s)))
    print("M2/M1: {:5.3f} +/- {:5.3f}".format(np.mean(M2s/M1s), np.std(M2s/M1s)))