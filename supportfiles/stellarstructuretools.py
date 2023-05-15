import numpy as np
from scipy.optimize import brentq
import supportfiles.config as config
#reading in Prof. Schlaufman's constants file
exec(open("./supportfiles/constants.py").read())

def eng_gen_pp(rho, T):
    """
    pp-chain energy generation rate
    inputs: rho (float): density (g cm^-3)
             T (float): temperature (K)
             
    returns: epsilon_pp (float): energy generation rate for the pp-chain (erg g^-1 s^-1)
                                 or np.nan if the rate is negative
    """
    f11 = np.exp(5.92*(10**-3) *2*4* ((rho* ((T/(10**7))**-3))**0.5))
    T9 = T/(10**9)
    g11 = 1 + (3.82*T9) + (1.51*(T9**2)) + (0.144*(T9**3)) - (0.0114*(T9**4))
    E_pp = (2.57 * (10**4) * f11 * g11 * rho * (config.X**2)
            *(T9**(-2/3)) * np.exp(-3.381*(T9**(-1/3))))
    return (2.57*(10**4)*f11*g11*rho*(config.X**2) 
            *(T9**(-2/3))*np.exp(-3.381*(T9**(-1/3))))

def eng_gen_cno(rho, T):    
    """
    CNO cycle energy generation rate
    inputs: rho (float): density (g cm^-3)
             T (float): temperature (K)
             
    returns: epsilon_CNO (float): energy generation rate for the CNO cycle (erg g^-1 s^-1) 
                                  or np.nan if the rate is negative
    """
    T9 = T/(10**9)
    g14_1 = 1 - (2.00*T9) + (3.41*(T9**2)) - (2.43*(T9**3))
    X_CNO = config.Z #(should be X_C + X_N + X_O instead of Z)
    E_CNO = (8.24*(10**25) * g14_1* X_CNO*config.X*rho*(T9**(-2/3))
            *np.exp(-15.231*(T9**(-1/3)) - ((T9/0.8)**2)))
    return E_CNO

def density(P, T):
    """
    find the density
    inputs: P (float): pressure (dyne cm^-2)
            T (float): temperature (K)
            
    returns: rho (float): density (g cm^-3)
    """
    rho = ((P-(a*(T**4)/3))*config.mu)/(Na*k*T) #from 3.104
    return np.where(rho<=0, 1e-1000, rho)

def pressure_opacity(rho, T, R):
    """
    find pressure from opacity
    inputs: rho (float): density (g cm^-3)
            T (float): temperature (K)
            R (float): total radius (cm)
            M (float): total mass (g)
            
    returns: P_kappa (float): pressure from opacity (dyne cm^-2)
    """
    g = (G*config.MASS)/(R**2)
    kappa = 10**(config.interp(np.log10(rho), np.log10(T)))
    return (2/3)*(g/kappa)

def pressure_total(rho, T):
    """
    find the total pressure
    inputs: rho (float): density (g cm^-3)
            T (float): temperature (K)
            
    returns: P (float): total pressure (dyne cm^-2)
    """
    return ((rho*Na*k*T)/config.mu)+((a*(T**4))/3)

def opacity_vs_total_pressure(rho, T, R):
    """
    find the contribution to total pressure from the ideal gas 
    (as opposed to the radiation/opacity)
    inputs: rho (float): density (g cm^-3)
            T (float): temperature (K)
            
    returns: fractional amount of pressure due to ideal gas (float) (unitless)
    """
    return (1 - (pressure_opacity(rho, T, R)/pressure_total(rho, T)))

def center_vals(m, P, T):
    """
    find the central values from the boundary conditions
    at a mass point m near the center of a star, find L, P, r, & T
    inputs: m (float): mass point (g)
            P (float): guess for the central pressure (dyne cm^-2)
            T (float): guess for the central temperature (K)
    returns: np.array of 
            L (float): luminosity emitted from a sphere enclosing m (ergs s^-1)
            P (float): pressure at the surface of a sphere enclosing m (dyne cm^-2)
            T (float): temperature at the surface of a sphere enclosing m (K)
            R (float): radius of a sphere enclosing m (cm)
    """
    rho_c = density(P, T)
    if rho_c <= 0:
        return np.array([np.nan, np.nan, np.nan, np.nan])
    
    r = (((3*m)/(4*np.pi*rho_c))**(1/3))
    P = ((-3*G)/(8*np.pi))*(((4/3)*np.pi*rho_c)**(4/3))*((m)**(2/3))+P
    energy_gen = eng_gen_pp(rho_c, T)+eng_gen_cno(rho_c, T)
    L = energy_gen * (m)
    
    try:
        print(rho_c, T)
        kappa = 10**(config.interp(np.log10(rho_c), np.log10(T)))
    except:
        print("Failed to calculate kappa. Are you trying to interpolate out of bounds?")
        return np.array([np.nan, np.nan, np.nan, np.nan])

    _, Delta_ad, Delta = Delta_finder(m, L, P, T, r, kappa)
    if Delta > Delta_ad:
        T = np.exp((((np.pi/6)**(1/3))*G*((Delta*(rho_c**(4/3)))/P)
            *(m**(2/3)))+np.log(T)) #CONVECTIVE
    else:
        T = (((-1/(2*a*c))*((3/(4*np.pi))**(2/3))
              *kappa*energy_gen*(rho_c**(4/3))*((m)**(2/3)))
              +T**4)**0.25 #RADIATIVE
    return np.array([L, P, T, r])

def surface_vals(m, L, R):
    """
    find the surface values from boundary conditions
    inputs: m (float): mass point (g)
            L (float): guess for the total luminosity (ergs s^-1)
            R (float): guess for the total radius (cm)
    returns: np.array of 
            L (float): total luminosity (ergs s^-1)
            P (float): pressure at the surface (dyne cm^-2)
            T (float): temperature at the surface (K)
            R (float): total radius (cm)
    """
    T = (L/(np.pi * (R**2) * a * c))**0.25
    #if you want to see info about the convergence, set full_output = True, add another variable to store
    #the ouptut in, and print that variable
    rho = brentq(opacity_vs_total_pressure, 1e-12, 1e-6, args=(T, R))
    P = pressure_opacity(rho, T, R)
    return np.array([L, P, T, R]) #just returning the same radius & luminosity you feed in...

def derivs(m, dep_vs):
    """
    find the derivatives for the four equations of state
    inputs: m (float): mass point (grams)
            dep_vs (list) of:
                L (float): luminosity emitted from a sphere enclosing m (ergs s^-1)
                P (float): pressure at the surface of a sphere enclosing m (dyne cm^-2)
                T (float): temperature at the surface of a sphere enclosing m (K)
                R (float): radius of a sphere enclosing m (cm)
    returns: list of:
                dL/dm (float): derivative of luminosity with respect to mass at mass point m
                dP/dm (float): derivative of pressure with respect to mass at mass point m
                dr/dm (float): derivative of radius with respect to mass at mass point m
                dT/dm (float): derivative of temperature with respect to mass at mass point m
    """
    L, P, T, r = dep_vs
    rho = density(P, T)
    kappa = 10**(config.interp(np.log10(rho), np.log10(T)))
    dLdm = eng_gen_pp(rho, T) + eng_gen_cno(rho, T)
    dPdm = -((G*m)/(4*np.pi*(r**4)))
    drdm = 1/(4*np.pi*(r**2)*rho)
    dTdm, _, _ = Delta_finder(m, L, P, T, r, kappa)

    return np.array([dLdm, dPdm, dTdm, drdm])

def Delta_finder(m, L, P, T, r, kappa, return_Delta_rad = False):
    """
    returns: dTdm
             Delta_ad
             actual Delta
    """
    #print(m, L, P, T, r, kappa)
    Delta_rad = (3/(16*np.pi*a*c)) * ((P*kappa)/(T**4)) * (L/(G*m))
    Delta = np.where(Delta_rad > config.Delta_ad, Delta_rad, config.Delta_ad)
    dTdm = -((G*m*T)/(4*np.pi*(r**4)*P))*Delta
    #if Delta_rad > config.Delta_ad:
    #    Delta = config.Delta_ad
    #    dTdm = -((G*m*T)/(4*np.pi*(r**4)*P))*Delta
    #else:
    #    Delta = Delta_rad
    #    dTdm = -((G*m*T)/(4*np.pi*(r**4)*P))*Delta
    if return_Delta_rad:
        return dTdm, config.Delta_ad, Delta_rad, Delta
    else:
        return dTdm, config.Delta_ad, Delta

"""    if Delta_rad > Delta_ad:
        dTdm = -((G*m*T)/(4*np.pi*(r**4)*P))*Delta_ad
        return dTdm, Delta_ad, Delta_ad
    else:
        return dTdm, Delta_ad, Delta_rad"""

