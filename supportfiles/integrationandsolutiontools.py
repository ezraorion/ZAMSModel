import numpy as np
from scipy.integrate import solve_ivp
from supportfiles.stellarstructuretools import *
import supportfiles.config as config
#reading in Prof. Schlaufman's constants file
from supportfiles.constants import *

def shootf(guess, mc = 1e-12, mf=0.2, ms = 0.9999, steps=1e6):
    """
    Implementation of the shooting method. Integrates from the center to a 
    fitting point and from the exterior to the same fitting point. Returns 
    a vector of modified residuals describing the disagreement at the 
    fitting point.

    inputs: guess (np.array): with
                L (float): luminosity emitted from a sphere enclosing m 
                           (ergs s^-1)
                P (float): pressure at the surface of a sphere enclosing m 
                           (dyne cm^-2)
                T (float): temperature at the surface of a sphere enclosing m 
                           (K)
                R (float): radius of a sphere enclosing m (cm)
            mc (float): central mass point (fractional)
            mf (float): fitting mass point (fractional)
            ms (float): surface mass point (fractional)
            steps (float): steps for integration. split equally between the 
                           interior and exterior solution

    returns: score (float): agreement between the interior and exterior 
                            solution. should be 0 for a perfect match, with 
                            increasing positive values representing 
                            increasingly poor matches.
    """
    L, P, T, R = guess

    center_guess = center_vals(mc * config.MASS, P, T)
    masses_cen = np.linspace(mc * config.MASS, mf * config.MASS, 
                             num=int(steps/2))
    solc = solve_ivp(derivs, (masses_cen[0], masses_cen[-1]), center_guess, 
                     t_eval = masses_cen)

    surface_guess = surface_vals(config.MASS * ms, L, R)
    masses_surf = np.linspace(config.MASS * ms, config.MASS * mf, 
                  num=int(steps/2))
    sols = solve_ivp(fun = derivs, t_span = (masses_surf[0], masses_surf[-1]),
                     y0 = surface_guess, t_eval = masses_surf)

    modresidual = (solc.y[:, -1] - sols.y[:, -1])/guess
    return modresidual

def solution(bestbc, mc = 1e-12, mf=0.2, ms = 0.9999, steps=1e6, 
             returnedsteps = 1e4):
    """
    Using the best boundary conditions from a run of shootf, integrate and 
    return parameters of interest.

    inputs: bestbc (np.array): with
                L (float): luminosity emitted from a sphere enclosing m 
                           (ergs s^-1)
                P (float): pressure at the surface of a sphere enclosing m 
                           (dyne cm^-2)
                T (float): temperature at the surface of a sphere enclosing m 
                           (K)
                R (float): radius of a sphere enclosing m (cm)
            mc (float): central mass point (fractional)
            mf (float): fitting mass point (fractional)
            ms (float): surface mass point (fractional)
            steps (float): steps for integration. split equally between the 
                           interior and exterior solution
            returnedsteps (float): size of the arrays returned by this 
                                   function (-1)

    returns: msol (array): log mass solution array (log g)
             Lsol (array): log luminosity solution array (log ergs s^-1)
             Psol (array): log pressure solution array (log dyne cm^-2)
             Tsol (array): log temperature solution array (log K)
             Rsol (array): log radius solution array (log cm)
             Esol (array): energy generation rate (erg s^-1)
             kappasol (array): opacity
             Delta_ad (array): adiabatic temperature gradient 
             Delt_rad (array): radiative temperature gradient 
             Delta (array): actual temperature gradient 
             nature (array): is this region convective or radiative?
    """
    L, P, T, R = bestbc

    center_guess = center_vals(mc*config.MASS, P, T)
    masses_cen = np.linspace(mc*config.MASS, mf*config.MASS, num=int(steps/2))
    solc = solve_ivp(derivs, (masses_cen[0], masses_cen[-1]), center_guess, 
                     t_eval = masses_cen)
    
    surface_guess = surface_vals(config.MASS*ms, L, R)
    masses_surf = np.linspace(config.MASS*ms, config.MASS*mf, num=int(steps/2))
    sols = solve_ivp(fun = derivs, t_span = (masses_surf[0], masses_surf[-1]), 
                     y0 = surface_guess, t_eval = masses_surf)
    
    msol0 = np.concatenate((solc.t, np.flip(sols.t)))
    Lsol0 = np.concatenate((solc.y[0], np.flip(sols.y[0])))
    Psol0 = np.concatenate((solc.y[1], np.flip(sols.y[1])))
    Tsol0 = np.concatenate((solc.y[2], np.flip(sols.y[2])))
    Rsol0 = np.concatenate((solc.y[3], np.flip(sols.y[3])))

    msol = np.append(msol0[0::int(steps/returnedsteps)], msol0[-1])
    Lsol = np.append(Lsol0[0::int(steps/returnedsteps)], Lsol0[-1])
    Psol = np.append(Psol0[0::int(steps/returnedsteps)], Psol0[-1])
    Tsol = np.append(Tsol0[0::int(steps/returnedsteps)], Tsol0[-1])
    Rsol = np.append(Rsol0[0::int(steps/returnedsteps)], Rsol0[-1])

    rhosol = density(Psol, Tsol)
    Esol = eng_gen_pp(rhosol, Tsol)+eng_gen_cno(rhosol, Tsol)
    kappasol = 10**config.interp(np.log10(rhosol), np.log10(Tsol))
    _, Delta_ad, Delta_rad, Delta = Delta_finder(msol, Lsol, Psol, Tsol, Rsol, 
                                                 kappasol, True)
    nature = np.where(Delta == Delta_ad, "radiative", "convective")
    return(msol, rhosol, Lsol, Psol, Tsol, Rsol, Esol, kappasol, Delta_ad, 
           Delta_rad, Delta, nature)