import numpy as np
from scipy.integrate import solve_ivp
from supportfiles.stellarstructuretools import *
import supportfiles.config as config
#reading in Prof. Schlaufman's constants file
exec(open("./supportfiles/constants.py").read())

def shootf(guess, mc = 1e-12, mf=0.2, ms = 0.9999, steps=1e6):
    """
    inputs: guess (np.array): with
                L (float): luminosity emitted from a sphere enclosing m (ergs s^-1)
                P (float): pressure at the surface of a sphere enclosing m (dyne cm^-2)
                T (float): temperature at the surface of a sphere enclosing m (K)
                R (float): radius of a sphere enclosing m (cm)
            mc (float): central mass point (fractional)
            mf (float): fitting mass point (fractional)
            ms (float): surface mass point (fractional)
            steps (float): steps for integration. split equally between the interior and exterior solution
    returns: score (float): agreement between the interior and exterior solution. should be 0 for a 
                            perfect match, with increasing positive values representing increasingly poor matches.
    """
    L, P, T, R = guess

    #if either solution fails to converge (will have nans), return np.inf so minimize knows this guess was bad
    center_guess = center_bc(mc*config.MASS, P, T)
    masses_cen = np.linspace(mc*config.MASS, mf*config.MASS, num=int(steps/2))
    solc = solve_ivp(derivs, (masses_cen[0], masses_cen[-1]), center_guess, t_eval = masses_cen)
    if np.isnan(np.sum(solc.y[:, -1])):
        return [np.inf, np.inf, np.inf, np.inf]

    surface_guess = surface_bc(config.MASS*ms, L, R)
    masses_surf = np.linspace(config.MASS*ms, config.MASS*mf, num=int(steps/2))
    sols = solve_ivp(fun = derivs, t_span = (masses_surf[0], masses_surf[-1]), y0 = surface_guess, t_eval = masses_surf)
    if np.isnan(np.sum(sols.y[:, -1])):
        return [np.inf, np.inf, np.inf, np.inf]

    residual = (solc.y[:, -1] - sols.y[:, -1])/guess
    return residual

def solution(guess, mc = 1e-12, mf=0.2, ms = 0.9999, steps=1e6):
    """
    inputs: guess (np.array): with
                L (float): luminosity emitted from a sphere enclosing m (ergs s^-1)
                P (float): pressure at the surface of a sphere enclosing m (dyne cm^-2)
                T (float): temperature at the surface of a sphere enclosing m (K)
                R (float): radius of a sphere enclosing m (cm)
            mc (float): central mass point (fractional)
            mf (float): fitting mass point (fractional)
            ms (float): surface mass point (fractional)
            steps (float): steps for integration. split equally between the interior and exterior solution.
    returns: msol (array): log mass solution array (log g)
             Lsol (array): log luminosity solution array (log ergs s^-1)
             Psol (array): log pressure solution array (log dyne cm^-2)
             Tsol (array): log temperature solution array (log K)
             Rsol (array): log radius solution array (log cm)
             Esol, kappasol, Deltaadsol, Deltasol, convectivesol
    """
    L, P, T, R = guess

    center_guess = center_bc(mc*config.MASS, P, T)
    masses_cen = np.linspace(mc*config.MASS, mf*config.MASS, num=int(steps/2))
    solc = solve_ivp(derivs, (masses_cen[0], masses_cen[-1]), center_guess, t_eval = masses_cen)
    
    surface_guess = surface_bc(config.MASS*ms, L, R)
    masses_surf = np.linspace(config.MASS*ms, config.MASS*mf, num=int(steps/2))
    sols = solve_ivp(fun = derivs, t_span = (masses_surf[0], masses_surf[-1]), y0 = surface_guess, t_eval = masses_surf)
    
    msol = np.concatenate((solc.t, np.flip(sols.t)))
    Lsol = np.concatenate((solc.y[0], np.flip(sols.y[0])))
    Psol = np.concatenate((solc.y[1], np.flip(sols.y[1])))
    Tsol = np.concatenate((solc.y[2], np.flip(sols.y[2])))
    Rsol = np.concatenate((solc.y[3], np.flip(sols.y[3])))

    rhosol = density(Psol, Tsol)
    Esol = eng_gen_pp(rhosol, Tsol)+eng_gen_cno(rhosol, Tsol)
    kappasol = 10**config.interp(np.log10(rhosol), np.log10(Tsol))
    _, Delta_ad, Delta_rad, Delta = Delta_finder(msol, Lsol, Psol, Tsol, Rsol, kappasol, True)
    nature = np.where(Delta == Delta_ad, "radiative", "convective")
    return msol, Lsol, Psol, Tsol, Rsol, Esol, kappasol, Delta_ad, Delta_rad, Delta, nature