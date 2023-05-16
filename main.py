import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
from scipy.optimize import root
from astropy.table import Table
from supportfiles.readingtablestools import getTRandK, getlogrho
from supportfiles.integrationandsolutiontools import shootf, solution
import supportfiles.config as config
#reading in Porf. Schlaufman's constants file
from supportfiles.constants import *

logR, logkappa, logT = getTRandK(config.fn, config.tablenum)
#make logT and logR arrays in the shape of the kappa array
logT_k = np.tile(logT[:, np.newaxis], (1, logkappa.shape[1]))
logR_k = np.tile(logR[np.newaxis, :], (logkappa.shape[0], 1))
logrho_k = getlogrho(logR_k, logT_k)
Xgrid, Ygrid = np.meshgrid(np.linspace(-9, 3, 1000), 
                           np.linspace(3.75, 7.5, 1000))
config.interp = LinearNDInterpolator((logrho_k.flatten(), logT_k.flatten()), 
                               logkappa.flatten())

Rt = ((config.MASS / Ms) ** 0.75) * Rs
Lt = ((config.MASS / Ms) ** 3.5) * Ls
rho_const = (3 * config.MASS) / (4 * np.pi * (Rt ** 3))
Pc = (3 * G * (config.MASS ** 2)) / (8 * np.pi * (Rt ** 4))
Tc = (G * config.MASS * config.mu) / (2 * Rt * Na * k)
#constant density model will lowball pressure, so bumping that up a bit
FUDGE_FACTORS = np.array([1, 1.1, 1, 1])
sol = root(shootf, np.array([Lt, Pc, Tc, Rt]) * FUDGE_FACTORS)
print(sol)

msol, rhosol, Lsol, Psol, Tsol, Rsol, Esol, kappasol, Deltaadsol, Deltaradsol, Deltasol, nature = solution(sol.x)
data = Table()
data["M"] = msol
data["rho"] = rhosol
data["R"] = Rsol
data["L"] = Lsol
data["P"] = Psol
data["T"] = Tsol
data["Epsilon"] = Esol
data["kappa"] = kappasol
data["Deltaad"] = Deltaadsol
data["Deltarad"] = Deltaradsol
data["Delta"] = Deltasol
data["nature"] = nature

data.write(f"results/{config.savename}Ms_star.mrt", format="mrt", 
    overwrite=True)
#data.write("results/{}Ms_star.csv".format(config.savename), 
#   overwrite=True) #used this to make my latex table with savedsteps = 50