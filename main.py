import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
from scipy.optimize import root
from astropy.io import ascii
from astropy.table import Table
import cdspyreadme
from supportfiles.readingtablestools import *
from supportfiles.stellarstructuretools import *
from supportfiles.integrationandsolutiontools import *
import supportfiles.config as config #global variables go here

logR, logkappa, logT = getTRandK(config.fn, config.tablenum) #get the solar values table, Table #73 

#make logT and logR arrays in the shape of the kappa array
logT_k = np.empty(logkappa.shape)
logR_k = np.empty(logkappa.shape)
for n in range(logkappa.shape[0]):
    logT_k[n,:] = logT[n]
for m in range(logkappa.shape[1]):
    logR_k[:,m] = logR[m]
logrho_k = getlogrho(logR_k, logT_k)

Xgrid, Ygrid = np.meshgrid(np.linspace(-9,3,1000), np.linspace(3.75,7.5,1000))

interplin = LinearNDInterpolator((logrho_k.flatten(), logT_k.flatten()), logkappa.flatten())
config.interp = interplin #could use this to switch the interpolation routine

#using homology relations (1.87 & 88) for L and R
Rt = ((config.MASS/Ms)**0.75)*Rs
Lt = ((config.MASS/Ms)**3.5)*Ls

#using a constant denisty model (1.41 & 56) as the starting point
rho_const = (3*config.MASS)/(4*np.pi*(Rt**3))
Pc = (3*G*(config.MASS**2))/(8*np.pi*(Rt**4))
Tc = (G*config.MASS*config.mu)/(2*Rt*Na*k) 

#FUDGEFACTORS = np.array([0.75, 1e2, 3, 1])
#FUDGEFACTORS = np.array([2.5, 100, 2, 1]) 
FUDGEFACTORS = np.array([1, 1.1, 1, 1]) #constant density model will lowball pressure
sol = root(shootf, np.array([Lt, Pc, Tc, Rt])*FUDGEFACTORS)
print(sol)

msol, Lsol, Psol, Tsol, Rsol, Esol, kappasol, Deltaadsol, Deltaradsol, Deltasol, nature = solution(sol.x)
data = Table()
data["M"] = msol
data["L"] = Lsol
data["P"] = Psol
data["T"] = Tsol
data["R"] = Rsol
data["Epsilon"] = Esol
data["kappa"] = kappasol
data["Deltaad"] = Deltaadsol
data["Deltarad"] = Deltaradsol
data["Delta"] = Deltasol
data["nature"] = nature
#ascii.write(data, "{}Ms_star.ecsv".format(config.savename), overwrite=True)
tablemaker = cdspyreadme.CDSTablesMaker()
tablemaker.addTable(data, name="table1")

# add an other local table (in VOTable) 
#table2 = Table.read("table.vot")
#tablemaker.addTable(table2, name="table2")

tablemaker.writeCDSTables()
tablemaker.makeReadMe()
