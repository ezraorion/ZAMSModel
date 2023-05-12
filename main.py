import numpy as np
import matplotlib.pyplot as plt
from readingtablestools import *
from stellarstructuretools import *
from integrationandsolutiontools import *
import config #global variables go here
from scipy.interpolate import LinearNDInterpolator #CloughTocher2DInterpolator, NearestNDInterpolator
from scipy.optimize import minimize
from astropy.io import ascii
from astropy.table import Table

#reading in Porf. Schlaufman's constants file
exec(open("./constants.py").read())

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
#interpcub = CloughTocher2DInterpolator((logrho_k.flatten(), logT_k.flatten()), logkappa.flatten())
#interpnear = NearestNDInterpolator((logrho_k.flatten(), logT_k.flatten()), logkappa.flatten())
config.interp = interplin #could use this to switch the interpolation routine

#using homology relations (1.87 & 88) for L and R
Rt = ((config.MASS/Ms)**0.75)*Rs
Lt = ((config.MASS/Ms)**3.5)*Ls

#using a constant denisty model (1.41 & 56) as the starting point
rho_const = (3*config.MASS)/(4*np.pi*(Rt**3))
Pc = (3*G*(config.MASS**2))/(8*np.pi*Rt**4)
Tc = (G*config.MASS*config.mu)/(2*Rt*Na*k) 

FUDGEFACTOR = 1.2 #constant density model will defintely lowball pressure and temp
bounds = ((Lt*0.5, Lt*2), (Pc*1e-1,Pc*1e+4),
          (Tc*1e-1, Tc*1e1), (Rt*0.7, Rt*1.5))

res = minimize(shootf, np.array([Lt, Pc*FUDGEFACTOR, Tc*FUDGEFACTOR, Rt]), method='TNC', 
	bounds = bounds, options ={"eta":0.08})
res2 = minimize(shootf, res.x, method='Nelder-Mead', bounds = bounds)
res3 = minimize(shootf, res2.x, method='Nelder-Mead')
res4 = minimize(shootf, res3.x, method='Nelder-Mead')
res5 = minimize(shootf, res4.x, method='Nelder-Mead')

msol, Lsol, Psol, Tsol, Rsol, Esol, kappasol, Deltaadsol, Deltasol, convectivesol = solution(res5.x)
#np.savez("{}solarmass_star.npz".format(config.savename), msol = msol, Lsol = Lsol, Psol = Psol, Tsol = Tsol, Rsol = Rsol)
data = Table()
data['mass'] = msol
data['L'] = Lsol
data['P'] = Psol
data['T'] = Tsol
data['R'] = Rsol
data['E'] = Esol
data['kappa'] = kappasol
data['Deltaad'] = Deltaadsol
data['Delta'] = Deltasol
data['convective'] = convectivesol
ascii.write(data, "{}solarmass_star.ecsv".format(config.savename), overwrite=True)