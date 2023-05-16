# import mesa_reader to make its classes accessible
import mesa_reader as mr
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import supportfiles.config as config
#reading in Prof. Schlaufman's constants file
from supportfiles.constants import *

path = '/Users/sukay/Desktop/stellar/tutorial/LOGS'
table = Table.read("results/{}Ms_star.mrt".format(config.savename), 
				   format="mrt")

# make a MesaData object from a history file
h = mr.MesaData(path+'/history.data')
L = 10 ** h.log_L

last = mr.MesaLogDir(path) #load the last profile saved (largest model number)
p_last = last.profile_data()
T = 10 ** p_last.logT
R = 10 ** p_last.logR
P = 10 ** p_last.logP

mesaresults = [L[0]*Ls, P[-1], T[-1], R[0]*Rs]
print("Age: {:.2E} years".format(p_last.time_step))

print("MESA \nL        P        T        R")
print("{:.2E} {:.2E} {:.2E} {:.2E}".format(L[0]*Ls, P[-1], T[-1], R[0]*Rs))

print("My Results \nL        P        T        R")
print("{:.2E} {:.2E} {:.2E} {:.2E}".format(table["L"][-1], table["P"][0], 
	  table["T"][0], table["R"][-1]))

myresults = [table["L"][-1], table["P"][0], table["T"][0], table["R"][-1]]

percentagree = (np.subtract(myresults, mesaresults) 
				/ (np.add(myresults, mesaresults) / 2) * 100)
print("Percent Agreement \nL    P    T   R")
print("{:.1f} {:.1f} {:.1f} {:.1f}".format(percentagree[0], percentagree[1], 
	  percentagree[2], percentagree[3]))