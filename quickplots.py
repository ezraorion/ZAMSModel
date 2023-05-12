import matplotlib.pyplot as plt
import numpy as np
import config

plt.style.use('dark_background')

npzfile = np.load("{}solarmass_star.npz".format(config.savename))
msol, Lsol, Psol, Tsol, Rsol = npzfile['msol'], npzfile['Lsol'], npzfile['Psol'], npzfile['Tsol'], npzfile['Rsol']

figsize=(10,10)

plt.plot(msol/msol[-1], Lsol/Lsol[-1], label="L/L$_*$")
plt.plot(msol/msol[-1], Psol/Psol[0], label="P/P$_c$")
plt.plot(msol/msol[-1], Tsol/Tsol[0], label="T/T$_c$")
plt.plot(msol/msol[-1], Rsol/Rsol[-1], label="R/R$_*$")

plt.ylim([0,1.125])
plt.legend(ncols=4,loc="upper left")
#axes[0,0].scatter(solc.t[-1], np.log10(solc.y[0, -1]), c='magenta')
#axes[0,0].scatter(sols.t[-1], np.log10(sols.y[0, -1]), c='magenta')
plt.xlabel("M (M$_*$)")

plt.savefig("internalstructure_{}solarmass.png".format(config.savename))