import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from astropy.table import Table
import supportfiles.config as config

matplotlib.use('Agg') #will save figures, but supresses plotting
#plt.style.use("dark_background")

table = ascii.read("{}Ms_star.ecsv".format(config.savename))
#msol, Lsol, Psol, Tsol, Rsol = npzfile['msol'], npzfile['Lsol'], npzfile['Psol'], npzfile['Tsol'], npzfile['Rsol']
#table.pprint()
figsize=(10,10)

plt.plot(table["M"]/table["M"][-1], table["L"]/table["L"][-1], label="L/L$_*$", c = "magenta")
plt.plot(table["M"]/table["M"][-1], table["P"]/table["P"][0], label="P/P$_c$", c = "dodgerblue")
plt.plot(table["M"]/table["M"][-1], table["T"]/table["T"][0], label="T/T$_c$", c = "limegreen")
plt.plot(table["M"]/table["M"][-1], table["R"]/table["R"][-1], label="R/R$_*$", c = "blueviolet")

plt.ylim([0,1.125])
plt.xlim([0,1])
plt.legend(ncols=4,loc="upper left", fontsize=12)
plt.xlabel("m/M$_*$", fontsize=12)
plt.tight_layout()
plt.savefig("./plots/internalstructure_{}Ms.png".format(config.savename))
plt.close()

figsize=(10,10)

plt.plot(table["R"]/table["R"][-1], table["L"]/table["L"][-1], label="L/L$_*$", c = "magenta")
plt.plot(table["R"]/table["R"][-1], table["P"]/table["P"][0], label="P/P$_c$", c = "dodgerblue")
plt.plot(table["R"]/table["R"][-1], table["T"]/table["T"][0], label="T/T$_c$", c = "limegreen")
plt.plot(table["R"]/table["R"][-1], table["M"]/table["M"][-1], label="M/M$_*$", c = "orange")

plt.ylim([0,1.125])
plt.xlim([0,1])
plt.legend(ncols=4,loc="upper left", fontsize=12)
plt.xlabel("r/R$_*$", fontsize=12)
plt.tight_layout()
plt.savefig("./plots/internalstructure_r_{}Ms.png".format(config.savename))
plt.close()

figsize=(10,10)
plt.plot(table["M"]/table["M"][-1], table["Deltaad"], label="$\Delta_{ad}$", c="dodgerblue")
plt.plot(table["M"]/table["M"][-1], table["Deltarad"], label="$\Delta_{rad}$", c = "magenta")
radiative = np.where(table["Deltarad"] < table["Deltaad"], True, False)
convective = np.where(table["Deltarad"] >= table["Deltaad"])
radiativemasses = (table["M"]/table["M"][-1])[radiative]
convectivemasses = (table["M"]/table["M"][-1])[convective]
if radiativemasses[0] != np.min(radiativemasses)\
	or radiativemasses[-1] != np.max(radiativemasses)\
	or convectivemasses[0] != np.min(convectivemasses)\
	or convectivemasses[-1] != np.max(convectivemasses):
		print("Can't do fill between plot because the situation is more commplicated "
			  +"that one convective zone and one radiative zone.")
else:
	ylim = [0, 1]
	plt.fill_betweenx(ylim, radiativemasses[0], radiativemasses[-1], label="Radiative Zone", alpha=0.2, zorder=0, color="dodgerblue")
	plt.fill_betweenx(ylim, convectivemasses[0], convectivemasses[-1], label="Convective Zone", alpha=0.2, zorder=0, color="magenta")
	plt.ylim(ylim)
plt.xlim([0,1])
plt.legend(fontsize=12)
plt.xlabel("m/M$_*$", fontsize=12)
plt.ylabel("$\Delta$", fontsize=12)
plt.tight_layout()
plt.savefig("./plots/Delta_{}Ms.png".format(config.savename))
plt.close()