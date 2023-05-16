import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import supportfiles.config as config

matplotlib.use('Agg') #will save figures, but suppresses plotting
#plt.style.use("dark_background")

table = Table.read("results/{}Ms_star.mrt".format(config.savename), 
				   format="mrt")

#internal structure plot with mass as the dependent variable
figsize=(10,10)
plt.plot(table["M"]/table["M"][-1], table["L"]/table["L"][-1], label="L/L$_*$",
		 c = "magenta", lw = 3)
plt.plot(table["M"]/table["M"][-1], table["P"]/table["P"][0], label="P/P$_c$",
		 c = "dodgerblue", ls = ":", lw = 3)
plt.plot(table["M"]/table["M"][-1], table["T"]/table["T"][0], label="T/T$_c$",
		 c = "limegreen", ls = "--", lw = 3)
plt.plot(table["M"]/table["M"][-1], table["R"]/table["R"][-1], label="R/R$_*$",
		 c = "blueviolet", ls ="-.", lw = 3)
plt.ylim([0,1.125])
plt.xlim([0,1])
plt.legend(ncols=4,loc="upper left", fontsize=12)
plt.xlabel("m/M$_*$", fontsize=12)
plt.tight_layout()
plt.savefig("./plots/internalstructure_{}Ms.png".format(config.savename))
plt.close()

#internal structure plot with radius as the dependent variable
figsize=(10,10)
plt.plot(table["R"]/table["R"][-1], table["L"]/table["L"][-1], label="L/L$_*$",
		 c = "magenta", lw = 3)
plt.plot(table["R"]/table["R"][-1], table["P"]/table["P"][0], label="P/P$_c$",
		 c = "dodgerblue", ls = ":", lw = 3)
plt.plot(table["R"]/table["R"][-1], table["T"]/table["T"][0], label="T/T$_c$",
		 c = "limegreen", ls = "--", lw = 3)
plt.plot(table["R"]/table["R"][-1], table["M"]/table["M"][-1], label="M/M$_*$",
		 c = "orange", ls =(0, (3, 5, 1, 5, 1, 5)), lw = 3)
plt.ylim([0,1.125])
plt.xlim([0,1])
plt.legend(ncols=4,loc="upper left", fontsize=12)
plt.xlabel("r/R$_*$", fontsize=12)
plt.tight_layout()
plt.savefig("./plots/internalstructure_r_{}Ms.png".format(config.savename))
plt.close()

#adiabatic and radiative temperature gradient plot
figsize=(10,10)
plt.plot(table["M"]/table["M"][-1], table["Deltaad"], label="$\Delta_{ad}$",
		 c="blueviolet", ls="--", lw = 3)
plt.plot(table["M"]/table["M"][-1], table["Deltarad"], label="$\Delta_{rad}$",
		 c = "limegreen", ls="-.", lw = 3)
plt.xlim([0,1])
plt.legend(fontsize=12)
plt.xlabel("m/M$_*$", fontsize=12)
plt.ylabel("$\Delta$", fontsize=12)
plt.tight_layout()
plt.savefig("./plots/Delta_{}Ms.png".format(config.savename))
plt.close()