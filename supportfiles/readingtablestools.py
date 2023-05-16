import numpy as np
import re
import pandas as pd
import supportfiles.config as config

def getTRandK(fn, tablenum):
    '''
    Get temperature, R, and kappa columns from an opacity table.

    inputs:
        fn: (string) filename of the opacities table
        tablenum: (int) the number of the table from the LLNL opacities table
        
    returns:
        logR: (np array, dim [20]) density/T6**3 (log g cm^-3 K^-3)
        kappa: (np array, dim [70, 20]) the opacities (g cm^-2)
        logT: (np array, dim [70]) temperatures (K)
    '''
    with open(fn) as f:
        lines = f.readlines()
    start_loc = re.compile("TABLE") #finding the start of each table
    #tableind will have the start of each table in it, with index = table # - 1
    tableind = []
    for n in range(len(lines)):
        if start_loc.match(lines[n]) != None:
            tableind.append(n)

    start = tableind[tablenum-1]
    end = tableind[tablenum]
    #making a list of the logRs & removing the logT string
    logR = lines[start + 4].split()[1:]
    logT = []
    kappa = []
    for i in range(start + 6, end - 1):
        logT.append(float(lines[i].split()[0]))
        prekappa = np.array(lines[i].split()[1:]).astype(float)
        #9.999 is their bad value, switching it to np.nan
        prekappa = np.where(prekappa == 9.999, np.nan, prekappa)
        #for some Rs and Ts there aren't values, 
        #they are at the end of the columns for Table 73, so adding np.nans
        KAPPA_LEN = 19
        a = np.empty((KAPPA_LEN - len(prekappa)))
        a[:] = np.nan
        kappa.append(np.append(prekappa,a))

    return np.array(logR).astype(float), np.array(kappa), np.array(logT)

def getlogrho(logR, logT):
    """
    Go from logR and logT to log rho.

    input: 
        logR: (float) log density/T6**3 (log g cm^-3 K^-3)
        logT: (float) log temperature (log K)
        
    output:
        logrho:(float) log density (log g cm^-3)
    """
    R = 10 ** logR
    T6 = (10 ** logT)/(10 ** 6)
    return np.log10(R* (T6 ** 3))

