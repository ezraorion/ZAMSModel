readme.txt this file
main.py the main script. saves a file with the profile of the star.
quickplots.py makes plots. 1. L, P, T, and R with respect to mass throughout the star. 2. L, P, T, and m with respect to radius. 3. delta_ad vs. delta_rad with respect to mass.

supportfiles/
config.py file with the parameters that can be changed by the user. (Stellar mass, opacity table, abundances...)
integrationandsolutiontools.py shooting method function and a solution function that returns the arrays to save as the final output.
stellarstructuretools.py functions that describe the interior of a star. (Density, ODEs, etc.)
readingtablestools.py function that gets T, R, and K from an opacity table and a function to convert from R to rho.
constants.py Prof. Schlaufman's astronomical constants file.
GN93hz.txt Grevesse & Noels opacities (1993) table from https://opalopacity.llnl.gov/existing.html#type1W

plots/
Delta_1_435Ms.png shows delta_ad and delta_rad over the mass profile of a 1_435 solar mass star.
internalstructure_1_435Ms.png shows how L, P, R, and T change over the mass profile of a 1_435 solar mass star.
internalstructure_r_1_435Ms.png shows how L, P, R, and T change over the mass profile of a 1_435 solar mass star.
KWW_SSE_2nd_ed_fig18_8_match.png My verion of Fig. 18.8 in Stellar Structure and Evolution (Kippenhahn, Weigert, & Weiss), verifying my energy generation rates from the PP-chain and CNO cycle.
