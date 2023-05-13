readme.txt this file
main.py the main script. saves a file with the profile of the star.
quickplots.py makes a plot of the pressure, luminosity, temperature, and radius with respect to mass throughout the star. 

supportfiles/
config.py file with the parameters that can be changed by the user. (Stellar mass, opacity table, abundances...)
integrationandsolutiontools.py shooting method function and a solution function that returns the arrays to save as the final output
stellarstructuretools.py functions the describe the interior of a star (density, ODEs, etc.)
readingtablestools.py function that gets T, R, and K from an opacity table (plus a function to convert from R to rho)
constants.py Prof. Schlaufman's astronomical constants file
GN93hz.txt opacities table

plots/
Delta_*Ms.png shows delta_ad and delta_rad over the mass profile of a * solar mass star
internalstructure_*Ms.png shows how L, P, R, and T change over the mass profile of a * solar mass star