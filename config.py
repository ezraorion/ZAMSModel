#reading in Porf. Schlaufman's constants file
exec(open("./constants.py").read())

MASS = Ms*2
savename = "2" #make this particular to this stars mass (i.e. 2 for a 2.00 solar mass star)
#Grevesse & Noels (1993) table from https://opalopacity.llnl.gov/existing.html#type1W
fn = '/Users/Sukay/Desktop/stellar/GN93hz.txt'
tablenum = 73
X=0.7000 #these are should be the values for the opacity table number you are using!! (these are for 73)
Y=0.2800
Z=0.0200
mu=4/(3+(5*X)) #from 1.55

Delta_ad = 2/5 #Delta_ad = Gamma_2-1 / Gamma_2, where Gamma_2 = 5/3 for an ideal gas
