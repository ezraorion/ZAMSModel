#reading in Prof. Schlaufman's constants file
exec(open("./supportfiles/constants.py").read())

MASS = Ms*1.435
savename = "1_435" #make this particular to this stars mass 
#(i.e. 2 for a 2.00 solar mass star)
fn = "./supportfiles/GN93hz.txt" ##Grevesse & Noels (1993) table from 
#https://opalopacity.llnl.gov/existing.html#type1W
tablenum = 73
#the abundances are should be the values for the opacity table number 
#you are using!! (these are for 73)
X=0.7000 
Y=0.2800
Z=0.0200
mu=4/(3+(5*X)) #from 1.55

#Delta_ad = Gamma_2-1 / Gamma_2, where Gamma_2 = 5/3 for an ideal gas
Delta_ad = 2/5 
