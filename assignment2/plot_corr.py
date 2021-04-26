## plot of correlation values for T = 0.5 and T = 1
## for L = 16


# Programme to read off txt files,
# one temperature
import numpy as np
import matplotlib.pyplot as plt


# take input
Ncycles = 10000;
L =  16;
rank = 0;

# Open file
infile = open("./Results/Data/Corr1D" + str(Ncycles) + \
                "-" + str(L) + "rank" + str(rank) +"-T1.txt", 'r')
infile2 = open("./Results/Data/Corr1D" + str(Ncycles) + \
                "-" + str(L) + "rank" + str(rank) +"-T05.txt", 'r')
# get the temperature
line = infile.readline()
T = line.split()[1]

line = infile2.readline()
T2 = line.split()[1]

# readlines
corr1 = [];
for line in infile:
    numbers = line.split()
    corr1.append(float(numbers[0]))
infile.close()

corr2 = [];
for line in infile2:
    numbers = line.split()
    corr2.append(float(numbers[0]))
infile2.close()

r = np.linspace(0,15,num = 16)

corr1 = np.array(corr1)
corr2 = np.array(corr2)

## plot analytical as well
def corr_func(T,L,r):
    beta = 1./T;
    lamda_p = 2*np.cosh(beta)
    lamda_n = 2*np.sinh(beta)
    c = np.zeros(L)
    factor = 1./((lamda_p**L + lamda_n**L)*L**2)
    for i in range(L):
        c[i] = factor*(lamda_p**(L-r[i])*lamda_n**(r[i]) + lamda_n**(L-r[i])*lamda_p**(r[i]))
    return c


# plot
plt.figure()
plt.title(r"$C(r)$ MC: {:d}".format(Ncycles), fontsize = 14)
plt.plot(r, corr1, "o-", label = "T = 1 [T/J]")
plt.plot(r, corr2, "o-",label = "T = 0.5 [T/J]")
## plot ana label T = 1 ana
plt.plot(r,corr_func(1.0, L,r), label = "ana T = 1 [T/J]" )
plt.plot(r,corr_func(0.5, L,r), label = "ana T = 0.5 [T/J]")
## plt ana label T = 0.5 ana
plt.legend(loc = "upper right")
plt.ylabel(r"$C(r)$")
plt.xlabel("r")
plt.tight_layout()
plt.savefig("./Results/Figs/corrfig.png")
plt.show()
