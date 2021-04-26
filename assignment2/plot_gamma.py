# programme to plot gamma
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


Ncycles = 50000      #number of monte carlo cycles
Tnum = 100              #number of temperature points for each core
Tpoints = 4*Tnum       #number of temperature points for the region T = [2.2,2.4]

L = np.array([8,16,32])    #different values of lattice size L

#storing temperatures, energy, mean magnetization, heat capacity and Susceptibility
#stored for each value of L
T = np.zeros((4,int(Tpoints)))
gamma = np.zeros((4,int(Tpoints)))

Lstr = ["L = 8", "L = 16", "L = 32"]
#Reading the files produced by the 4 different cores, rank0,rank1,rank2,rank3
for i in range(len(L)):
    infile_rank0 = open("./Results/Data/expvalues" + str(Ncycles) + \
    "-" + str(L[i]) + "by" + str(L[i]) + "rank" + str(0) + ".txt", "r")
    infile_rank0.readline()
    rank0 = np.loadtxt(infile_rank0)
    infile_rank1 = open("./Results/Data/expvalues" + str(Ncycles) + \
    "-" + str(L[i]) + "by" + str(L[i]) + "rank" + str(1) + ".txt", "r")
    infile_rank1.readline()
    rank1 = np.loadtxt(infile_rank1)
    infile_rank2 = open("./Results/Data/expvalues" + str(Ncycles) + \
    "-" + str(L[i]) + "by" + str(L[i]) + "rank" + str(2) + ".txt", "r")
    infile_rank2.readline()
    rank2 = np.loadtxt(infile_rank2)
    infile_rank3 = open("./Results/Data/expvalues" + str(Ncycles) + \
    "-" + str(L[i]) + "by" + str(L[i]) + "rank" + str(3) + ".txt", "r")
    infile_rank3.readline()
    rank3 = np.loadtxt(infile_rank3)
    #rank = np.array((rank0[:,0],rank1[:,0],rank2[:,0],rank3[:,0]))

    for a in range(Tnum):
        T[i,a] = rank0[a,0]
        gamma[i,a] = rank0[a,-1]         #Filling values for exp values from the first core

    a += 1
    for b in range(Tnum):
        T[i,a+b] = rank1[b,0]
        gamma[i,a+b] = rank1[b,-1]        #Filling values for exp values from the second core

    b += 1
    for c in range(Tnum):
        T[i,a+b+c] = rank2[c,0]
        gamma[i,a+b+c] = rank2[c,-1]

    c += 1
    for d in range(Tnum):
        T[i,a+b+c+d] = rank3[d,0]
        gamma[i,a+b+c+d] = rank3[d,-1]      #Filling values for exp values from the fourth core

# plot
for i in range(len(L)):    #plotting the various expectation values for different values of L
    plt.plot(T[i,:],gamma[i,:], label = Lstr[i])

plt.xlabel(r'Temperature $[T/J]$',fontsize=13)
plt.ylabel(r'$\Gamma(T)$',fontsize=13)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.title(r'$\Gamma(T) = <M^{4}>$/${<M^{2}>}^{2}$ (L = 8, 16, 32)',fontsize=14)
plt.axvline(x=2.2691853, color = "r")
plt.text(2.27, 1.0, r"$T_{c} \approx 2.2691853$ ", fontsize=12,color='red')
plt.ylim([0.5,2.0])
plt.xlim([2.20,2.40])
plt.tight_layout()
plt.legend(loc = "best",fontsize = 13)
plt.savefig("./Results/Figs/gamma-TC-L-8-16-32.png")
plt.show()
