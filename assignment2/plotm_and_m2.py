# plot magnetization per site
#solving task f) and g) in the prosject description
import numpy as np
import matplotlib.pyplot as plt


Ncycles = 10000         #number of monte carlo cycles
Tnum = 20               #number of temperature points for each core
Tpoints = 4*Tnum        #number of temperature points for the region

L = 16     # Lattice size L

#storing temperatures, mean magnetization
#stored for each value of L
T = np.zeros(int(Tpoints))
M = np.zeros(int(Tpoints))
M2 = np.zeros(int(Tpoints))


#Reading the files produced by the 4 different cores, rank0,rank1,rank2,rank3
infile_rank0 = open("./Results/Data/expvalues" + str(Ncycles) + \
"-" + str(L) + "by" + str(L) + "rank" + str(0) + "-mm2.txt", "r")
infile_rank0.readline()
rank0 = np.loadtxt(infile_rank0)
infile_rank1 = open("./Results/Data/expvalues" + str(Ncycles) + \
"-" + str(L) + "by" + str(L) + "rank" + str(1) + "-mm2.txt", "r")
infile_rank1.readline()
rank1 = np.loadtxt(infile_rank1)
infile_rank2 = open("./Results/Data/expvalues" + str(Ncycles) + \
"-" + str(L) + "by" + str(L) + "rank" + str(2) + "-mm2.txt", "r")
infile_rank2.readline()
rank2 = np.loadtxt(infile_rank2)
infile_rank3 = open("./Results/Data/expvalues" + str(Ncycles) + \
"-" + str(L) + "by" + str(L) + "rank" + str(3) + "-mm2.txt", "r")
infile_rank3.readline()
rank3 = np.loadtxt(infile_rank3)
#rank = np.array((rank0[:,0],rank1[:,0],rank2[:,0],rank3[:,0]))

for a in range(Tnum):
    T[a] = rank0[a,0]
    M[a] = rank0[a,3]         #Filling values for exp values from the first core
    M2[a]= rank0[a,4]

a += 1
for b in range(Tnum):
    T[a+b] = rank1[b,0]
    M[a+b] = rank1[b,3]        #Filling values for exp values from the second core
    M2[a+b] = rank1[b,4]

b += 1
for c in range(Tnum):
    T[a+b+c] = rank2[c,0]
    M[a+b+c] = rank2[c,3]      #Filling values for exp values from the third core
    M2[a+b+c] = rank2[c,4]

c += 1
for d in range(Tnum):
    T[a+b+c+d] = rank3[d,0]
    M[a+b+c+d] = rank3[d,3]      #Filling values for exp values from the fourth core
    M2[a+b+c+d] = rank3[d,4]


plt.figure()
plt.title(r'$<M>$ and $<M^{2}>$ per site')
plt.plot(T, M,  label = r"$<M>$")
plt.plot(T, M2, label = r"$<M^{2}>$")
plt.xlabel('Temperature T')
plt.ylabel('Magnetic moment')
plt.legend(loc = "best",fontsize = 13)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.tight_layout()
plt.savefig("./Results/Figs/m-m2.png")
plt.show()
