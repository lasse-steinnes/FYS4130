# numerical evaluation of standard set
import numpy as np
import matplotlib.pyplot as plt

#global variables
N = 10;
a = 0.1;
b = 2.0;
k = 1.0;

factor = 1.01
T_c = 8/(27*k)*a/b
T = factor*T_c
V = np.linspace(10,25,1000)
P_c = a/(27*b**2)

## defining functions
def Pvdw(T,V):
    """
    Van der Waal gas
    """
    return N*k*T/(V-b*N)-a*N**2/V**2

def kt_inv(T,V):
    return (k*N*T*V)/(V-b*N)**2 - 2*a*N**2/V**2

def cv(T,V):
    const = np.zeros(len(V)) + k
    return 3/2*const

def alpha(T,V):
    return N*k*1/kt_inv(T,V)*1/(V-b*N)

def cp(T,V):
    return 3/2*k + (T*V*kt_inv(T,V)*alpha(T,V)**2)/N

x = (Pvdw(T,V)-P_c)


plt.figure()
ax1 = plt.subplot(421)
plt.title(r'Termodynamic standard set, T = {:.2f} $T_c$'.format(factor))
plt.plot(x,1./alpha(T,V))
plt.ylabel(r'1/$\alpha$')
#ax1.set_xlim(-0.2,0.2)
ax1.set_ylim(-1,1)

ax2 = plt.subplot(423, sharex = ax1)
plt.plot(x, kt_inv(T,V))
plt.ylabel(r'1/$\kappa_{t}$')

ax3 = plt.subplot(427, sharex = ax1)
plt.plot(x,1/cv(T,V))
plt.ylabel(r'1/$c_{v}$')

ax4 = plt.subplot(425, sharex = ax1)
plt.plot(x,1/cp(T,V))
plt.ylabel(r'1/$c_{p}$')
plt.xlabel(r'P(V)-P_{c}')

ax5 = plt.subplot(422)
plt.title('Termodynamic standard set')
plt.plot(x,alpha(T,V))
plt.ylabel(r'$\alpha$')

ax6 = plt.subplot(424, sharex = ax5)
plt.plot(x,1./kt_inv(T,V))
plt.ylabel(r'$\kappa_{t}$')

ax7 = plt.subplot(428, sharex = ax5)
plt.plot(x,cv(T,V))
plt.ylabel('$c_{v}$')

ax8 = plt.subplot(426, sharex = ax5)
plt.plot(x,cp(T,V))
plt.ylabel('$c_{p}$')
plt.xlabel('P(V)-$P_{c}$')
plt.tight_layout()
plt.show()
