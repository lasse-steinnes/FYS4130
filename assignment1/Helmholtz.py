## python programme relevant for 3b)
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

## define a function to be used by
## scipy optimize minimize

def F_hat(Ns):
    """
    Input:
    - Ns: int Array of length 3 (Nx,Ny,Nz)
    """
    gamma = 30.0; alpha = 1.0;
    p1 =  Ns[0]*np.log(alpha*Ns[0]/V_hat)
    p2 =  Ns[1]*np.log(alpha*Ns[1]/V_hat)
    p3 = Ns[2]*np.log(alpha*Ns[2]/V_hat)
    p4 = gamma*(Ns[0]*Ns[1] + Ns[1]*Ns[2]  + Ns[2]*Ns[0])/V_hat
    return p1 + p2 + p3 + p4;

## minimize with constraint N = Nx + Ny + Nz while keeping V constant
N0 = np.array([3,3,4],dtype = int) # initial guess
N_tot = 1000; # number of particles


V = np.array([500,800,1000,2000,3000]);
f_vals = np.zeros((N_tot));
x_vals = np.zeros((N_tot,3));


for j in range(len(V)):
    V_hat = V[j];
    for i in range(1,N_tot):
        ## constraints
        N = i;
        bnds = ((0, N), (0, N),(0,N)) ## Ns must be positive
        cons = ({'type': 'eq', 'fun': lambda Ns: N - (Ns[0] + Ns[1] + Ns[2])})
        result = minimize(F_hat, N0, bounds = bnds, constraints = cons, tol = 1e-6) # method = "CG" let it choose method itself?
        f_vals[i] = result.fun;
        x_vals[i,:] = result.x


    N = np.linspace(1,N_tot,N_tot)
    plt.plot(N,f_vals, label = "V: {:}".format(V_hat))

plt.ylabel("F")
plt.xlabel("N")
plt.show()

#print(result.fun) # we find that nx,ny,nz should be equal to each other (gamma =1)
# at equilibrium
## The result of x tells us what the equilibrium values are for Nx, Ny, Nz
## which strongly depend on gamma

"""
Results for gamma = 30.0
"""
# We can find critical pressure from Gibbs duheim equation (extensive),
# dP = dmu N/v -> P  = mu N/V
N_c = [51,84, 106, 225, 362] # critical N from graph
n_c = np.mean(N_c/V)
print("critical concentration {:.2f} ".format(n_c))


## numerical derivative of F wrt N
mu = np.array(np.gradient(f_vals)[362])
P_c = mu*n_c; ## Because we are looking at lst volume and N = i
print("Critical pressure,",P_c)

## orientation of rods at n_c
V_hat = 3000
## minimize with constraint N = Nx + Ny + Nz while keeping V constant
N0 = np.array([3,3,4],dtype = int) # initial guess
N_tot = 1000; # number of particles
V = np.array([500,800,1000,2000,3000]);

N = 362
bnds = ((0, N), (0, N),(0,N)) ## Ns must be positive
cons = ({'type': 'eq', 'fun': lambda Ns: N - (Ns[0] + Ns[1] + Ns[2])})
result = minimize(F_hat, N0, bounds = bnds, constraints = cons, tol = 1e-6) # method = "CG" let it choose method itself?
x_vals = result.x

print("critial orientation:", x_vals)


### task c
## finding out what happens when P changes
mu = np.array(np.gradient(f_vals))
N = np.linspace(1,N_tot,N_tot)
plt.plot(N,mu*N)
plt.ylabel("G")
plt.xlabel("P")
plt.show()cl
